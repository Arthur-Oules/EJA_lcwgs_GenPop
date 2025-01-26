order_descencing <- function(df) {
  chromosomes_desc_length <- df |>
    count(CHROM) |>
    arrange(desc(n)) |>
    pull(CHROM)
  
  df |>
    mutate(CHROM = factor(CHROM, levels = chromosomes_desc_length)) |>
    arrange(CHROM)
}

order_map <- function(df, map) {
  df |> mutate(CHROM = factor(CHROM, levels = map)) |> arrange(CHROM)
}

Get_window <- function(chromosome, position) {
  chromosome_length <- chromosomes_length |>
    filter(accession == chromosome) |>
    pull(length)
  if (chromosome_length - position < 500) {
    (chromosome_length - 1000):chromosome_length
  } else if (position < 500) {
     0:1000
  } else {
    (position - 500):(position + 500)
  }
}

library("xml2")

XML_to_df <- function(xml_file) {
  xml_tib <- xml_file |>
    as_list() |>
    as_tibble() |> 
    mutate(BlastXML2 = map(BlastXML2, ~if(is.list(.x)) .x else list(.x))) |> 
    unnest_longer(BlastXML2) |>
    filter(BlastXML2_id == "report") |>
    unnest_wider(BlastXML2) |> 
    unnest_longer(Report) |>
    filter(Report_id == "results") |> 
    unnest_wider(Report) |> 
    unnest_longer(Results) |>
    unnest_wider(Results) |>
    unnest_longer(Search) |>
    filter((Search_id == "query-title") | (Search_id == "hits")) |>
    select(c(Search, Search_id)) |>
    mutate(
      Rank = ifelse(
        Search_id == "hits",
        rep(Search[Search_id == "query-title"], each = 2),
        NA
      )
    ) |>
    filter(Search_id == "hits") |>
    select(-Search_id) |>
    unnest(Rank) |> unnest(c(Search, Rank)) |>
    unnest_longer(Search, names_repair = "universal") |>
    filter(Search_id != "len") |> 
    mutate(
      Hit = rep(Search[Search_id == "num"], each = 3)
    ) |> 
    filter(Search_id != "num") |> 
    unnest(Hit) |> unnest(Hit) |>
    unnest_longer(Search, names_repair = "universal")
  
  metadata <- xml_tib |>
    filter(Search_id...2 == "HitDescr") |>
    unnest_wider(Search) |>
    select(c(accession_number = accession, title, sciname, Rank, Hit)) |>
    unnest(c(accession_number, title, sciname))
  
  xml_tib <- xml_tib |>
    filter(Search_id...2 == "Hsp") |> 
    unnest_wider(Search) |>
    select(num, bit_score = "bit-score", evalue, Rank, Hit) |>
    unnest(cols = c(num, bit_score, evalue)) |> 
    left_join(
      x = _,
      y = metadata,
      by = join_by(Rank == Rank, Hit == Hit),
      relationship = "many-to-many"
    ) |> 
    select(
      Rank,
      Hit,
      Hits = num,
      accession_number,
      title,
      sciname,
      bit_score,
      evalue
    ) |>
    mutate_at(vars(Hit, Hits, bit_score, evalue), as.numeric) |> 
    unnest(cols = c(accession_number, title, sciname))
  
  xml_tib
}

library("rentrez")

Get_protein_sequences <- function(df) {
  df |>
    select(accession_number) |>
    mutate(
      Sequences = accession_number |>
        map(
          \(x) entrez_fetch(
            db      = "nuccore",
            id      = x,
            rettype = "fasta_cds_aa",
            retmode = "text"
          )
        ) |> 
        map(
          \(str) {
            str |> 
              strsplit("]") |>
              unlist() |>
              tail(1) |>
              gsub(pattern = "\n", replacement = "",  x = _)
          }
        )
    ) |>
    unnest(Sequences) |>
    filter(!Sequences == "")
}

Get_pvalues <- function(pcadapt) {
  -as.numeric( # Gets non adjusted -log10(pvalues) without infinite values
    pchisq(
      pcadapt$chi2.stat,
      df         = attr(pcadapt, "K"),
      lower.tail = FALSE,
      log.p      = TRUE
    )/log(10)
  )
}

Get_outliers <- function(pcadapt,
                         method           = "bonferroni",
                         pvalue_threshold = .1,
                         minimize_length  = TRUE) {
  
  if (method == "bonferroni") {
    outliers_positions <- readRDS(file = here("data", "SNP_positions.rds")) |> 
      mutate(
        logpvalues = Get_pvalues(pcadapt),
        padj       = p.adjust(pcadapt$pvalues, method = "bonferroni") # Extract adjusted pvalues
      ) |>
      filter(padj < pvalue_threshold)  # Get outliers from pvalues
    
  } else if (method == "manual") {
    outliers_positions <- readRDS(file = here("data", "SNP_positions.rds")) |> 
      mutate( 
        pvalues    = pcadapt$pvalues,
        logpvalues = Get_pvalues(pcadapt)
      ) |> # Extract pvalues
      filter(is.na(pvalues) == FALSE) |> # Remove NAs
      filter(logpvalues >= pvalue_threshold)
    
  } else {
    stop("Method not recognized.")
  }
  
  print(paste0("There are ", as.character(dim(outliers_positions)[1]), " outliers SNPs."))
  
  if (isTRUE(minimize_length)&(dim(outliers_positions)[1]) > 10000) {
    outliers_positions <- outliers_positions |>
      arrange(desc(logpvalues), CHROM, POS) |> 
      filter(row_number() <= 10000) |> 
      arrange(CHROM, POS)
    
    print(paste0("Outliers number was reduced to 10.000 lowest pvalues."))
  }
  
  outliers_positions
}

Get_annotations <- function(positions, chrom_info, gff_annotation) {
  
  positions_annotations <- positions |>
    rename(`GenBank seq accession` = CHROM) |> 
    left_join(chrom_info, by = join_by(`GenBank seq accession`)) |> 
    rowwise() |> 
    mutate(
      annotation = (\(x, y) {
        gff_annotation |>
          filter(seqid == x) |> 
          filter(start <= y & y <= end) |>
          pull(attributes)
      })(`RefSeq seq accession`, POS) |>
        list()
    ) |> 
    filter(length(annotation) != 0)
  
  print(paste0(as.character(dim(positions_annotations)[1]), " outliers were annotated."))
  
  positions_annotations
}

Clean_annotations <- function(annotations) {
  print(paste("There are", as.character(dim(annotations)[1]), "annotations to clean."))
  
  annotations_genes <- annotations |>
    unnest(annotation) |>
    filter(grepl("ID=gene", annotation)) |>
    mutate(
      annotation_ID          = annotation |>
        str_split_i(pattern = ";", i = 1) |>
        str_remove(pattern = "ID="),
      annotation_description = annotation |>
        str_split_i(pattern = ";", i = 4) |>
        str_remove(pattern = "description=")
    ) |> 
    select(-c(annotation)) |>
    distinct(annotation_ID, .keep_all = TRUE)
  
  print(paste("There are", as.character(dim(annotations_genes)[1]), "genes."))
  
  annotations_rnas <- annotations |>
    unnest(annotation) |>
    filter(grepl("ID=rna", annotation)) |>
    mutate(
      annotation_ID      = annotation |>
        str_split_i(pattern = ";", i = 1) |>
        str_remove(pattern = "ID="),
      annotation_product = annotation |>
        str_match("product=\\s*(.*?)\\s*;") |>
        _[, 2]
    ) |>
    select(-c(annotation)) |>
    distinct(annotation_ID, .keep_all = TRUE)
  
  print(paste("There are", as.character(dim(annotations_rnas)[1]), "rnas."))
  
  annotations_exons <- annotations |>
    unnest(annotation) |>
    filter(grepl("ID=exon", annotation)) |>
    mutate(
      annotation_ID      = annotation |>
        str_split_i(pattern = ";", i = 1) |>
        str_remove(pattern = "ID="),
      annotation_gene    = annotation |>
        str_match("gene=\\s*(.*?)\\s*;") |>
        _[, 2],
      annotation_product = annotation |>
        str_match("product=\\s*(.*?)\\s*;") |>
        _[, 2]
    ) |>
    select(-c(annotation)) |>
    distinct(annotation_ID, .keep_all = TRUE)
  
  print(paste("There are", as.character(dim(annotations_exons)[1]), "exons."))
  
  annotations_cds <- annotations |>
    unnest(annotation) |>
    filter(grepl("ID=cds", annotation)) |>
    mutate(
      annotation_ID      = annotation |>
        str_split_i(pattern = ";", i = 1) |>
        str_remove(pattern = "ID="),
      annotation_gene    = annotation |>
        str_match("gene=\\s*(.*?)\\s*;") |>
        _[, 2],
      annotation_product = annotation |>
        str_match("product=\\s*(.*?)\\s*;") |>
        _[, 2]
    ) |>
    select(-c(annotation)) |>
    distinct(annotation_ID, .keep_all = TRUE)
  
  print(paste("There are", as.character(dim(annotations_cds)[1]), "cds."))
  
  annotations_clean <- bind_rows(
    annotations_genes,
    annotations_rnas,
    annotations_exons,
    annotations_cds
  ) |> 
    mutate(
      annotation_type        = annotation_ID |> str_split_i(pattern = "-", i = 1),
      annotation_ID          = annotation_ID |> str_split_i(pattern = "-", i = 2),
      annotation_description = annotation_description |>
        str_remove_all(pattern = "%2C"),
      annotation_product     = annotation_product |>
        str_remove_all(pattern = "%2C")
    )
  
  annotations_clean_len <- annotations_clean |> 
    distinct(`GenBank seq accession`, POS) |>
    dim() |>
    _[1]
  
  print(
    paste(
      "There are",
      as.character(annotations_clean |>
                     dim() |>
                     _[1]),
      "annotations with",
      as.character(annotations_clean_len),
      "unique positions.",
      as.character(
        dim(annotations)[1] - annotations_clean_len
      ),
      "were bad annotations."
    )
  )
  
  if ("padj" %in% colnames(annotations_clean)) {
    annotations_clean |>
      select(
        c("GenBank seq accession", "POS",
          "RefSeq seq accession", "annotation_type", "annotation_ID",
          "annotation_description", "annotation_product", "annotation_gene",
          "padj", "logpvalues")
      ) |>
      arrange(padj, `GenBank seq accession`, POS)
  } else {
    annotations_clean |>
      select(
        c("GenBank seq accession", "POS",
          "RefSeq seq accession", "annotation_type", "annotation_ID",
          "annotation_description", "annotation_product", "annotation_gene",
          "pvalues", "logpvalues")
      ) |>
      arrange(desc(logpvalues), `GenBank seq accession`, POS)
  }
}

lm_eqn <- function(df, r = manteltest$statistic, pp = manteltest$signif) {
  m <- lm(Dgen ~ Dgeo, df)
  eq <- substitute(
    italic(y) == a + b %.% italic(x) * "," ~ ~italic(R)^2 ~ "=" ~ r2 * "," ~ ~italic(p) ~ "=" ~ pp,
    list(
      a = format(unname(coef(m)[1]), digits = 2),
      b = format(unname(coef(m)[2]), digits = 2),
      r2 = format(summary(m)$r.squared, digits = 3),
      pp = format(pp, digits = 3)
    )
  )
  as.character(as.expression(eq))
}
