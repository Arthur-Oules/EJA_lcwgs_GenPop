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

read_admixture <- function(path, pops) {
  path |>
    read_table(col_names = paste0("Pop ", seq(1, pops))) |> 
    mutate(
      Individuals = read_csv2(here("data", "Individuals_coordinates.csv")) |>
        pull(Indiv) |>
        factor(
          levels = here("data", "Individuals_plot_order.txt") |> 
            read.table() |>
            unlist() |>
            as.vector() 
        )
    ) |> 
    relocate(Individuals) |> 
    pivot_longer(
      -c(Individuals),
      names_to  = "Populations",
      values_to = "Ancestry proportions"
    )
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
      annotation = map2(`RefSeq seq accession`, POS,
                        \(x, y) {
                          gff_annotation |>
                            filter(seqid == x) |> 
                            filter(start <= y & y <= end) |>
                            pull(attributes)
                        }) |> list()
    ) |> 
    unnest(annotation) |> 
    filter(annotation |> map(\(x) length(unlist(x))) |> unlist() > 0)
  
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

aggregate_by_population <- function(matrix, populations) {
  unique_pops <- unique(populations)
  pop_dist <- matrix(0, nrow = length(unique_pops), ncol = length(unique_pops))
  rownames(pop_dist) <- colnames(pop_dist) <- unique_pops
  
  for (i in seq_along(unique_pops)) {
    for (j in seq_along(unique_pops)) {
      group1 <- which(populations == unique_pops[i])
      group2 <- which(populations == unique_pops[j])
      if (i == j) {
        pop_dist[i, j] <- mean(matrix[group1, group2])/2
      } else {
        pop_dist[i, j] <- mean(matrix[group1, group2])
      }
    }
  }
  
  return(pop_dist)
}

Substract_diagonal <- function(matrix) {
  for (i in 1:dim(matrix)[1]) {
    for (j in 1:dim(matrix)[1]) {
      if (i == j) {
        matrix[i, j] <- matrix[i, j] - matrix[i, i]
      } else {
        matrix[i, j] <- matrix[i, j] - (matrix[i, i] + matrix[j, j])
      }
      
    }
  }
  
  return(matrix)
}

lm_eqn <- function(df, r, pp) {
  m <- lm(Dgen ~ Dgeo, df)
  
  eq <- substitute(
    atop(italic(y) == a + b %.% italic(x) * "," ~ ~italic(R)^2 ~ "=" ~ r2 * "," ~ ~italic(p) ~ "=" ~ pp, italic(R)^2 ~ "=" ~ mr2 * "," ~ ~italic(p) ~ "=" ~ mpp),
    list(
      a   = format(unname(coef(m)[1]), digits = 2),
      b   = format(unname(coef(m)[2]), digits = 2),
      r2  = format(summary(m)$r.squared, digits = 3),
      pp  = format(pp, digits = 3),
      mr2 = format(r, digits  = 3),
      mpp = format(pp, digits = 3)
    )
  )
  as.character(as.expression(eq))
}

convert_to_outflank <- function(vcf_path, popmap, path_save = NULL) {
  if (!file.exists(path_save)) {
    print(paste0("Converting ", vcf_path))
    
    vcfR <- read.vcfR(vcf_path)
    geno <- extract.gt(vcfR) # Character matrix containing the genotypes
    position <- getPOS(vcfR) # Positions in bp
    chromosome <- getCHROM(vcfR) # Chromosome information
    
    G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
    
    G[geno %in% c("0/0", "0|0")] <- 0
    G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
    G[geno %in% c("1/1", "1|1")] <- 2
    G[is.na(G)] <- 9
    
    EJA_outFLANK <- list("position" = position, "chromosome" = chromosome, "G" = G, "pop" = popmap)
    
    saveRDS(EJA_outFLANK, file = path_save)
  }
}