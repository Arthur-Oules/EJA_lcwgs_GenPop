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

Get_annotations <- function(positions, chrom_info, gff_annotation) {
  positions |>
    rename(`GenBank seq accession` = CHROM) |> 
    left_join(chrom_info) |> 
    rowwise() |> 
    mutate(
      annotation = (\(x, y) {
        gff_annotation |>
          filter(seqid == x) |> 
          filter(start <= y & y <= end) |>
          pull(attributes)
      })(`RefSeq seq accession`, POS) |>
        list()
    )
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