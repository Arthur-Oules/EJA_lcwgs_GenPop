reduce_all_numeric <- function(liste) {
  for (i in 1:length(liste)) {
    liste[i] <- as.numeric(liste[[i]][1])
  }
  liste
}

Get_outliers <- function(vcf_path, pcadapt, filename = NULL) {
  vcf_to_catalog_ranks <- paste0(vcf_path, "/populations.snps.vcf") |> # Vcf file from ustacks populations output is expected
    read.vcfR() |>
    getID() |>
    strsplit(split = ":") |>
    reduce_all_numeric()
  
  #Outliers selection by pcadapt pvalues
  padj <- pcadapt$pvalues |> p.adjust(method = "bonferroni")
  outliers_vcf_ranks <- which(padj < .1) #alpha = .1 chosen following https://bcm-uga.github.io/pcadapt/articles/pcadapt.html recommendations
  
  outliers_catalog_ranks <- vcf_to_catalog_ranks[outliers_vcf_ranks]
  
  outliers_phi_st_filtered <- read.table( # tsv file from ustacks populations output is expected
    file             = paste0(vcf_path, "/populations.phistats.tsv"),
    header           = FALSE,
    sep              = "\t",
    skip             = 9,
    stringsAsFactors = FALSE
  ) |>
    select(c(1, 4)) |>
    rename("catalog_rank" = V1, "phi_st" = V4) |> 
    filter(catalog_rank %in% outliers_catalog_ranks) |>
    filter(phi_st > .1)
  paste0(
    "There are ", dim(outliers_phi_st_filtered)[1], " outliers with a phi_st > 0.1."
  ) |> print()
  
  #Saving outliers list as txt file for stacks2 whitelist analysis
  if (!is.null(filename)) {
    sink(here("output", paste0(filename, "_outliers_catalog_ranks.txt")))
    for (outlier in outliers_catalog_ranks){
      cat(paste0(outlier, "\n"), append = TRUE)
    }
    sink()
    
    paste0(
      "Outliers list saved as ",
      here("output", paste0(filename, "_outliers_catalog_ranks.txt"))
    ) |> print()
  }
  data.frame(
    "vcf_ranks"     = outliers_vcf_ranks,
    "catalog_ranks" = as.numeric(outliers_catalog_ranks)
  )
}

Get_catalog_sequences <- function(catalog_path,
                                  outliers_ranks,
                                  filename = NULL) {
  catalog_sequences <- read.csv(# fasta file from ?output is expected
    catalog_path,
    col.names = c("Sequences"),
    stringsAsFactors = FALSE
  ) |>
    filter(row_number() %% 2 != 0) |>
    mutate("catalog_rank" = row_number()) |>
    select(catalog_rank, Sequences)
  
  outliers_sequences <- catalog_sequences |>
    filter(row_number() %in% outliers_ranks$catalog_ranks)
  
  if (!is.null(filename)) {
    write.csv(
      outliers_sequences,
      file = here("output", paste0(filename, ".csv"))
    )
    
    #Write as .fa file
    sink(here("output", paste0(filename, ".fa")))
    for (i in 1:dim(outliers_sequences)[1]){
      paste0(">", outliers_sequences$catalog_rank[i], "\n") |> cat(append = TRUE)
      paste0(outliers_sequences$Sequences[i], "\n") |> cat(append = TRUE)
    }
    sink()
  }
  outliers_sequences
}

XML_to_df_old <- function(XML) {
  # Initialisation
  query_df <- data.frame(
    matrix(
      NA,
      nrow = length(xml_find_all(XML, ".//Hit")),
      ncol = 7
    )
  )
  colnames(query_df) <- c(
    "query_number",
    "hit_number",
    "accession_number",
    "title",
    "organism",
    "bit_score",
    "eval_number"
  )
  queries <- xml_find_all(XML, ".//BlastOutput2")
  end <- length(queries)
  counter <- 0
  # Main loop
  for (query in queries) { # Reads XML BLAST results as data frame
    counter <- counter + 1
    print(paste0("Extracting query ", counter, "/", end))
    print("Parsing hits")
    for (Hit in xml_find_all(query, ".//Hit")) { # For each locus get all hits
      match_1 <- xml_find_all(Hit, ".//Hsp") |> _[1] # Keep best match
      query_df[which(rowSums(is.na(query_df)) > 0)[1], ] <- data.frame(
        "query_number"     = xml_find_all(query, ".//query-title") |>
          xml_text() |> 
          as.numeric(),
        "hit_number"       = xml_find_all(Hit, ".//num") |>
          xml_text() |>
          _[1] |> 
          as.numeric(),
        "accession_number" = xml_find_all(Hit, ".//accession") |>
          xml_text() |>
          _[1],
        "title"            = xml_find_all(Hit, ".//title") |>
          xml_text() |>
          _[1],
        "organism"         = xml_find_all(Hit, ".//sciname") |>
          xml_text() |>
          _[1],
        "bit_score"        = xml_find_all(match_1, ".//bit-score") |>
          xml_text() |>
          as.numeric(),
        "eval_number"      = xml_find_all(match_1, ".//evalue") |>
          xml_text() |>
          as.numeric()
      )
    }
    print("End of parsing")
  }
  query_df
}

XML_to_df <- function(xml_file){
  xml_tib <- xml_file |> as_list() |> as_tibble() |> 
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
      Rank = ifelse(Search_id == "hits", rep(Search[Search_id == "query-title"], each = 2), NA)
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
    left_join(x = _, y = metadata, by = join_by(Rank == Rank, Hit == Hit), relationship = "many-to-many") |> 
    select(Rank, Hit, Hits = num, accession_number, title, sciname, bit_score, evalue) |>
    mutate_at(vars(Rank, Hit, Hits, bit_score, evalue), as.numeric) |> 
    unnest(cols = c(accession_number, title, sciname))
  xml_tib
}

Get_protein_sequences <- function(df, silent = TRUE) {
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