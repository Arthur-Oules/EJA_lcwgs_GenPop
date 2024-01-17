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
    
    sink(here("output", paste0(filename, ".fa")))
    for (i in 1:dim(outliers_sequences)[1]){
      paste0(">", outliers_sequences$catalog_rank[i], "\n") |> cat(append = TRUE)
      paste0(outliers_sequences$Sequences[i], "\n") |> cat(append = TRUE)
    }
    sink()
  }
  outliers_sequences
}

XML_to_DF <- function(XML) {
  # Initialisation
  end <- length(xml_find_all(XML, ".//BlastOutput2"))
  counter <- 0
  query_df <- data.frame()
  # Main loop
  for (query in xml_find_all(XML, ".//BlastOutput2")) { # Reads XML BLAST results as data frame
    counter <- counter + 1
    # query_id <- 
    print(paste0("Extracting query ", counter, "/", end))
    print("Parsing hits")
    for (Hit in xml_find_all(query, ".//Hit")) { # For each locus get all hits
      match_1 <- xml_find_all(Hit, ".//Hsp") |> _[1] # Keep best match
      new_row <- data.frame(
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
      query_df <- query_df |>
        bind_rows(new_row)
    }
    print("End of parsing")
  }
  query_df
}

Get_protein_sequences <- function(query) {
  # Initialisation
  n_protein <- dim(query)[1]
  
  prot_sequences <- data.frame(
    "accession_number"  = query$accession_number,
    "sequences" = rep(NA, n_protein)
  )
  
  pb <- txtProgressBar(min = 0, max = n_protein, initial = 0)
  for (i in 1:n_protein) {
    sequence_info <- entrez_fetch( # Fetch protein info from GenBank
      "nuccore",
      id      = prot_sequences$accession_number[i],
      rettype = "fasta_cds_aa",
      retmode = "text"
    )
    prot_sequences$sequences[i] <- # Extract protein sequences from fetch result
      sequence_info |>
      strsplit("]") |>
      unlist()|>
      tail(1) |> 
      gsub(pattern = "\n", replacement = "",  x = _)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  prot_sequences
}