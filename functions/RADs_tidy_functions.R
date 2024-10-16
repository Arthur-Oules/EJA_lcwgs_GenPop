Get_outliers <- function(vcf_path,
                         phistats_path = NULL,
                         pcadapt,
                         filename = NULL) {
  vcf_to_catalog_ranks <- vcf_path |> # Vcf file from ustacks populations output is expected
    read.vcfR() |>
    getID() |>
    strsplit(split = ":") |>
    lapply(\(x) x[1] |> as.numeric())
  
  #Outliers selection by pcadapt pvalues
  padj <- pcadapt$pvalues |> p.adjust(method = "bonferroni")
  outliers_vcf_ranks <- which(padj < .1) #alpha = .1 chosen following https://bcm-uga.github.io/pcadapt/articles/pcadapt.html recommendations
  
  outliers_catalog_ranks <- vcf_to_catalog_ranks[outliers_vcf_ranks]
  
  if (!is.null(phistats_path)) {
    outliers_phi_st_filtered <- read.table( # tsv file from ustacks populations output is expected
      file             = phistats_path,
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
  }
  
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

nc_crop <- function(var,
                    lon_min, lon_max,
                    lat_min, lat_max) {
  # Coordinates conversion
  lon_min <- lon_min * 10
  lon_max <- lon_max * 10
  if (lon_min < 0) {
    lon_min <- 3600 + lon_min
  }
  if (lon_max < 0) {
    lon_max <- 3600 + lon_max
  }
  lat_min <- lat_min * 10 + 750
  lat_max <- lat_max * 10 + 750
  
  # Get and crop data
  var_crop <- ncvar_get(var) |>
    _[, , 1] |>
    as.matrix() |>
    _[(lon_min):(lon_max), (lat_min):(lat_max)] |>
    as.vector()
  
  var_crop
}

reduce_density <- function(high_density_vector, factor) {
  high_density_vector |>
    data.frame("value" = _) |>
    mutate(value = ifelse(row_number() %% factor == 1, value, NA)) |> # Reduce vector density
    pull(value)
}