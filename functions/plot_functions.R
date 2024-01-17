save_open_plot <- function(path, plot, width, height) {
  ggsave(
    path,
    plot   = plot,
    width  = width,
    height = height
  )
  system2(
    "open",
    args = c("-a Preview.app", path),
    wait = FALSE
  )
}

manhattan_plot_custom <- function( # Draws a manhattan plot using ggplot()
  vcf_file, # vcf file from ustacks population analysis
  pcadapt_analysis, # Output from the pcadapt() function
  n_chromosome, # Number of effective chromosomes
  outliers_ranks = NULL # Output from the Get_outliers() function
) {
  chromosome_names <- vcf_file |> # Extract chromosome number map from vcf info
    getFIX() |>
    _[, "CHROM"] |>
    lapply(FUN = gsub, pattern = "SCAF_", replacement = "") |>
    as.integer()
  
  # Getting x axis coordinates for chromosome id positions
  mid_pos <- rep(NA, n_chromosome) # Initialize midrank vector
  count_pos <- table(chromosome_names) |> as.data.frame() # Number of locus per chromosome
  
  for (i in 1:n_chromosome){
    mid_pos[i] <-
      Position(function(x) x == i, chromosome_names) + count_pos$Freq[i]/2
  }
  
  notNA.idx <- !is.na(pcadapt_analysis$pvalues)
  
  manhattan_df <- data.frame(
    x   = which(notNA.idx), # Catalog rank of locus with non-NA pcadapt p-value 
    y   = -as.numeric( # Gets non adjusted -log10(pvalues) without infinite values
      pchisq(
        pcadapt_analysis$chi2.stat[notNA.idx],
        df         = attr(pcadapt_analysis, "K"),
        lower.tail = FALSE,
        log.p      = TRUE
      )/log(10)
    ), 
    chr = chromosome_names%%2 |> _[notNA.idx] # Chromosome colour mapping
  )
  
  manhattan_plot <- ggplot(manhattan_df, aes(x = x, y = y)) +
    geom_point(aes(colour = factor(chr))) +
    guides(colour = FALSE) +
    scale_color_manual(values = c("black", "grey")) +
    scale_x_continuous(breaks = mid_pos, labels = 1:24) +
    labs(x = "chromosome", y = "-log10(p-value)") +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  if (!is.null(outliers_ranks)) {
    outliers_pvalues_rank <- data.frame(
      "vcf_rank"                 = outliers_ranks$vcf_rank,
      "catalog_rank"             = as.character(outliers_ranks$catalog_ranks),
      "outliers_-log10(pvalues)" = -log10(
        pcadapt_analysis$pvalues[outliers_ranks$vcf_ranks]
      ),
      "outliers_pvalue_plot"     =
        manhattan_df$y[manhattan_df$x %in% outliers_ranks$vcf_ranks]
    )
    
    manhattan_plot_labels <- manhattan_plot +
      geom_label_repel(
        data    = outliers_pvalues_rank,
        mapping = aes(x = vcf_rank, y = outliers_pvalue_plot, label = catalog_rank)
      ) #Manhattan plot with outliers catalog ranks
    manhattan_plot_labels
  } else {
    manhattan_plot
  }
}

PCA_plot <- function(pcadapt_output, popmap, x_offsets, y_offsets) {
  #Format data frame
  PCA_df <- data.frame(
    "x"   = -pcadapt_output$scores[, 1],
    "y"   = -pcadapt_output$scores[, 2],
    "pop" = popmap
  )
  
  populations <- popmap |> as.factor() |> levels()
  
  PCA_average <- data.frame("pop" = populations)
  for (i in 1:length(PCA_average$pop)){
    PCA_average$average_x[i] <- PCA_df |> 
      filter(pop == PCA_average$pop[i]) |> 
      pull(x) |>
      mean()
    PCA_average$average_y[i] <- PCA_df |> 
      filter(pop == PCA_average$pop[i]) |> 
      pull(y) |>
      mean()
  }
  
  PCA_average$pop <- populations |>
    str_to_title() |>
    gsub(pattern = "_", x = _, replacement = " ")
  
  # Add offsets
  PCA_average$average_x <- PCA_average$average_x + x_offsets
  
  PCA_average$average_y <- PCA_average$average_y - y_offsets
  
  PCA_percentages <- pcadapt_output$singular.values^2 |> percent()
  
  # ggplot prompt
  PCA <- ggplot() +
    geom_point(
      PCA_df,
      mapping = aes(x = y, y = x, fill = pop),
      size = 3,
      shape = 21
    ) +
    geom_text(
      PCA_average,
      mapping = aes(x = average_y, y = average_x, label = pop)
    ) +
    labs(
      x = paste0("-PC2: ", PCA_percentages[2]),
      y = paste0("-PC1: ", PCA_percentages[1])
    ) +
    coord_fixed(ratio = 1.2) +
    theme(
      legend.position  = "none",
      panel.background = element_rect(
        fill = "white",
        color = "black"
      ),
      panel.grid.major = element_line(
        color = "black",
        linetype = "dashed",
        linewidth = .12
      )
    )
  PCA
}

nc_crop <- function(var, lon_min, lon_max, lat_min, lat_max) {
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