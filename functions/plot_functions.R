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

manhattan_plot_custom <- function(vcf_file,                # vcf file from ustacks population analysis
                                  pcadapt_analysis,        # Output from the pcadapt() function
                                  n_chromosome,            # Number of effective chromosomes
                                  outliers_ranks = NULL) { # Output from the Get_outliers() function
  chromosome_names <- vcf_file |> # Extract chromosome number map from vcf info
    getCHROM() |>
    lapply(FUN = gsub, pattern = "SCAF_", replacement = "") |>
    as.integer()
  
  # Getting x axis coordinates for chromosome id positions
  mid_pos <- rep(NA, n_chromosome) # Initialize midrank vector
  count_pos <- table(chromosome_names) |> as.data.frame() # Number of locus per chromosome
  
  for (i in 1:n_chromosome){
    mid_pos[i] <-
      Position(\(x) x == i, chromosome_names) + count_pos$Freq[i]/2
  }
  
  notNA.idx <- !is.na(pcadapt_analysis$pvalues)
  
  manhattan_df <- tibble(
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

manhattan_plot_custom_2 <- function(df, pvalues, outliers = NULL) {
  p <- df |> ggplot() +
    geom_point(aes(x = rank, y = {{pvalues}}, colour = factor(coloured)))
  
  if (!is.null(outliers)) {
    p <- p + geom_point(data = outliers, aes(x = rank, y = logpvalues), colour = "red")
  }
    
  p + scale_colour_manual(values = c("black", "grey")) +
    scale_x_continuous(
      breaks = df$rank[df$lab_position],
      labels = df$CHROM[df$lab_position]
    ) +
    labs(x = "chromosome", y = "-log10(p-value)") +
    theme_classic() +
    theme(
      axis.text.x     = element_text(angle = 90),
      legend.position = "none"
    )
}

PCA_plot <- function(pcadapt_output,
                     popmap,
                     x_offsets = NULL,
                     y_offsets = NULL) {
  # Format data frame
  PCA_df <- tibble(
    x   = -pcadapt_output$scores[, 1],
    y   = -pcadapt_output$scores[, 2],
    pop = popmap
  )
  
  # Compute coords averages by populations
  PCA_average <- PCA_df |>
    group_by(pop) |>
    summarise(average_x = mean(x), average_y = mean(y)) |>
    mutate(
      pop = pop |>
        gsub(pattern = "_", x = _, replacement = " ") |>
        str_to_title()
    )
  
  # Compute axes percentages
  PCA_percentages <- pcadapt_output$singular.values^2 |> percent()
  
  PCA <- ggplot() +
    geom_point(
      PCA_df,
      mapping = aes(x = y, y = x, fill = pop),
      size    = 3,
      shape   = 21
    ) +
    labs(
      x = paste0("-PC2: ", PCA_percentages[2]),
      y = paste0("-PC1: ", PCA_percentages[1])
    ) +
    coord_fixed(ratio = 1.2) +
    theme_minimal()
  
  if (is.null(x_offsets) | is.null(y_offsets) == TRUE) {
    PCA + geom_text_repel(
      PCA_average,
      mapping = aes(x = average_y, y = average_x, label = pop)
    )
  } else {
    # Add offsets
    PCA_average <- PCA_average |> mutate(average_x = average_x + x_offsets)
    PCA_average <- PCA_average |> mutate(average_y = average_y + y_offsets)
    
    PCA +
    geom_text(
      PCA_average,
      mapping = aes(x = average_y, y = average_x, label = pop)
    ) +
    theme(legend.position = "none")
  }
}

mantel_plot <- function(Matx, Maty) {
  
  mantel_result <- vegan::mantel(Matx, Maty)
  
  IBD_PCA_tb <- tibble(
    Dgeo = as.numeric(Matx),
    Dgen = as.numeric(Maty)
  )
  
  ggplot(IBD_PCA_tb) +
    geom_smooth(aes(x = Dgeo, y = Dgen), method = "lm", se = TRUE) +
    geom_point(aes(x = Dgeo, y = Dgen)) +
    annotate(
      "text",
      label = lm_eqn(IBD_PCA_tb, r = mantel_result$statistic, pp = mantel_result$signif),
      x     = Inf,
      y     = -Inf,
      parse = TRUE,
      hjust = 1.05,
      vjust = 0
    ) +
    theme_classic() +
    theme(aspect.ratio = 1)
}