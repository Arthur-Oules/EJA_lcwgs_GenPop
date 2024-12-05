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

manhattan_plot_custom_2 <- function(genome_mapping,
                                    pvalues,
                                    outliers_cutoff = NULL) {
  
  labs_rank <- genome_mapping |>
    mutate(labs_rank = row_number()) |>
    distinct(CHROM, .keep_all = TRUE) |>
    left_join(genome_mapping |> count(CHROM)) |> 
    mutate(midranks = round(labs_rank + n/2)) |> 
    filter(row_number() <= 25) |> 
    pull(midranks)
  
  df <- genome_mapping |>
    mutate(
      rank         = row_number(),
      lab_position = if_else(row_number() %in% labs_rank, TRUE, FALSE),
      coloured     = CHROM |>
        {\(x) factor(x, levels = unique(x))}() |>
        as.integer() |>
        (`%%`)(2)
    )
  
  p <- df |> ggplot() +
    geom_point(aes(x = rank, y = {{pvalues}}, colour = factor(coloured)))
  
  if (!is.null(outliers_cutoff)) {
    p <- p + geom_point(
      data   = df |> filter({{pvalues}} >= outliers_cutoff),
      aes(x = rank, y = {{pvalues}}),
      colour = "red"
    ) +
      geom_hline(
        yintercept = outliers_cutoff,
        colour     = "red",
        linetype   = "dashed"
      )
  }
    
  p + scale_colour_manual(values = c("black", "grey")) +
    scale_x_continuous(
      breaks = df$rank[df$lab_position],
      labels = df$CHROM[df$lab_position],
      guide  = guide_axis(angle = 45)
    ) +
    labs(x = "chromosome", y = bquote(-log[10]~("p-value"))) +
    theme_classic() +
    theme(legend.position = "none")
}

PCA_plot <- function(pcadapt_output,
                     popmap,
                     axis_one  = 1,
                     axis_two  = 2,
                     x_offsets = NULL,
                     y_offsets = NULL) {
  # Format data frame
  PCA_df <- tibble(
    x      = -pcadapt_output$scores[, axis_one],
    y      = -pcadapt_output$scores[, axis_two],
    labels = popmap$long_names
  ) |> 
    arrange(labels)
  
  # Compute coords averages by populations
  PCA_average <- PCA_df |>
    group_by(labels) |>
    summarise(average_x = mean(x), average_y = mean(y)) |> 
    ungroup() |>
    mutate(
      labels = labels |>
        gsub(pattern = "_", x = _, replacement = " ") |>
        str_to_title()
    )
  
  # Compute axes percentages
  PCA_percentages <- percent(pcadapt_output$singular.values^2)
  
  PCA <- PCA_df |> ggplot() +
    geom_point(
      mapping = aes(x = y, y = x, fill = labels),
      size    = 4,
      shape   = 21
    ) +
    scale_fill_hue() +
    labs(
      x = paste0("-PC", as.character(axis_two), ": ", PCA_percentages[axis_two]),
      y = paste0("-PC", as.character(axis_one)," : ", PCA_percentages[axis_one])
    ) +
    coord_fixed(ratio = 1.2) +
    theme_minimal() +
    theme(legend.position = "none")
  
  PCA
  if (is.null(x_offsets) || is.null(y_offsets) == TRUE) {
    PCA + geom_text_repel(
      PCA_average,
      mapping = aes(x = average_y, y = average_x, label = labels, fontface = "bold"),
      point.size = 10,
      min.segment.length = 0.3
    )
  } else {
    # Add offsets
    PCA_average <- mutate(
      PCA_average,
      average_x = average_x + x_offsets,
      average_y = average_y + y_offsets
    )

    PCA +
    geom_text(
      PCA_average,
      mapping = aes(x = average_y, y = average_x, label = labels)
    ) +
    theme(legend.position = "none")
  }
}

mantel_plot <- function(Matx, Maty, xlab = NULL, ylab = NULL) {
  
  mantel_result <- vegan::mantel(Matx, Maty)
  
  IBD_PCA_tb <- tibble(
    Dgeo = as.numeric(Matx),
    Dgen = as.numeric(Maty)
  )
  
  p <- ggplot(IBD_PCA_tb) +
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
  
  if (!is.null(xlab) & !is.null(ylab)) {
    p <- p + labs(
      x = xlab,
      y = ylab
    )
  }
  
  p
}