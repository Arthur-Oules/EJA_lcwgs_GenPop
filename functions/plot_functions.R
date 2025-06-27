library("scales")

save_open_plot <- function(path, plot, ...) {
  ggsave(
    path,
    plot   = plot,
    ... = ...
  )
  system2(
    "open",
    args = c("-a Preview.app", path),
    wait = FALSE
  )
}

pcadapt_screeplot <- function(pcadapt) {
  tibble(
    components  = 1:length(pcadapt$singular.values),
    proportions = 100*pcadapt$singular.values^2
  ) |> 
    ggplot(mapping = aes(x = components, y = proportions)) +
    geom_col() +
    geom_label(
      mapping = aes(x = components, y = proportions, label = round(proportions, 2)),
      label.r = unit(0, "lines"),
      vjust   = 1.1
    ) +
    labs(x = "Component rank", y = "Proportion of Variance (%)") +
    scale_x_continuous(breaks = length(pcadapt$singular.values)) +
    scale_y_continuous(breaks = seq(0, 100, 5), minor_breaks = seq(0, 100)) +
    theme(
      panel.background   = element_rect(fill = "white", color = "black"),
      panel.grid.major.y = element_line(
        color = "black",
        linetype = "solid",
        linewidth = .12
      ),
      panel.grid.minor.y = element_line(
        color = "black",
        linetype = "dashed",
        linewidth = .08
      ),
      panel.grid.major.x = element_line(linetype = "blank"),
      panel.grid.minor.x = element_line(linetype = "blank")
    )
}

PCA_plot <- function(pcadapt_output,   popmap,
                     axis_one  = 1,    axis_two  = 2,
                     x_offsets = NULL, y_offsets = NULL) {
  # Format data frame
  PCA_df <- tibble(
    x           = axis_one/abs(axis_one)*pcadapt_output$scores[, abs(axis_one)],
    y           = axis_two/abs(axis_two)*pcadapt_output$scores[, abs(axis_two)],
    populations = popmap
  ) |> 
    arrange(populations)
  
  # Compute coords averages by populations
  PCA_average <- PCA_df |>
    group_by(populations) |>
    summarise(average_x = mean(x), average_y = mean(y)) |> 
    ungroup() |>
    mutate(
      labels = populations |>
        gsub(pattern = "_", x = _, replacement = " ") |>
        str_to_title()
    )
  
  # Compute axes percentages
  PCA_percentages <- percent(pcadapt_output$singular.values^2)
  
  PCA <- PCA_df |> ggplot() +
    geom_point(
      mapping = aes(x = x, y = y, fill = populations),
      size    = 4,
      shape   = 21
    ) +
    scale_fill_hue() +
    labs(
      x = paste0("PC", as.character(axis_one), ": ", PCA_percentages[abs(axis_one)]),
      y = paste0("PC", as.character(axis_two)," : ", PCA_percentages[abs(axis_two)])
    ) +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    theme(legend.position = "none")
  
  PCA
  if (is.null(x_offsets) || is.null(y_offsets) == TRUE) {
    PCA + geom_text_repel(
      PCA_average,
      mapping = aes(x = average_x, y = average_y,
                    label = populations, fontface = "bold"),
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
      mapping = aes(x = average_x, y = average_y, label = populations)
    ) +
    theme(legend.position = "none")
  }
}

manhattan_plot_custom_2 <- function(pcadapt,
                                    SNP_positions, chromosome_map,
                                    outliers_positions = NULL,
                                    outliers_match     = NULL,
                                    data.only = FALSE) {
  
  genome_mapping <- SNP_positions |>
    mutate(logpvalues = Get_pvalues(pcadapt)) |> # Gets log10 of pvalues
    order_map(map = chromosome_map$`GenBank seq accession`)
  
  labs_rank <- genome_mapping |> # Get positions of chromosome names in axis 1
    mutate(labs_rank = row_number()) |>
    distinct(CHROM, .keep_all = TRUE) |>
    left_join(genome_mapping |> count(CHROM)) |> 
    mutate(midranks = round(labs_rank + n/2)) |> 
    filter(row_number() <= 25) |> 
    pull(midranks)
  
  df <- genome_mapping |>
    mutate(
      rank         = row_number(),
      lab_position = if_else(row_number() %in% labs_rank, TRUE, FALSE), # Map chromosome names positions
      coloured     = CHROM |>                                           # Map grey/black colours
        {\(x) factor(x, levels = unique(x))}() |>
        as.integer() |>
        (`%%`)(2)
    )
  
  if (data.only == TRUE) return(df |> select(-c(lab_position, coloured)))
  
  p <- df |> ggplot() +
    geom_point(aes(x = rank, y = logpvalues, colour = factor(coloured))) # Baseplot
  
  if (!is.null(outliers_positions)) {
    p <- p + geom_point(
      data   = df |> inner_join(outliers_positions, by = c("CHROM", "POS")),
      aes(x = rank, y = logpvalues),
      colour = "red"
    ) # +
    #   geom_hline(
    #     yintercept = outliers_cutoff,
    #     colour     = "red",
    #     linetype   = "dashed"
    #   )
  }
  
  if (!is.null(outliers_match)) {
    p <- p + geom_point(
      data   = df |> inner_join(outliers_match, by = c("CHROM", "POS")),
      aes(x = rank, y = logpvalues),
      colour = "green"
    )
  }
  
  p + scale_colour_manual(values = c("black", "grey")) +
    scale_x_continuous(
      breaks = df$rank[df$lab_position],
      labels = df$CHROM[df$lab_position],
      guide  = guide_axis(angle = 45)
    ) +
    labs(x = "Genome Position", y = bquote(-log[10]~("p-value"))) +
    theme_classic() +
    theme(legend.position = "none")
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

admixture_plot <- function(data, pop_map = NULL) {
  plot <- data |>
    ggplot(
      aes(
        x    = Individuals,
        y    = `Ancestry proportions`,
        fill = Populations
      )
    ) +
    geom_col(position = "stack", width = 1) +
    scale_fill_viridis(name = "Ancestral\nPopulations", discrete = TRUE, option = "turbo") +
    labs(x = "Individuals") +
    scale_y_continuous(expand = c(0,0)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid  = element_blank()
    )
  
  if (!is.null(pop_map)) {
    axis <- data |>
      arrange(Individuals) |> 
      left_join(pop_map, by = join_by(Individuals == Indiv)) |> 
      mutate(rank = row_number()) |> 
      group_by(popmap) |> 
      mutate(mean_rank = as.integer(mean(rank))) |> 
      ungroup() |> 
      mutate(mean_rank = Individuals[mean_rank]) |> 
      distinct(popmap, mean_rank) |> 
      mutate(popmap = popmap |> gsub("_", " ", x = _) |> str_to_title())
    
    plot <- plot + scale_x_discrete(breaks = axis$mean_rank, labels = axis$popmap)
  }
  
  plot
}

plot_matrix <- function(df) {
  df |>
    ggplot(aes(x = Var1, y = Var2)) + 
    geom_raster(aes(fill = dist/1000000)) +
    geom_text(aes(label = round(dist/1000000, digits = 2))) +
    scale_fill_viridis(name = "Distance in Mb", option = "turbo") +
    labs(x = "", y = "") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev(levels(df$Var2)))
}

plot_matrix_notext <- function(df) {
  df |> 
    ggplot(aes(x = Var1, y = Var2)) + 
    geom_raster(aes(fill = dist/1000000)) +
    scale_fill_viridis(name = "Distance in Mb", option = "turbo") +
    labs(x = "", y = "") +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev(levels(df$Var2)))
}