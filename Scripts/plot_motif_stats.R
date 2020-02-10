e93sens_cumdist <- function(peaks, all_motifs, motif_name){
  
  dist_data <- peaks %>% 
    mutate(start = start + summitPos,
           width = 1,
           end = start) %>% 
    GRanges %>% 
    distanceToNearest(., GRanges(all_motifs), ignore.strand = T, select = "arbitrary") %>% 
    data.frame %>% 
    dplyr::mutate(dist_from_summit = distance)
  
  peaks %>% 
    tibble::rowid_to_column("row") %>% 
    dplyr::full_join(., dist_data, by = c("row" = "queryHits")) %>% 
    dplyr::filter(!is.na(e93_sensitive_behavior)) %>% 
    ggplot(aes(dist_from_summit)) +
      stat_ecdf(aes(color = e93_sensitive_behavior), pad = F) +
      theme_cowplot() +
      theme(legend.position = c(0.4, 0.5)) +
      scale_color_manual(values = behavior_colors) +
      labs(x = glue::glue("{motif_name} motif distance from summit (bp)"),
           y = "Fraction of Sites",
           color = NULL) +
      guides(color = guide_legend(override.aes = list(size = 1),
                                  keyheight = unit(0.01, "lines"), 
                                  label.theme = element_text(size = 8, family = "Arial"))) +
      theme_development
}

e93sens_motif_count_pie_factor <- function(peak_motif_anno, motif_name){
  peak_motif_anno %>% 
    dplyr::mutate(nHits = map_dbl(motif_qVal, length)) %>% 
    dplyr::filter(!is.na(e93_sensitive_behavior)) %>% 
    ggplot(aes(factor("Peaks"))) +
    geom_bar(aes(fill = factor(nHits),
                 color = factor(nHits)), position = "fill") +
    facet_wrap(~e93_sensitive_behavior) +
    coord_polar(theta = "y", start = 0) +
    labs(x = NULL, 
         y = NULL,
         fill = glue::glue("{motif_name} Matches")) +
    theme_minimal() +
    theme_development +
    theme(axis.text = element_blank(), 
          panel.grid = element_blank()) +
    guides(color = F)
  
}

e93sens_motif_count_pie_binary <- function(peak_motif_anno, motif_name){
  peak_motif_anno %>% 
    dplyr::mutate(nHits = map_dbl(motif_qVal, length)) %>% 
    dplyr::filter(!is.na(e93_sensitive_behavior)) %>% 
    ggplot(aes(factor("Peaks"))) +
    geom_bar(aes(fill = nHits > 0,
                 color = nHits > 0), position = "fill") +
    facet_wrap(~e93_sensitive_behavior) +
    coord_polar(theta = "y", start = 0) +
    labs(x = NULL, 
         y = NULL,
         fill = glue::glue("{motif_name} Matches")) +
    theme_minimal() +
    theme_development +
    theme(axis.text = element_blank(), 
          panel.grid = element_blank()) +
    guides(color = F) +
    scale_fill_manual(values = c("TRUE" = "Firebrick",
                                 "FALSE" = "Black"),
                      labels = c("TRUE" = "Match",
                                 "FALSE" = "No Match")) +
    scale_color_manual(values = c("TRUE" = "Firebrick",
                                 "FALSE" = "Black"),
                      labels = c("TRUE" = "Match",
                                 "FALSE" = "No Match"))
  
}

e93sens_motif_quality_boxplots <- function(peak_motif_anno, motif_name){
  
  motif_df <- peak_motif_anno %>% 
    dplyr::mutate(nHits = map_dbl(motif_qVal, length)) %>% 
    # first drop any site without a hit, because the q-value is Inf anyway
    dplyr::filter(nHits > 0) %>% 
    dplyr::mutate(minQ = map_dbl(motif_qVal, min, na.rm = T))
  
  # fig.height = 3.5, fig.width =2
  best <- motif_df %>% 
    dplyr::filter(!is.na(e93_sensitive_behavior)) %>% 
    ggplot(aes(e93_sensitive_behavior, -log10(minQ))) +
    geom_boxplot(aes(fill = e93_sensitive_behavior)) +
    scale_fill_manual(values = behavior_colors) +
    guides(fill = F) +
    labs(x = NULL,
         y = "-log10(q-value)",
         title = glue::glue("Best {motif_name} Motif Match")) +
    theme_development +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 8, face = "plain"))
  
  # fig.height = 3.5, fig.width = 2
  all <- motif_df %>% 
    tidyr::unnest(motif_qVal) %>% 
    dplyr::filter(!is.na(e93_sensitive_behavior)) %>% 
    ggplot(aes(e93_sensitive_behavior, -log10(motif_qVal))) +
    geom_boxplot(aes(fill = e93_sensitive_behavior)) +
    scale_fill_manual(values = behavior_colors) +
    guides(fill = F) +
    labs(x = NULL,
         y = "-log10(q-value)",
         title = glue::glue("All {motif_name} Motif Matches")) +
    theme_development +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 8, face = "plain"))
  
  return(list("all" = all,
              "best" = best))
}

e93sens_plot_motif_stats <- function(peak_motif_anno, all_motifs, motif_name, ...){
  
  pie <- list("count" = e93sens_motif_count_pie_factor(peak_motif_anno, motif_name),
              "binary" = e93sens_motif_count_pie_binary(peak_motif_anno, motif_name))
  
  quality_boxplots <- e93sens_motif_quality_boxplots(peak_motif_anno, motif_name)
  
  cumdist <- e93sens_cumdist(peak_motif_anno, all_motifs, motif_name)
  #cumdist <- NULL
  
  return(list("pie" = pie,
              "best_boxplot" = quality_boxplots$best,
              "all_boxplot" = quality_boxplots$all,
              "cumdist" = cumdist))
  
}

cow_arrange_motif_stats <- function(motif_plots){
  h <- 3.5
  W <- 3.5 + 2 + 2
  
  top <- cowplot::plot_grid(
    motif_plots$cumdist,
    motif_plots$all_boxplot,
    motif_plots$best_boxplot,
    nrow = 1,
    rel_widths = c(3.5/W, 2/W, 2/W),
    labels = c("A", "B", "C"),
    label_fontfamily = "Arial",
    label_size = 12
  )
  
  H <- h + 3
  
  plot <- cowplot::plot_grid(
    top, 
    motif_plots$pie$binary, 
    ncol = 1, 
    rel_heights = c(h/H, 3/H),
    labels = c("", "D"),
    label_fontfamily = "Arial",
    label_size = 12
  )
  
  return(list("plot" = plot,
              "height" = H,
              "width" = W))
}