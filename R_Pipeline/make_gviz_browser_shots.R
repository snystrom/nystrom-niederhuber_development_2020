###
# Methods
###

gviz_faire_tracks <- function(range, ylim, ...){
  faire_tracks <- map(faire_info %>% 
                        dplyr::mutate(Sample = factor(Sample, levels = unique(Sample))) %>% 
                        split(.$Sample), ~{
                          Gviz::DataTrack(range = .$BigwigZnorm, genome = "dm3", type = "histogram", 
                                          name = .$id, 
                                          group = .$grp,
                                          fill.histogram = .$color,
                                          col.histogram = .$color, 
                                          #background.title = .$color,
                                          background.title = "white",
                                          fontcolor.title = "black",
                                          col.axis = "black",
                                          rotation.title = 0,
                                          #title.width = 1000,
                                          #cex.title = 0.2,
                                          #cex.axis = axis_cex,
                                          ylim = ylim)
                        })
  return(faire_tracks)
}

gviz_chip_tracks <- function(range, ylim, ...){
  chip_tracks <- map(chip_info %>% 
                       dplyr::mutate(Sample = factor(Sample, levels = unique(Sample))) %>% 
                       split(.$Sample), ~{
                         Gviz::DataTrack(range = .$BigwigZnorm, genome = "dm3", type = "histogram", 
                                         name = .$id, 
                                         group = .$grp,
                                         fill.histogram = .$color,
                                         col.histogram = .$color, 
                                         #background.title = .$color,
                                         background.title = "white",
                                         fontcolor.title = "black",
                                         col.axis = "black",
                                         rotation.title = 0,
                                         #title.width = 1000,
                                         #cex.title = 0.2,
                                         ylim = ylim)
                       })
  
  return(chip_tracks)
}

plotGViz_ChipFaire <- function(range, 
                               ext = c(10000,10000), 
                               scale = sum(ext), 
                               faire_ylim = c(-1, 15), 
                               chip_ylim = c(-1, 15),
                               anno_fill = "black",
                               anno_fontcolor = "white",
                               anno_font.cex = 0.7, 
                               cex.axis = 0.3, 
                               ...){
  # range must be granges object w/ "name" column
  faire_tracks <- gviz_faire_tracks(range, faire_ylim)
  
  chip_tracks <- gviz_chip_tracks(range, chip_ylim)
  
  annotation <- Gviz::AnnotationTrack(range = range, fill = anno_fill, col = NA, 
                                      #featureAnnotation = mcols(range)$name, 
                                      featureAnnotation = "id", 
                                      id = mcols(range)$name,
                                      #feature = mcols(range)$name,
                                      showFeatureId = T,
                                      #group = mcols(range)$name,
                                      #groupAnnotation = "id",
                                      #group = mcols(range)$name,
                                      #just.group = "below", 
                                      #fontcolor.group = "black",
                                      fontcolor.feature = anno_fontcolor,
                                      fontface.feature = "bold",
                                      cex.feature = anno_font.cex,
                                      background.title = NA)
  
  # make empty annotation track so there's a space between FAIRE & ChIP data tracks
  spacer <- Gviz::AnnotationTrack(background.title = "white")
  
  highlight_tracks <- Gviz::HighlightTrack(c(faire_tracks, spacer, chip_tracks), 
                                           range = range, 
                                           fill = genome_highlight_color, 
                                           alpha = 1/2, 
                                           col = NA)
  
  axis <- Gviz::GenomeAxisTrack(scale = scale, col = "grey30")
  
  Gviz::plotTracks(c(axis, highlight_tracks, annotation), 
                   extend.right = ext[2], extend.left = ext[1], 
                   group_annotation = "group", 
                   showTitle = F, 
                   sizes = c(0.1, rep(0.2, length(faire_tracks)), 0.1, rep(0.2, length(chip_tracks)), 0.2),
                   cex.axis = cex.axis, ...) 
  
}


label_gviz_faireChip_plot <- function(range, ...){
  # Take range, plot with plotGViz_ChipFaire then use cowplot to add labels.
  # Either add sample labels or sample + assay labels
 
  # Make browser shot 
  gviz_plot <- grid::grid.grabExpr(
    plotGViz_ChipFaire(range = range, ...)
  )
  
  # generate sample labels
  lab_str <- c("", faire_info$id, "",  chip_info$id, " ") 
  lab_heights <- ifelse(lab_str == "", 1/length(lab_str), 2/length(lab_str)) 
  
  sample_lab <- map(lab_str, ~{
    ggdraw() +
      cowplot::draw_text(.x, hjust = 1, fill = NA, x = 1)
  }) %>% 
    cowplot::plot_grid(plotlist = ., ncol = 1, rel_heights = lab_heights) 
 
  # Generate assay labels 
  assay_str <- c("", " ", "FAIRE", " ", "", "ChIP", " ", " ")
  assay_lab <- map(assay_str, ~{
    ggdraw() +
      cowplot::draw_text(.x, hjust = 1, fill = NA, x = 0.8)
  }) %>% 
    cowplot::plot_grid(plotlist = ., ncol = 1, rel_heights = lab_heights) 
  
  labeled_plot <- list("gviz" = gviz_plot,
                       "sample_labels" = sample_lab,
                       "assay_labels" = assay_lab)
  
  return(labeled_plot)
}

label_gviz_faireChip_plot_horiz <- function(range, ...){
  # Take range, plot with plotGViz_ChipFaire then use cowplot to add labels.
  # Either add sample labels or sample + assay labels
 
  # Make browser shot 
  gviz_plot <- grid::grid.grabExpr(
    plotGViz_ChipFaire(range = range, ...)
  )
  
  # generate sample labels
  lab_str <- c("", faire_info$id, "",  chip_info$id, " ") 
  lab_heights <- ifelse(lab_str == "", 1/length(lab_str), 2/length(lab_str)) 
  
  sample_lab <- map(lab_str, ~{
    ggdraw() +
      cowplot::draw_text(.x, hjust = 1, fill = NA, x = 1)
  }) %>% 
    cowplot::plot_grid(plotlist = ., ncol = 1, rel_heights = lab_heights) 
 
  # Generate assay labels 
  assay_str <- c("FAIRE-seq", "ChIP-seq")
  
  spacer_lab <- ggdraw() +
      cowplot::draw_label("", angle = 90)
  faire_lab <- ggdraw() +
      cowplot::draw_label(expression(underline("FAIRE-seq")), angle = 90)
  chip_lab <- ggdraw() +
      cowplot::draw_label(expression(underline("ChIP-seq")), angle = 90)
  
  assay_lab <- cowplot::plot_grid(plotlist = list(spacer_lab, faire_lab, chip_lab, spacer_lab), 
                                  ncol = 1, 
                                  rel_heights = c(0.2, 1, 1, 0.2)) 
  
  labeled_plot <- list("gviz" = gviz_plot,
                       "sample_labels" = sample_lab,
                       "assay_labels" = assay_lab)
  
  return(labeled_plot)
}

label_gviz_chipOnly_plot <- function(range, ...){
  # Take range, plot with plotGViz_ChipFaire then use cowplot to add labels.
  # Either add sample labels or sample + assay labels
 
  # Make browser shot 
  gviz_plot <- grid::grid.grabExpr(
    plotGViz_ChipFaire(range = range, ...)
  )
  
  # generate sample labels
  lab_str <- c("",  chip_info$id, " ") 
  lab_heights <- ifelse(lab_str == "", 1/length(lab_str), 2/length(lab_str)) 
  
  sample_lab <- map(lab_str, ~{
    ggdraw() +
      cowplot::draw_text(.x, hjust = 1, fill = NA, x = 1)
  }) %>% 
    cowplot::plot_grid(plotlist = ., ncol = 1, rel_heights = lab_heights) 
 
  # Generate assay labels 
  assay_str <- c( " ", "", "ChIP", " ", " ")
  assay_lab <- map(assay_str, ~{
    ggdraw() +
      cowplot::draw_text(.x, hjust = 1, fill = NA, x = 0.8)
  }) %>% 
    cowplot::plot_grid(plotlist = ., ncol = 1, rel_heights = lab_heights) 
  
  labeled_plot <- list("gviz" = gviz_plot,
                       "sample_labels" = sample_lab,
                       "assay_labels" = assay_lab)
  
  return(labeled_plot)
}

plot_label_chip_faire <- function(range, width = 11, type = "chip_faire_horiz", ...){
  # width is absolute plot width in inches
  
  #plots <- label_gviz_faireChip_plot(range, ...)
  plots <- switch(type,
                  chip_faire_horiz = label_gviz_faireChip_plot_horiz(range, ...),
                  chip_faire = label_gviz_faireChip_plot(range, ...),
                  chip_only = label_gviz_chipOnly_plot(range, ...))
  
  # arrange plot
  # a&s correspond to width of assay/sample
  a <- switch(type,
              chip_faire = 1,
              chip_faire_horiz = 0.5,
              chip_only = 1)
  s <- 2
  W <- width
  rel_widths <- c(a/W, s/W, (W-(a+s))/W)
  cowplot::plot_grid(plots$assay_labels, plots$sample_labels, plots$gviz,
                     nrow = 1, rel_widths = rel_widths)
}
