# Utilities for plotting genomic signal with ggplot2

viz_getSignal <- function(region, bigwig, labels, binsize = 10, ...){
  # region = GRanges object
  # bigwig = vector of paths to bigwigs
  # labels = vector of names for each bigwig file
  # binsize = binsize of each bigwig interval in which to compute signal
  #
  # output: data.frame coercible to GRanges with columns named "labels" with genomic signal
  #
  # This is to overcome the issue of variable-width regions in bigwig files for compression purposes, 
  # because this makes variable length outputs, which is not good for visualization. 
  # This fixes the issue by smoothing between intervals (set to binsize of bigwigs) so that all input bigwigs
  # are comparable at a position-by-position basis.
  
  region_tile <- GenomicRanges::slidingWindows(region, width = binsize, step = binsize) %>% 
    unlist(.)
  
  signal <- purrr::map(bigwig, ~{rtracklayer::import.bw(., which = region_tile, as = "NumericList")}) 
  
  signal_mean <- purrr::map(signal, mean)
  names(signal_mean) <- labels
  signal_mean

  region_df <- region_tile %>% data.frame
  
  lapply(labels, function(name){
    name <- quo_name(name)
    region_df <<- region_df %>% 
      dplyr::mutate(!!name := signal_mean[[name]])
  }) 
  
 return(region_df)
}

# write function for reshaping this output to ggplot tidy format,
# reorder variable names by factor order
# write plotting function/ theme for genome viz
viz_meltSignal <- function(vizSignal, value.name = "score", variable.name = "sample"){
  vizSignal %>% 
    reshape2::melt(id.vars = c("seqnames", "start", "end", "width", "strand"), value.name = value.name, variable.name = variable.name)
}

viz_plotTracks <- function(vizSignal, color_by = "sample", facet_by = "sample",
                           ylim = c(min(vizSignal$score), max(vizSignal$score)), breaks = round(ylim), ylabels = breaks){
  # vizSignal is output of viz_meltSignal
  # by is column name of each sample ID to color & facet by
  # alternately can pass color_by & facet_by
  
  #color_by <- quo_name(color_by)
  
  vizSignal %>% 
    ggplot(aes(start, score)) +
      geom_line(aes_string(color = color_by)) +
      geom_area(aes_string(fill = color_by)) +
      facet_wrap(facet_by, ncol = 1, strip.position = "left") +
      #facet_grid(grp ~ ., switch = "both") +
      #theme_minimal() +
      theme(legend.position = "left", 
            strip.text.y = element_text(angle = 180, vjust = 0.5, hjust = 1, color = "black", size = 12), 
            strip.background = element_blank(), 
            strip.placement = "outside", 
            plot.background = element_rect(fill = "white", color = "white"),
            panel.background = element_rect(fill = "white", color = "white"), 
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            axis.line = element_line(color = "black"), 
            axis.text.y = element_text(color = "black", size = 10), 
            axis.text.x = element_text(color = "black")) +
      scale_x_continuous(expand = c(0,0)) +
      # y scaling is broken in ggplot2v3.2.1.9000
      #scale_y_continuous(breaks = breaks, limits = ylim, labels = ylabels) +
      ylab(NULL) +
      xlab(NULL)
}

viz_annotation <- function(title, position = "right", ...){
  # add annotation for panel of tracks
  # exploits lack of ggplot guide which plots externally
  # this is super hacky but it works.
  # ... passes to guide_legend()
  guides(color = F, 
         fill = guide_legend(title = title,
         label = F, 
         title.position = position,
         keywidth = unit(0, "lines"), 
         keyheight = unit(0, "lines"), 
         override.aes = list(alpha = 0), ...))
}

# Highlight peak within a vizualization on all tracks
viz_highlight <- function(region, fill = "green", alpha = 1/3, ymin = -Inf, ymax = Inf){
  # region = granges of peak
  # can be multiple peaks, will highlight each one with same settings
  region_df <- data.frame(region)
  #geom_rect(xmin = min(region_df$start), xmax = max(region_df$end), ymin = ymin, ymax = ymax, fill = fill, alpha = alpha)
  annotate("rect", xmin = region_df$start, xmax = region_df$end, ymin = ymin, ymax = ymax, fill = fill, alpha = alpha)
}

# Highlight peak within a vizualization on all tracks
viz_highlight_aes <- function(regions, fill = NULL, alpha = 1/3, ymin = -Inf, ymax = Inf){
  stop("ERROR: This function doesn't work")
  # region = granges of peak
  # can be multiple peaks, will highlight each one with same settings
  regions_df <- data.frame(regions)
  
  #if (!(fill %in% names(regions_df))) { stop("error: fill must be set to a column of regions") }
  
  #geom_rect(xmin = min(region_df$start), xmax = max(region_df$end), ymin = ymin, ymax = ymax, fill = fill, alpha = alpha)
  #annotate("rect", data = regions_df, xmin = regions_df$start, xmax = regions_df$end, ymin = ymin, ymax = ymax, aes_string(fill = fill), alpha = alpha)
  annotate("rect", data = regions_df, aes_string(xmin = "start", xmax = "end", fill = fill), ymin = ymin, ymax = ymax, alpha = alpha)
}

viz_highlight_range <- function(region, fill = "green", alpha = 1/500, ymin = -Inf, ymax = Inf){
  # region = granges of peak
  # will highlight entire region covered inside "region" if it contains multiple ranges
  region_df <- data.frame(region)
  geom_rect(xmin = min(region_df$start), xmax = max(region_df$end), ymin = ymin, ymax = ymax, fill = fill, alpha = alpha)
}

viz_trackRange <- function(region, mod = 1, units = "bp", round_nits = 1){
  # region = granges of entire plotting region
  # mod = modifier for unit
  # units = "string for unit"
  # round_nits = number of decimals to leave when rounding
  #
  # Returns x axis scale that only shows number of bp 
  # in plotting area
    
  trackRange <- (width(region) / mod) %>% 
    round(., round_nits) %>% 
    glue::glue("{range} {units}", range = .)
  
  midpoint <- resize(region, 1, "center") %>% 
    start
  
  scale_x_continuous(breaks = midpoint, labels = trackRange, expand = c(0,0)) 
    #theme(axis.ticks.x = element_blank())
}

