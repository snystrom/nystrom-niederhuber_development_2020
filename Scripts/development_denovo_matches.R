
plot_denovo_matches_development <- function(dreme_out, database = Sys.getenv("TOMTOM_DB"), method = "bits", 
                                            titles = F, 
                                            base_size = 12, 
                                            base_family = "", 
                                            ...){
  # returns cowplot grids of de-novo matches | tomtom bestmatch motif.
  # can be input to cowplot::plot_grid(.$denovo, .$match, ncol = 2) for plotting
  info <- dreme_out$info
  
  denovo_motifPlots <- map(info$seq, ~{
    seq <- .
    motifPlot <- ggseqlogo::ggseqlogo(dreme_out$motifs[[seq]]@mat, method = method, ...) +
      ggseqlogo::theme_logo(base_size = base_size,
                            base_family = base_family) +
      theme_development
    
    if (titles == T){
      motifPlot <- motifPlot + 
        ggtitle(seq) +
        theme(plot.title = element_text(hjust = 0.5, size = 8))
    }
    
    return(motifPlot)
  })
  
  motif_db <- motifStack::importMatrix(database) 
  
  match_motifPlots <- map(info$bestMatch, ~{
    match_id <- .
    
    if (match_id == "None") { return(NULL) } # if match is None then plot empty motif
     
    match_index <- grep(match_id %>% gsub("[-\\(\\)]", ".", .), names(motif_db)) # need to replace "-" and "()" for "." because importMatrix does.
    #TODO: Fix this for more sophisticated solution?
    match_index <- match_index[1] # return 1st match
    
    match_motif <- motif_db[[match_index]]@mat
    
    motifPlot <- ggseqlogo::ggseqlogo(match_motif, method = method, ...) +
      ggseqlogo::theme_logo(base_size = base_size,
                            base_family = base_family) +
      theme_development
    
    if (titles == T){
      motifPlot <- motifPlot + 
        ggtitle(match_id) +
        theme_development +
        theme(plot.title = element_text(hjust = 0.5, size = 8)) 
    }
    
    return(motifPlot)
  })
 
  denovo <- cowplot::plot_grid(plotlist = denovo_motifPlots, ncol = 1) 
  
  match <- cowplot::plot_grid(plotlist = match_motifPlots, ncol = 1) 
 
  plot_stacks <- list(denovo = denovo,
                      match = match) 
  #cowplot::plot_grid(denovo, match, ncol = 2, labels = c("A", "B"))
  return(plot_stacks)
}