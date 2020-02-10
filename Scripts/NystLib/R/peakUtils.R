# readNarrowpeak
readNarrowPeak <- function(path){
  # read MACS2 narrowpeak file
  bedDf <- readr::read_tsv(path, col_names = F)
  names(bedDf) <- c("chr", "start", "end", "id", "score", "strand", "fold-change", "pValue", "qValue", "summitPos")
  #bedDf %<>% dplyr::select(-strand)
  bed <- GenomicRanges::GRanges(bedDf)
  return(bed) 
}

getPeakData <- function(peakTable, by = "sample_id", narrowPeak_colname = "peakFile"){
  # read each peak dataset, annotate each peak with sample metadata,
  # return dataframe of all peaks
  allPeaks <- lapply(split(peakTable, peakTable[[by]]), function(x){
    peaks <- readNarrowPeak(x[[narrowPeak_colname]])
    #peaks$sample_id <- x$sample_id
    for (name in names(x)){
      # add all metadata associated with each dataset to each peakset
      peaks@elementMetadata[name] <- unlist(x[name])
    }
    
    return(peaks)
  }) %>% 
    GRangesList(.)
  
  allPeakData <- unlist(allPeaks) %>% 
    data.frame(.)
  
  return(allPeakData)
}

get_bw_score <- function(bwFile, peaks, type = "mean", ...){
  # input: bwFile = path to bigwig file,
  # peaks = GRanges object in which to summarize signal
  # type = any type valid to pass to rtracklayer::summary() of bigwigFile object
  # ... passed to summary
  
  bwFileConnection <- rtracklayer::BigWigFile(bwFile)
  
  bw_signal <- GenomicRanges::summary(bwFileConnection, peaks, type = type, ...)
  
  score <- unlist(bw_signal) %>% 
    data.frame %>% 
    .$score
  
  return(score)
}

###
# Peak overlaps ----------
###

peak_overlaps_by_condition <- function(allPeaksList, by = "Rep", ...){
  # input: allPeaksList = GRanges list split by sample group, contains mcol with replicate id (or other col to split)
  # where by = colname to split on 
  # returns: ChIPpeakAnno overlap object
  ol_list <- lapply(names(allPeaksList), function(group_name){
    peaks <- allPeaksList[[group_name]]
    peaks_by_condition <- peaks %>% 
      split(., f = mcols(.)[[by]])
    
    ol <- ChIPpeakAnno::findOverlapsOfPeaks(peaks_by_condition, ...)
    return(ol)
  })
 names(ol_list) <- names(allPeaksList) 
 return(ol_list)
}

plot_overlaps_by_condition <- function(peakOverlaps){
  # input: list of chippeakanno overlap objects from 'peak_overlaps_by_condition()'
  # where names of list objects are the group names (will become the plot titles)
  # output: venn diagrams of all conditions in the list.
  # function is designed for EDA, not publication
  lapply(names(peakOverlaps), function(name){
    ol <- peakOverlaps[[name]]
    ChIPpeakAnno::makeVennDiagram(ol, main = name)
  })
}

annotate_overlap_status <- function(overlappingPeaks, by, anno_colName = "peak_status"){
  # input: overlappingPeaks object from ChiPPeakAnno
  # by = original column name by which peaks were divided
  # output: original peak list with annotated column whether peak is unique or shared by a group
  
  if (class(overlappingPeaks) != "overlappingPeaks"){ stop("Error: overlappingPeaks must be of class 'overlappingPeaks'")} 
   
  mcols(overlappingPeaks$peaksInMergedPeaks)[[anno_colName]] <- "shared" 
  mcols(overlappingPeaks$uniquePeaks)[[anno_colName]] <- paste0(mcols(overlappingPeaks$uniquePeaks)[[by]], "_unique")
  annotated_peaks <- GRangesList(overlappingPeaks$peaksInMergedPeaks,
                                      overlappingPeaks$uniquePeaks) %>% 
    unlist() 
  names(annotated_peaks) <- NULL
  annotated_peaks_List <- annotated_peaks %>% 
    split(., f = mcols(.)[[by]])
  annotated_peaks_List_sorted <- lapply(annotated_peaks_List, sort)
  return(annotated_peaks_List_sorted)
}

###
# Computing Summits from data ----------
###

#resummit <- function(peaks, bam, d, summit_pos_colname = "summitPos"){
#  # call new peak summit from 
#  # resummit peaks:
#    # get pileup of reads across peak with FunChIP
#    # find site of max signal
#    # if there are ties, take mean position between each position and round to nearest bp
#  
#  pileup <- FunChIP::pileup_peak(peaks, bamf = bam, d = d)
#  
#  summitPos <- pileup %>%
#    mcols(.) %>% 
#    data.frame %>% 
#    dplyr::mutate(unique_peak_id =  seq(1,length(peaks))) %>% 
#    split(., .$unique_peak_id) %>% 
#    lapply(., function(x){
#      counts <- x$counts %>% unlist
#      
#      summitPos <- which(counts==max(counts)) %>% 
#        mean %>% 
#        round
#      
#      return(summitPos)
#    }) %>% 
#    do.call("rbind", .)
#  
#  mcols(peaks)[[summit_pos_colname]] <- as.numeric(summitPos)
#  return(peaks)
#}

resummit <- function(peaks, bam, d, ...){
  # outputs vector in order of peaks of summitPos
  summits <- peaks %>% 
    FunChIP::pileup_peak(., bamf = bam, d = d, ...) %>% 
    FunChIP::summit_peak(.) %>% 
    data.frame %>% 
    .$summit_spline
  
  return(summits)
}










###
# QC Plotting ----------
###

peakNumBarPlot <- function(peakData){ 
  # barplot of peak numbers
  barCount <- ggplot(peakData) +
    geom_bar(aes(sample_id, fill = sample_id), color = "black") +
    ggtitle("Total Peaks") +
    xlab(NULL) +
    ylab("Number of Peaks") + 
    theme(legend.position = "none")
  
  return(barCount)
} 

peakQvalHist <- function(peakData, ...){ 
  # histogram of qValues
  peakHist <- ggplot(peakData) + 
    #geom_histogram(aes(x = qValue, fill = sample_id), color = "black" binwidth = 1) +
    geom_histogram(aes(x = qValue, fill = sample_id, color = group_id), binwidth = 1) +
    facet_wrap(~sample_id) +
    ylab("Number of Peaks") + 
    theme(legend.position = "none")
  return(peakHist)
}  

peakCumDistQval <- function(peakData, ...){
  # cumulative distribution of qValue
  peakCumDist <- ggplot(peakData, aes(x = qValue)) + 
    stat_ecdf(geom = "step", aes(color = sample_id)) +
    ggtitle("CumDist of peaks by qValue") +
    ylab("Fraction of Peaks")
  return(peakCumDist)
}







# TODO:
# This is too big & doesn't work great. break up into smaller functions.
# End.
peakListQC <- function(peakTable, by = "sample_id") {
  # read each peak dataset
  allPeaks <- lapply(split(peakTable, peakTable[by]), function(x){
    peaks <- readNarrowPeak(x$peakFile)
    #peaks$sample_id <- x$sample_id
    for (name in names(x)){
      # add all metadata associated with each dataset to each peakset
      peaks@elementMetadata[name] <- unlist(x[name])
    }
    
    return(peaks)
  }) %>% 
    GRangesList(.)
  
  allPeakData <- unlist(allPeaks) %>% 
    data.frame(.)
  
  # barplot of peak numbers
  barCount <- ggplot(allPeakData) +
    geom_bar(aes(sample_id, fill = sample_id), color = "black") +
    ggtitle("Total Peaks") +
    xlab(NULL) +
    ylab("Number of Peaks") + 
    theme(legend.position = "none")
  
  # histogram of qValues
  peakHist <- ggplot(allPeakData) + 
    geom_histogram(aes(x = qValue, fill = sample_id), color = "black", binwidth = 1) +
    facet_wrap(~sample_id) +
    ylab("Number of Peaks") + 
    theme(legend.position = "none")
  
  # cumulative distribution of qValue
  peakCumDist <- ggplot(allPeakData, aes(x = qValue)) + 
    stat_ecdf(geom = "step", aes(color = sample_id)) +
    ggtitle("CumDist of peaks by qValue") +
    ylab("Fraction of Peaks")
  
   
  #plot_grid(barCount, peakHist, peakCumDist) 
  
  gridExtra::grid.arrange(barCount)
  gridExtra::grid.arrange(peakHist)
  gridExtra::grid.arrange(peakCumDist)
}


communities_by_distance <- function(granges, output = "granges"){
  # Input GRanges, output GRanges with "community" column,
  # Or return everything with "all"
  
  if (class(granges) != "GRanges") { stop("input must be GRanges object") }
  
  # Check validity of output type
  valid_output = c("granges", "all")
  output <- tolower(output)
  match.arg(output, valid_output, several.ok = F)
  
  # annotate positionid for community assignment later
  mcols(granges)$graph_rank_id <- seq_along(granges) 
  
  # Create unique edgelist matrix (remove duplicate edges)
  nearestHits <- GenomicRanges::distanceToNearest(granges) %>% 
    data.frame %>% 
    dplyr::filter(queryHits != dplyr::lead(subjectHits)) %>% 
    #dplyr::filter(subjectHits > queryHits) %>% 
    dplyr::select(-distance) %>% 
    as.matrix
  
  # Create graph & find communities with depth-first algorithm
  graph <- graph::ftM2graphNEL(nearestHits, edgemode = "undirected")
  communities <- RBGL::connectedComp(graph) %>% 
    tibble::enframe(., name = "community", value = "graph_rank_id") %>% 
    tidyr::unnest() %>% 
    dplyr::mutate(community = as.integer(community),
                  graph_rank_id = as.integer(graph_rank_id))
 
  # Annotate community membership for each range 
  mcols(granges) <- mcols(granges) %>% 
    data.frame %>% 
    dplyr::left_join(., communities, by = "graph_rank_id") %>% 
    dplyr::select(-graph_rank_id)
  
  switch(output, 
         granges = return(granges),
         all = return(list(ranges = granges, nearestHits = nearestHits, graph = graph, communities = communities))
         )
}


getDistanceCommunityRanges <- function(granges){
  # For each community return the full region as granges
  granges %>% 
    data.frame %>% 
    dplyr::filter(!is.na(community)) %>% 
    dplyr::group_by(community) %>% 
    tidyr::nest() %>% 
    dplyr::mutate(seqnames = purrr::map_chr(data, ~{unique(.$seqnames) %>% as.character()})) %>% 
    dplyr::mutate(start = purrr::map_int(data, ~{min(.$start)})) %>% 
    dplyr::mutate(end = purrr::map_int(data, ~{max(.$end)})) %>% 
    dplyr::select(-data) %>% 
    GRanges 
}
