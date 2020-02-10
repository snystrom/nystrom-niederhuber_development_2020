# Utilities for dealing with unionPeak list (GRanges)
# Idea here is to build a unionPeak list that contains all sample-level information, 
# then merge all peaks in the union set for DEseq. These functions facilitate extracting sample-specific 
# peak information from the merged peak DEseq results for plotting.

# If you do the feature counting before subsetting the peaklist, 
# then you can go back and dynamically resubset from the count matrix. 
# This could provide a better framework for testing which cutoff to use (ie so computation is only done once.)
# EXAMPLE:
# countMatrix_top7500 <- topQPeakNames(unionPeakList, 7500) %>% subsetCountMatrix(countMatrix, .)


computeSummits.GRanges <- function(peaks, summitColumn = "summitPos"){
  # input GRanges object & summit position column, return Granges with summits
  
  if (!(summitColumn %in% names(mcols(peaks)))) {stop(glue::glue("{summitColumn} does not exist as a column in input ranges"))}
  
  summitpos <- mcols(peaks)[[summitColumn]]
  GenomicRanges::start(peaks) <- GenomicRanges::start(peaks) + summitpos
  GenomicRanges::width(peaks) <- 1
  return(peaks)
}

dropMetadata <- function(granges, cols_to_keep = NULL) {
  # utility for quickly dropping metadata before combining GRanges objects with differing columns
  # cols_to_keep is character vector of column names to keep between objects
  
  #mcols(granges) <- NULL
  grCols <- c("seqnames", "start", "end", "width", "strand")
  keep <- c(grCols, cols_to_keep) 
  granges %<>% 
    data.frame %>% 
    select_(., .dots = keep) %>% 
    GRanges()
  
  return(granges)
}

subset_by_condition.GRanges <- function(results, features_metadata, by = "grp"){
  # Select GRanges regions that match condition1 or condition2 in the 'by' column
  colname <- by
  conditions <- c(unique(results$condition1), unique(results$condition2))
  subset <- features_metadata[mcols(features_metadata)[[colname]] %in% conditions]
  return(SubsetByOverlaps(results,subset))
}

topQPeakNames.GRanges <- function(peakList, nPeaks) {
    # subset descending by qvalue
    # XXX:
    # Rewrite this to actually just take the raw values (colName) as input so the function does the sorting itself if
    # the user wants.
    # End.
      peakList[peakList@elementMetadata$qValueSort <= nPeaks]@elementMetadata$id
}

writeGRanges <- function(GR, filename, colNames = F, bed = T){
   
    if (bed == T) {
        outFile <- GR %>% 
          data.frame %>% 
          dplyr::select(seqnames, start, end)
        outFile$start <- outFile$start - 1
        # avoid scientific notation appearing as coordinates
        outFile$start %<>% format(., scientific = F)
        outFile$end %<>% format(., scientific = F)
        write.table(outFile, file=filename, quote=F, sep="\t", row.names=F, col.names=F)  
    } else {
       outFile <- GR %>% data.frame
       write.table(outFile, file=filename, quote=F, sep="\t", row.names=F, col.names=T)  
    }
    
}

get_mcol_of_overlaps <- function(query, subject, colName, naVal = "NA", fun = return, ...){
  # input: 2 granges objects, 
  # output: vector of values along query objects with values from colName of subjects as vector if > 1 value
  # fun = lapply function wrapping function to run on output list (alternately, could pipe output to lapply to unlist)
  
  # TODO: check colName is found in names(subject)
  # TODO: check query & subject are GRanges
  ov <- findOverlaps(query, subject)
  
  ov_df <- ov %>% 
    data.frame %>% 
    split(.$queryHits)
  
  out <- rep(naVal, length(query)) 
    
  
  lapply(ov_df, function(x){
    i <- unique(x$queryHits)
    
    out[i] <<- mcols(subject[x$subjectHits,])[[colName]] %>% list
  })
  
  fun(out)
}
