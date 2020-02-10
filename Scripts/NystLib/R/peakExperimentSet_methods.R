# Accessor Functions ----------  

sampleInfo.peakExperimentSet <- function(peakExperimentSet){
	# Returns coldata from deseq as data.frame
  colData(peakExperimentSet$DESeqDataSet) %>% data.frame
}

peakInfo.peakExperimentSet <- function(peakExperimentSet){
  peakExperimentSet$peakInfo %>% data.frame
}

'peakInfo<-.peakExperimentSet' <- function(x, value){
  x$peakInfo <- value
  return(x)
}

dds.peakExperimentSet <- function(peakExperimentSet){
	# Returns coldata from deseq as data.frame
  peakExperimentSet$DESeqDataSet
}

setInfo.peakExperimentSet <- function(peakExperimentSet){
  peakExperimentSet$setInfo
}


res.peakExperimentSet <- function(peakExperimentSet){
  peakExperimentSet$allResults 
}

'res<-.peakExperimentSet' <- function(x, value){
  x$allResults <- value
  return(x)
}

meta.peakExperimentSet <- function(peakExperimentSet) {
	mcols(peakExperimentSet$metaFeatures) %>% data.frame
}

'meta<-.peakExperimentSet' <- function(x, value) {
  old_meta <- metadata(x$metaFeatures) %>% data.frame 
  mcols(x$metaFeatures) <- c(old_meta, value)

  return(x)
}

misc.peakExperimentSet <- function(x) {
  return(x$misc)
}

# XXX:
# Broken:
# Doesn't reassign values in expected way.
# to append/set misc use atomic vectors
#'misc<-.peakExperimentSet' <- function(x, value) {
#  if (is.null(value)){ 
#    x$misc <- NULL
#    } else {
#    x$misc <- c(x$misc, value)
#  }
#  return(x)
#}

metaRanges.peakExperimentSet <- function(peakExperimentSet) {
	peakExperimentSet$metaFeatures
}

sampleRanges.peakExperimentSet <- function(peakExperimentSet) {
	peakExperimentSet$sampleFeatures
}

summits.peakExperimentSet <- function(x) {
  x$resultsSummits
}


# DESeq2 Functions ----------  

generate_DESeq2_contrast_strings.peakExperimentSet <- function(data, grp = "grp", orderParams = "Time"){
    # Generate all unique pairwise comparisons for contrasting DESeq2 results
    # sampleSheet = metafile
    # grp = column in sampleSheet by which to stratify contrasts
    # orderParams = columns in sampleSheet which provide ordering to which component is var1 or var2 in the contrast.
    #   ex. Default is to order temporally, where later timepoint is 2nd.
    
    # Sort to ensure later condition(s) is second, grab unique grps
    # be sure to convert to character in case `grp` is a factor, which will cause numeric values during combn
		
    sampleSheet <- sampleInfo(data)
		contrast_strings <- internal_generate_DESeq2_contrast_strings(sampleSheet, 
																																	grp = grp,
																																 	orderParams = orderParams)
		#data$DESeqDataSet@metadata$contrast_strings <- contrast_strings
    data$setInfo$contrast_strings <- contrast_strings
		return(data)
}


generate_DESeq2_contrasts.peakExperimentSet <- function(peakExperimentSet, lfcCutoff, padjCutoff){
  if (length(setInfo(peakExperimentSet)$contrast_strings) == 0) {
    stop("ERROR: no contrast_strings have been generated for this peakExperimentSet. Try running `generate_DESeq2_contrast_strings`")
  }

  dds <- dds(peakExperimentSet)
  contrasts <- setInfo(peakExperimentSet)$contrast_strings
  allResults <- internal_generate_DESeq2_contrasts(dds, contrasts, lfcCutoff, padjCutoff)
  peakExperimentSet$allResults <- allResults
  return(peakExperimentSet)
}

annotate_results_subsets.peakExperimentSet <- function(data){
  # input: data
  # will go through and add 'subset' column to each results metadata which says which peaks overlap the peaks from the input dataset.
  # can filter based on this value for plotting.
   if(!class(data) == "peakExperimentSet"){
     stop(paste0("ERROR: data is not of class peakExperimentSet. Instead is: ", class(data)))
   }
   
   if (length(res(data)) == 0){
     stop("ERROR: allResults slot is empty. Try running `generate_DESeq2_contrasts`")
   }

  data$allResults <- lapply(data$allResults, function(res){
      c1 <- unique(mcols(res)$condition1)
      c2 <- unique(mcols(res)$condition2)
      
      sampleRanges <- data$sampleFeatures[c(c1, c2)] %>% unlist(.)
      
      bool_overlap <- (res %over% sampleRanges)
      
      mcols(res)$subset <- bool_overlap
      
      return(res)
    })
return(data)
}

computeContrastSummits.peakExperimentSet <- function(data, summitColumn = "summitPos", 
                                                     static_condition = "condition1" , 
                                                     static_summit_extend = 100){
    # compute summits within each contrast based on behavior
    # where increasing peaks use summit from condition2 (ie the dataset in which the peak exists), 
    # decreasing use summit from condition1, 
    # and static use midpoint between condition1 or condition2, or center of condition1&condition2 based on static_condition
  
    # static_condition values:
    # condition1 = use static peaks from condition1
    # condition2 = use static peaks from condition2
    # both = use static peaks from mean of condtion1 & condition2 static peak center (after extend & merge)
  
    # Returns: copy of data with new $resultsSummits slot, with mcols()$behavior_type & mcols()$condition,
    # which describe the peak behavior & condition from which the summit was computed
   
    ########## 
    # Local Functions
    ########## 
  
    
    find_static_summit_mean_position <- function(peaks, c1_summits, c2_summits, extend = static_summit_extend){
      # Extend summits by some number of bases, reduce ranges & take center of new merged peaks
       
      if (!is.numeric(static_summit_extend)){
        stop(paste0("ERROR: static_summit_extend must be numeric. Currently: ", class(static_summit_extend)))
      }
      
      c1_static_summits <- subsetOverlaps_and_flag_nonunique(c1_summits, peaks)
      c2_static_summits <- subsetOverlaps_and_flag_nonunique(c2_summits, peaks)
      
      all_static_summits <- c(c1_static_summits, c2_static_summits) %>% 
                            resize(., width = extend, fix = "center") %>% 
                            reduce %>% 
                            resize(., width = 1, fix = "center") %>% 
                            subsetOverlaps_and_flag_nonunique(., peaks)
      
      return(all_static_summits)
    }
    
    return_static_summits <- function(peaks, c1_summits, c2_summits, static_condition, static_summit_extend){
      
      static_summits <- switch(static_condition,
             condition1 = subsetOverlaps_and_flag_nonunique(c1_summits, peaks),
             condition2 = subsetOverlaps_and_flag_nonunique(c2_summits, peaks),
             both = find_static_summit_mean_position(peaks, c1_summits, c2_summits, static_summit_extend))
      
      return(static_summits)
    }
    
    
    ##########
    # Checks
    ##########
    if (!tolower(static_condition) %in% c("condition1", "condition2", "both")) {
      stop(paste0("ERROR: static_condition must be either condition1, condition2, or both. Current: ", static_condition))
    } 
    
    ##########
    # Compute summits
    ##########
     
    data$resultsSummits <- lapply(data$allResults, function(res){
        if (!"behavior_type" %in% names(mcols(res))){
          stop("ERROR: there is no column named 'behavior_type' in allResults metacolumns. Did you annotate behaviors?")
        }
         
        # select only peaks that overlap peaks from the samples in contrast & that have a defined (non-NA) behavior_type
        resSubset <- res[mcols(res)$subset == T & !is.na(mcols(res)$behavior_type)]
         
        condition1 <- unique(mcols(res)$condition1)
        condition2 <- unique(mcols(res)$condition2)
        
        c1_sampleRangesSummits <- data$sampleFeatures[condition1] %>% unlist %>% computeSummits(., summitColumn = summitColumn)
        c2_sampleRangesSummits <- data$sampleFeatures[condition2] %>% unlist %>% computeSummits(., summitColumn = summitColumn)
        
        resSubsetList <- resSubset %>% split(., f = mcols(.)$behavior_type) 
    
         
        sample_behavior_summits <- lapply(names(resSubsetList), function(behavior){
           peaks <- resSubsetList[behavior] %>% unlist(.) 
          
           # Grab the summit from a given dataset depending on behavior
           # Static peaks can either be returned by condition or by center of a merged summit between 2 conditions
           # annotate behavior_type of summit & from which dataset it originated

           summits <- switch(behavior,
                             Increasing = subsetOverlaps_and_flag_nonunique(c2_sampleRangesSummits, peaks),
                             Decreasing = subsetOverlaps_and_flag_nonunique(c1_sampleRangesSummits, peaks),
                             Static = return_static_summits(peaks,
                                                            c1_sampleRangesSummits,
                                                            c2_sampleRangesSummits, 
                                                            static_condition, 
                                                            static_summit_extend))
           
           mcols(summits)$behavior_type <- behavior
           
           mcols(summits)$condition <- switch(behavior,
                                              Increasing = condition2,
                                              Decreasing = condition1,
                                              Static = switch(static_condition,
                                                              condition1 = condition1,
                                                              condition2 = condition2,
                                                              both = "both"))
           return(summits) 
           })
        
    sample_behavior_summits_GR <- sample_behavior_summits %>% 
                                  GRangesList() %>%
                                  unlist
    
    mcols(sample_behavior_summits_GR)$id <- sample_behavior_summits_GR %>% 
                                         paste0("summit_", seq(1, length(.)))
    
    return(sample_behavior_summits_GR)
    })
    
return(data)
}


get_sample_summits_of_affected_peaks.peakExperimentSet <- function(data, affected_bool_column, sample_summits, dropMetadata = F){
  # return the sample-specific summit positions of peaks affected by a given condition (defined by affected_bool_column)
  # will optionally preserve all metadata (default is to preserve), which can be useful if the affected ranges contain additional information the user would
  # like to keep in order to better understand the sample-specific peaks.
  # Reutrns behavior_type of sample
  affected_ranges <- metaRanges(data)[meta(data)[[affected_bool_column]] == T]
 
  affected_ov <- subsetOverlaps_and_flag_nonunique(affected_ranges, sample_summits, dropMetadata = F) %>% 
                  data.frame() %>% 
                  dplyr::select(-unique, -unique.1) %>% 
                  GRanges
  
  sample_ov <- subsetOverlaps_and_flag_nonunique(sample_summits, affected_ranges, dropMetadata = dropMetadata)
  
  metadata_merge <- c(mcols(sample_ov), mcols(affected_ov))
  mcols(sample_ov) <- NULL
  mcols(sample_ov) <- metadata_merge
  return(sample_ov)
}
