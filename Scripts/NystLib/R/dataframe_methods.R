generate_DESeq2_contrast_strings.data.frame <- function(sampleSheet, grp = "grp", orderParams = "Time"){
  # Generate all unique pairwise comparisons for contrasting DESeq2 results
  # sampleSheet = metafile
  # grp = column in sampleSheet by which to stratify contrasts
  # orderParams = columns in sampleSheet which provide ordering to which component is var1 or var2 in the contrast.
  #   ex. Default is to order temporally, where later timepoint is 2nd.
  
  # Sort to ensure later condition(s) is second, grab unique grps
  # be sure to convert to character in case `grp` is a factor, which will cause numeric values during combn
  
  contrast_strings <- internal_generate_DESeq2_contrast_strings(sampleSheet, 
                                                                grp = grp,
                                                                orderParams = orderParams)
  
  return(contrast_strings)
}


generate_DESeq2_contrasts.DESeqDataSet <- function(dds, contrasts, lfcCutoff, padjCutoff, format = "data.frame"){
  
  # XXX:
  # Check that contrasts is contrast_string_object
  #stop("ERROR: no contrast_strings have been generated for this peakExperimentSet. Try running `generate_DESeq2_contrast_strings`")
  # End.

  allResults <- internal_generate_DESeq2_contrasts(dds, contrasts, lfcCutoff, padjCutoff, format = format)
  
  return(allResults)
}


pairwise_comparisons <- function(sampleSheet, grp = "grp", orderParams = "Time"){
    # Generate all unique pairwise comparisons for samples in sampleSheet at the grp level
    # sampleSheet = metafile
    # grp = column in sampleSheet by which to stratify contrasts
    # orderParams = columns in sampleSheet which provide ordering to which component is var1 or var2 in the contrast.
    #   ex. Default is to order temporally, where later timepoint is 2nd.
    
    # Sort to ensure later condition(s) is second, grab unique grps
    # be sure to convert to character in case `grp` is a factor, which will cause numeric values during combn
    grps <- sampleSheet %>% 
            dplyr::arrange_(orderParams) %>% 
            dplyr::select_(grp) %>%
            .[[grp]] %>% as.character() %>% 
            unique() 
    
    grpCombn <- c(combn(grps, m = 2))
    
    resStr <- list()
    for (i in seq(1, length(grpCombn), by = 2)){
        condition1 <- grpCombn[i]
        condition2 <- grpCombn[i + 1]
        
        comparison <- paste0(condition1, "vs", condition2)
        resStr[[comparison]]$contrast <- c(grp, condition1, condition2)
        resStr[[comparison]]$comparison <- comparison
        resStr[[comparison]]$condition1 <- condition1
        resStr[[comparison]]$condition2 <- condition2
    }
    return(resStr)
}
