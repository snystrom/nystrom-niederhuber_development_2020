internal_generate_DESeq2_contrast_strings <- function(sampleSheet, grp = "grp", orderParams = "Time"){
    # Generate all unique pairwise comparisons for contrasting DESeq2 results
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

    
annotate_DESeq2_behavior <- function(dfResults, lfcCutoff = 1, padjCutoff = 0.05, infer_annotation = T){
    # Input: DataFrame of DEseqResults
    # Output: $behavior column with Higher_Condition1 / Higher_Condition2 / Static / NA values
    # Where NA are all features that don't meet lfcCutoff but are statistically significant (meet padjCutoff)
    # infer_annotation = T will check dfResults for condition1 & condition2 columns and replace their values for
    # Condition1/2 respectively.
    #
    # Will also return $behavior_type which is the interpreted behavior ("Increasing", "Decreasing", "Static", "NA")
    # useful in some cases for plotting, or comparison between contrasts of different conditions.
    # (ex. Do the features that increase in WT-Time1 vs WT-Time2 also increase in Mut-Time1 vs Mut-Time2?)
    
    # USAGE: dds_res <- annotate_DESeq2_behavior(dds_res, lfcCUtoff = 1, padjCutoff = 0.05)
    require(magrittr)
  
    ####################  
    # Configure parameters:
    ####################  
    
    # if condition1 and condition2 columns exist, change string to correct annotation, unless infer_annotation = F 
    if (infer_annotation == T) {
        if (("condition1" %in% names(dfResults)) & ("condition2" %in% names(dfResults))){
            c1 <- unique(dfResults$condition1)
            c2 <- unique(dfResults$condition2)
            
            if (length(c1) == 1 & length(c2) == 1) {
            higher_c1_string = paste0("Higher_", c1)
            higher_c2_string = paste0("Higher_", c2)
            }
        }
    } else {
        higher_c1_string = "Higher_Condition1"
        higher_c2_string = "Higher_Condition2"
    }
    
    #################### 
    # Annotate Behaviors:
    #################### 
  # XXX: WARNING:
    dfResults$behavior <- NA # SETTING THIS TO NA is dangerous, consider "NA"
    # just don't want to break any existing code so postponing this change
    
    # Things more closed in 1st dataset 
    dfResults <- dplyr::mutate(dfResults, 
                               behavior = 
                               replace(behavior, (log2FoldChange < -lfcCutoff & padj < padjCutoff), higher_c2_string))
    
    # Things more open in 1st dataset 
    dfResults <- dplyr::mutate(dfResults,
                               behavior = 
                               replace(behavior, (log2FoldChange > lfcCutoff & padj < padjCutoff), higher_c1_string))
    
    # Things that are not statistically significant are static 
    dfResults <- dplyr::mutate(dfResults,
                               behavior = replace(behavior, (padj > padjCutoff & !is.na(padj)), "Static"))
    
    #################### 
    # Annotate Behavior Type:
    #################### 
    dfResults <- dfResults %>% dplyr::mutate(behavior_type = case_when(
                                             grepl(higher_c1_string, .$behavior) ~ "Decreasing",
                                             grepl(higher_c2_string, .$behavior) ~ "Increasing",
                                             .$behavior == "Static" ~ "Static"))
    
    return(dfResults)
} 


internal_generate_DESeq2_contrasts <- function(dds, contrast_string_object, lfcCutoff = 1, padjCutoff = 0.05,
                                               format = "GRanges"){
  # XXX:
  # Include parameter that allows parallel = F & doesn't require BiocParallel.
  require(BiocParallel)
  # End.
  
    # Input: DEseqDataset & the constrast_string_object from `generate_DESeq_contrast_strings()`
    # Output: list of each contrast w/ annotated behavior, id, comparison, condition1/2 columns
    
   # XXX:
   # ERROR: This fails if contrast_string_object is not a list (contains single object)
   # End.
    resList <- bplapply(contrast_string_object, function(x){
      # TODO:
      # add format option so that GRanges or DataFrame could be used (be sure when implement to fix return statment as well!)
      # "DataFrame" option
        #results <- DESeq2::results(dds, contrast = x$contrast, format = "DataFrame") %>% 
        #           data.frame(.)
      # End.
        if (tolower(format) == "granges"){
          results_GR <- DESeq2::results(dds, contrast = x$contrast, format = "GRanges") 
          results <- mcols(results_GR) %>% data.frame()
          
        } 
        if (tolower(format) == "dataframe"){
          results <- DESeq2::results(dds, contrast = x$contrast, format = "DataFrame") %>% data.frame()
        }
        
        if (!tolower(format) %in% c("granges", "dataframe")) {
          stop("ERROR: format must be either GRanges or dataframe")
        }
        
        # Add annotations describing contrast to results object  
        results$comparison <- x$comparison
        results$condition1 <- x$condition1
        results$condition2 <- x$condition2
        
        # Add feature names (ex. genes or peaks) to 'id' column 
        results$id <- rownames(dds)
        
        results <- annotate_DESeq2_behavior(results, lfcCutoff = lfcCutoff, padjCutoff = padjCutoff, infer_annotation = T)
        # TODO:
        # replace to return to using dataframe, not GRanges if format == "DataFrame"
        #return(results)
        # End.
        
        if (tolower(format) == "granges"){
          mcols(results_GR) <- results
          return(results_GR)
          
        } 
        if (tolower(format) == "dataframe"){
          return(results)
        }
    })
    
    return(resList)
}

compare_behavior_between_contrasts <- function(reference, test, feature_metadata, 
                                               annotation_basename, 
                                               idcol = "id",
                                               affected_prefix = "affected_",
                                               unaffected_prefix = "unaffected_"){
	    # Where reference & test are dataframes generated by generate_DESeq2_contrasts 
	    # Purpose is to compare for example, differences in behavior between two conditions.
	    # for example, peaks that increase in wild-type are annotated in WT-Time1vsWT-Time2, 
	    # while peaks that are affected by mutating a gene are annotated in mut-Time1vsmut-Time2.
	    # By comparing the differences in behavior between these two conditions, we can determine which 
	    # regions are affected by the mutation, ergo that gene is required for the wild-type changes to occurr.
	    #
	    # Usage: `reference` is the "wild-type", `test` is the "mutant" in the above example.
	    #
	    # features are annotated as <affected | unaffected>_<behavior_in_results1>
        # Two columns are added to the feature_metadata file, annotation_basename_full & annotation_basename_simple
        # simple is annnoatted as "(un)affected_<behavior_in_reference>"
        # full is annnoatted as "(un)affected_<behavior_in_reference>.<behavior_in_test>"
        # The idea here is to allow grouping by many parameters (complex & simple) 
        # and to abstract as many of those manipulations away from the user as possible for ease of use and code readability.
    
        # feature_metadata is at minimum a dataframe with a single 'id' column (matching the features used in DESeq)
        #       can be easily generated by feature_metadata <- data.frame("id" = rownames(dds))
        # annotation_basename is the base string for the output columns (<base>_full, <base>_simple)
        # reference & test are two DEseq2 contrasts. 
            # WARNING: Ensure their direcionality makes sense so that the behaviors are comparable!
            # (ie ensure that for both contrasts "increasing" means the same thing in each individual dataset, 
            # otherwise output will be exactly opposite of expected)
    
    ########
	# Checks:
    ########
  
    feature_metadata <- switch(class(feature_metadata),
                               GRanges = mcols(feature_metadata) %>% data.frame,
                               data.frame = feature_metadata)
    
    if (!("id" %in% names(feature_metadata))) {
        # TODO:
        # Allow nonstandard evaluation of 'id' column.
        # End.
        stop("ERROR: feature_metadata does not have an 'id' column.")
    }
  
    test <- switch(class(test),
           GRanges = mcols(test) %>% data.frame,
           data.frame = test)
    
    reference <- switch(class(reference),
           GRanges = mcols(reference) %>% data.frame,
           data.frame = reference)
    
    # Parse colnames:
    bool_annotation_colname <- paste0(annotation_basename, "_bool")
    simple_annotation_colname <- paste0(annotation_basename, "_simple")
    full_annotation_colname <- paste0(annotation_basename, "_full")
    
    # Grab data
    # Replace NA with "NA" (string) so that NA values will be correctly compared
    ref_output <- reference %>% 
        dplyr::select(id, behavior_type) %>% 
        tidyr::replace_na(., replace = list(behavior_type = "NA"))

    # XXX:            
    # need a way to figure out how to handle this line if using NSE to get the idcol
    names(ref_output) <- paste0("reference_", names(ref_output))
    ref_output <- rename(ref_output, reference_id = "id")
    # End.
    
    test_output <- test %>% 
        dplyr::select(id, behavior_type) %>% 
        tidyr::replace_na(., replace = list(behavior_type = "NA"))
    names(test_output) <- paste0("test_", names(test_output))
    test_output <- rename(test_output, test_id = "id")
    
    # merge comparison data with feature information:
    behavior_comparison <- dplyr::select(feature_metadata, id)
    behavior_comparison <- behavior_comparison %>%
          dplyr::left_join(., ref_output, by = "id") %>% 
          dplyr::left_join(., test_output, by = "id")
    
    # First assign simple & full behavior descriptors: 
    behavior_comparison[[simple_annotation_colname]] <- behavior_comparison$reference_behavior_type
    
    behavior_comparison[[full_annotation_colname]] <- paste0(behavior_comparison$reference_behavior_type,
                                                             ".", behavior_comparison$test_behavior_type)
    
    # Decide which features are affected & unaffected (store in $effect_str) 
    behavior_comparison[[bool_annotation_colname]] <- (behavior_comparison$reference_behavior_type != behavior_comparison$test_behavior_type)
    
    behavior_comparison <- behavior_comparison %>% 
                           mutate(effect_str = case_when(
                               .[[bool_annotation_colname]] == FALSE ~ unaffected_prefix,
                               .[[bool_annotation_colname]] == TRUE ~ affected_prefix)) 
                               #.[[bool_annotation_colname]] == FALSE ~ "unaffected_",
                               #.[[bool_annotation_colname]] == TRUE ~ "affected_")) 
    
    # Add affected/unaffected strings to simple & full behavior descriptors
    behavior_comparison[[simple_annotation_colname]] <- paste0(behavior_comparison$effect_str,
                                                               behavior_comparison[[simple_annotation_colname]])
    behavior_comparison[[full_annotation_colname]] <- paste0(behavior_comparison$effect_str,
                                                             behavior_comparison[[full_annotation_colname]])
    
    # Drop the extraneous columns used in processing, remerge with feature_metadata & return
    comparison_data <- behavior_comparison %>% select_("id",
                                                       simple_annotation_colname, 
                                                       full_annotation_colname,
                                                       bool_annotation_colname)
    feature_metadata <- dplyr::left_join(feature_metadata, comparison_data, by = "id")
    
    return(feature_metadata)
}

annotateMAplotR <- function(results, ylim = c(-15,15), xlim = c(-2, 25), subset_features = F) {
    # Plot smoothscatter-like MA plot w/ red dots as significant.
    # Filter out features that are not from comparison (OPIONAL)
    #   Pass column name (quoted) by which to subset to subset_features (will filter out matches to Condtion1/2), 
    #   else no subsetting will occurr
    # calculate counts of each peak behavior, and output to plot using gridExtra table
    # Requires plotting with grid.arrange() if not immediately output
    
    require(magrittr)
    require(ggplot2)
    require(gridExtra)
    
    # Checks:
    # TODO:
    # Check that condition1 & condition2, behavior, behavior_type columns exist
    
    # End.
  
    results <- switch(class(results), 
                      GRanges = mcols(results) %>% data.frame,
                      data.frame = results)
    
    plotName <- unique(results$comparison)
    plotTitle <- plotName %>% gsub("vs", " vs ", .)
    
    if (subset_features == TRUE) {
            plotRes <- results %>% dplyr::filter(subset == TRUE)
        
    } else {
        plotRes <- results
    } 
    
    # Increasing & decreasing features have already been assigned based on statistical cutoff
    significantResults <- dplyr::filter(plotRes, behavior_type == "Increasing" | behavior_type == "Decreasing")
    
    MAplot <-  ggplot(plotRes, aes(log2(baseMean), log2FoldChange)) +
        stat_density_2d(aes(fill = ..density..^0.25), geom = "tile", contour = F, n = 200) +
        geom_point(data = significantResults, aes(x = log2(baseMean), y = log2FoldChange), color = "red", size = 0.1) +
        scale_fill_continuous(low = "grey94", high = "dodgerblue4") +
        ylim(ylim[1], ylim[2]) + xlim(xlim[1], xlim[2]) + theme(legend.position = "none") +
        ggtitle(plotTitle)
    
    countSummary <- dplyr::count(plotRes, behavior) %>% 
                    gridExtra::tableGrob(., rows = NULL, theme = ttheme_default(base_size = 10))
    layout <- rbind(c(1,1,1,1,2), c(1,1,1,1,2))
    AnnotatedPlot <- gridExtra::arrangeGrob(MAplot, countSummary, ncol = 5, nrow = 5, layout_matrix = layout, name = plotName)
    return(AnnotatedPlot)
}


