subset_features_by_condition <- function(results, features){
  # output Granges object with metadata of results by "id" (or other identifier column once implemented).
  ######
  # TODO:
  # Requires id column in features (update to use NSE)
  condition1 <- unique(results$condition1) 
  condition2 <- unique(results$condition2) 

  matchCols<- c(condition1, condition2) %>% 
              paste(., collapse="|")

  featureIDs <- mcols(features) %>%
             data.frame(.) %>% 
             dplyr::filter(grepl(matchCols, id)) %>% 
             .$id

  subset_features <- features[mcols(features)$id %in% featureIDs]

  mcols(subset_features) <- mcols(subset_features) %>%
                          data.frame(.) %>%
                          dplyr::left_join(., results, by = "id")
  return(subset_features)
}

get_heatmap_from_seqplots <- function(signal_data, feature_id_colname){
  # return melted dataframe containing positional signal for each feature
  melted_list <- lapply(signal_data$data, function(feature){
    data <- lapply(feature, function(sample){
      
      featureIDs <- mcols(sample$anno)[[feature_id_colname]]
      
      colnames(sample$heatmap) <- sample$all_ind
      rownames(sample$heatmap) <- featureIDs
      countData_melt <- reshape2::melt(sample$heatmap, value.name = "signal", varnames = c("feature", "position")) 
      countData_melt$sample.desc <- sample$desc 
      
      return(countData_melt)
    })
    return(dplyr::bind_rows(data))
  })
  
  melted_data <- dplyr::bind_rows(melted_list)
  return(melted_data)
}

get_avgSignal_from_seqplots <- function(signal_data){
  # Return average-signal style data from plotSeqArray object
  melted_list <- lapply(signal_data$data, function(feature){
    data <- lapply(feature, function(sample){
      countData <- data.frame(sample$mean, sample$all_ind, sample$desc)
      names(countData) <- c("mean", "position", "sample.desc")
      return(countData)
    })
    return(dplyr::bind_rows(data))
  })
  
  melted_data <- dplyr::bind_rows(melted_list)
  return(melted_data)
}

get_featureSignalSummary_from_seqplots <- function(signal_data, feature_id_colname){
  # Extract signal information across features for each subset of data
  # Return as melted dataframe
  melted_list <- lapply(signal_data$data, function(feature){
    data <- lapply(feature, function(sample){
      signal_sum <- rowSums(sample$heatmap)
      signal_mean <- rowMeans(sample$heatmap)
      sample.desc <- rep(sample$desc, length(signal_sum))
      
      # Using actual feature ID from $anno object. 
      # My understanding is that order is always preserved in this case. 
      # If region is out of bounds after extending by xmin & max, then signal = NA
      # (ie if signal can't be calculated, the region is not dropped from output, but given NA value)
      # so, this should be a robust solution.
      feature <- mcols(sample$anno)[[feature_id_colname]]
      
      countData <- data.frame(signal_sum, signal_mean, feature, sample.desc)
      return(countData)
    })
    return(dplyr::bind_rows(data))
  })
  
 melted_data <- dplyr::bind_rows(melted_list)
 return(melted_data) 
}

get_enrichedHeatmap_from_seqplots <- function(signal_data, feature_id_colname, xmin, xmax, binsize){
  # Extract signal information across features for each subset of data
  # Return as melted dataframe
  
  if (length(signal_data$data) > 1) {
    stop("ERROR: input ranges must not be split ranges (not GRangesList) for enrichedHeatmap output. Instead, split features after computation.")
  }
  
  data <- lapply(signal_data$data$feature_1, function(sample){
    ranges <- sample$anno 
    matrix <- sample$heatmap 
    
    mcols(ranges)$matrix_signal.sum <- rowSums(matrix) 
    mcols(ranges)$matrix_signal.mean <- rowMeans(matrix) 
   
    
    # TODO:
    # convert type to normalized matrix
    
    normMatrix <- matrix
    attributes(normMatrix) <- NULL
    attr(normMatrix, "upstream_index") <- 1:(xmin/binsize)
    attr(normMatrix, "target_index") <- integer(0)
    attr(normMatrix, "downstream_index") <- 1:((xmax/binsize) + 1) + (xmin/binsize)
    attr(normMatrix, "extend") <- c(xmin, xmax)
    attr(normMatrix, "target_name") <- "center"
    # signal_name = track
    attr(normMatrix, "signal_name") <- sample$desc
    dim(normMatrix) <- dim(matrix)
    class(normMatrix) <- c("normalizedMatrix", "matrix")

    # End.
    
     
    # My understanding is that order is always preserved in this case. 
    # If region is out of bounds after extending by xmin & max, then signal = NA
    # (ie if signal can't be calculated, the region is not dropped from output, but given NA value)
    # so, by sorting the ranges and applying that sorted vector to the matrix, it can be sorted/split in various ways, etc.
    sample_rangedMatrix <- list("ranges" = data.frame(ranges), "matrix" = normMatrix) 
    
    return(sample_rangedMatrix)
    #return(dplyr::bind_rows(data))
  })
  
 #melted_data <- dplyr::bind_rows(melted_list)
 #return(melted_data) 
 return(data)
}

get_signal_data <- function(sampleData, tracks, features, feature_id_colname = "id", output_type = "avgSignal",
                                     refgenome = 'dm3', xmin = 1000, xmax = 1000, binsize = 10, type = "mf"){
  # This is a wrapper around the seqplot::getPlotSetArray() function to return meaningful subsets of the data
  # which are annotated with any metadata contained in sampleData (a metafile).
  # output_type can be one of: "avgSignal", "featureSummary", "Heatmap" (case insensitive).
  # Where:
  # output_type = "avgSignal": 
    # return melted dataframe of sum(signal) and mean(signal) across features for each dataset.
  # output_type = "featureSummary":
    # This is not for average Signal plots, but instead for comparing total signal within a featureset across many datasets.
  # output_type = "Heatmap":
    # This returns position-specific signal information for each feature.
  
  # each feature gets unique ID so signal in features can be compared across datasets. 
  # where feature_id_colname is the column name of the GRanges metadata that stores the feature id (must be unique)
  
  # where tracks is the vectorized column of sampleData (vector of bigwigs, typically)
  # so usage should be: sampleData = metafile, tracks = metafile$bigwigs
  
  # TODO:
  #  Check: tracks is vector contained within sampleData
  #  FUNCTION: parse refgenome & load correct library
  # End.
  
  
  # Checks:
  if (class(sampleData) != "data.frame"){
    stopString <- paste0("ERROR: sampleData must be a data.frame! sampleData class is currently: ", class(sampleData), ".")
    stop(stopString)
  }
  
  if (class(unlist(features)) != "GRanges"){
    stop("ERROR: features is not a GRanges object or GRangesList!")
  }
  
  if (!(feature_id_colname  %in% names(mcols(unlist(features))))) {
    stopString <- paste0("ERROR: feature_id_colname (", feature_id_colname, ") does not exist in features metadata.")
    stop(stopString)
  }
 
  feature_id_check <- mcols(unlist(features))[[feature_id_colname]]  
  if (length(feature_id_check) != length(unique(feature_id_check))){
    stop("ERROR: featureIDs (in: ", feature_id_colname, " of `features`) are non-unique.")
  }
  
  if (!(tolower(output_type) %in% c("avgsignal", "featuresummary", "heatmap", "enrichedheatmap"))){
    stop("ERROR: unknown output_type: ", output_type)
  }
  
  if (!(tolower(type) %in% c("mf", "af", "pf", "ef"))){
    stop(paste0("ERROR: unknown type: ", type, " ! Must be one of: mf, af, pf, ef."))
  }
  
  # XXX:
  # TODO:
  # Support all organism genomes
  library(BSgenome.Dmelanogaster.UCSC.dm3)
  
  signal_data <- seqplots::getPlotSetArray(tracks = tracks, features = features, refgenome = refgenome,
                  xmin = xmin, xmax = xmax, bin = binsize, type = type)
  detach("package:BSgenome.Dmelanogaster.UCSC.dm3", unload = TRUE)
  # End.
   
  parsed_signal_data <- switch(tolower(output_type),
                       avgsignal = get_avgSignal_from_seqplots(signal_data),
                       featuresummary = get_featureSignalSummary_from_seqplots(signal_data, feature_id_colname),
                       heatmap = get_heatmap_from_seqplots(signal_data, feature_id_colname),
                       enrichedheatmap = get_enrichedHeatmap_from_seqplots(signal_data, feature_id_colname, xmin, xmax, binsize),
                       stop("ERROR: unknown output_type: ", output_type, " escapes checks, but has no corresponding function.")
                       )

 # Now add sample metadata to each signal entry by merging on signal-file (track) used to generate counts. 
 # Requires first removing path and filetype information of signal file in `tracks` variable 
 # and merging this into sampleData stored in new variable $track.
 sampleData$track <- tools::file_path_sans_ext(basename(tracks))

 # sample.desc contains <track>\n@<feature_description>, 
 # so these need to be split out into two columns (track & desc, respectively) before merging Sample metadata by $track
 
 # TODO:
 # Append feature-level metadata to this as well? 
 # End.
 
 melt_final_data <- function(parsed_signal_data, sampleData){
  finalData <- parsed_signal_data %>% 
               mutate(track = stringr::str_split_fixed(sample.desc, '\n@', 2)[,1]) %>% 
               mutate(desc = stringr::str_split_fixed(sample.desc, '\n@', 2)[,2]) %>% 
               dplyr::left_join(., sampleData, by = "track")
  return(finalData)
 }
 
 rename_tracks_to_sample <- function(parsed_signal_data, sampleData){
   # enrichedHeatmap output is in form: data$trackName$list(ranges, matrix)
   # change trackName to sampleName
   if (! identical(names(parsed_signal_data), sampleData$track)){
     stop("ERROR: heatmap sample names do not match track names")
   }
   
   names(parsed_signal_data) <- sampleData$Sample
   return(parsed_signal_data)
 }
 
 finalData <- switch(class(parsed_signal_data),
                     data.frame = melt_final_data(parsed_signal_data, sampleData),
                     list = rename_tracks_to_sample(parsed_signal_data, sampleData)) 
 
 
 return(finalData)
}

toScatterPlot <- function(featureSummary, key = "Genotype", value = "signal_sum", keep = "desc"){
  # featureSummary is the output of `get_signal_data(output_type = "featureSummary")`
  # will use tidyr to pivot table for pairwise comparisons between conditions (where key = condition)
 
  # keep = column names to include in summary table (can be NA defaults to desc)
  # TODO:
  # add checks to make sure this doesn't generate bogus tables? ie confounding
  # variables that are diferent between conditions will result in unusable melts
  # because 'NA' will be introduced into each cast column because features will
  # separate by confoudning varible.
  #
  # In short, the only data that should be in the 'keep' variable are columns that describe FEATURES, not conditions (ie higher_condition1/2)
  if (!is.na(keep[1])){
    subset <- featureSummary %>% dplyr::select_("feature", key, value, .dots = keep)
  } else {
    subset <- featureSummary %>% dplyr::select_("feature", key, value)
  }
    
  scatterData <- tidyr::spread_(subset, key = key, value = value)
  
  return(scatterData)
}

split_sample.desc <- function(avgSignal, n = 2, on = "///", col_names = paste0("desc", seq_len(n))){
  # if want avgSignal plots by multiple categories,
  # ex behavior_bound:
  # then make GRList in the form `by = behavior///bound` (GRanges %>% split(., mcols(.)$by))
  # run get_signal_data
  # then run split_sample.desc on the resulting data.frame (will return new data frame)
  # 
  # Input: avgSignal = data.frame ouput from get_signal_data()
  # n = number of categories stored in sample.desc to split
  # on = delimiter of categories
  # col_names = names of each column to be created (default `desc<n_index>`)
  descCols <- avgSignal$sample.desc %>% 
    stringr::str_split_fixed("\n@", n = 2) %>% 
    .[,2] %>% 
    stringr::str_split_fixed(on, n = n) %>% 
    data.frame
    
  names(descCols) <- col_names
  avgSignal <- cbind(avgSignal, descCols)
  return(avgSignal)
}
