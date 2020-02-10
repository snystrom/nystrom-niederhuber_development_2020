# Pipeline Functions
setup_packages <- function(n_threads = 4){
  library(magrittr)
  library(dplyr)
  library(modelr)
  library(purrr)
  library(ggplot2)
  library(BiocParallel)
  register(MulticoreParam(n_threads))
  library(GenomicRanges)
  
  nyLibDir <- "Scripts/Nystlib/R/"
  nyLibsrc <- paste0(nyLibDir, list.files(nyLibDir))
  sapply(nyLibsrc, source, .GlobalEnv)
}

prepare_sample_sheet <- function(path = "SampleSheet.tsv"){

  # Prepare sample sheet
  sampleSheet <- read.delim(path, stringsAsFactors = F)
  sampleSheet$Sample <- paste(sampleSheet$Genotype, sampleSheet$Time, sampleSheet$Tissue, sampleSheet$Rep, sep = ".")
  sampleSheet$Peaks <- paste(sampleSheet$Genotype, sampleSheet$Time, sampleSheet$Tissue, sep = ".")
  #-----
  # I'm hard-coding the timecourse order to avoid any weird behavior later.
  # Could also do levels = as.vector(sampleSheet$Time) if samples were in correct order in sampleSheet
  sampleSheet$Time <- factor(sampleSheet$Time, levels = c("3LW", "24APF", "44APF"))
  sampleSheet$Genotype <- factor(sampleSheet$Genotype, levels = c("OR", "E93mut", "E93GOF")) 
  sampleSheet <- dplyr::arrange(sampleSheet, Genotype, Time, Rep)
  
  #sampleSheet$peaks <- paste0() 
   
  #sampleSheet$Sample <- factor(sampleSheet$Sample, levels = sampleSheet$Sample)
  sampleSheet$grp <- paste0(sampleSheet$Genotype, ".", sampleSheet$Time)
  sampleSheet$grp <- factor(sampleSheet$grp, levels = unique(sampleSheet$grp))
  
  return(sampleSheet)
}

prepare_peak_table <- function(path = "PeakTable.tsv", sampleSheet, grp_level = "OR.3LW"){
  peakTable <- read.delim(path, stringsAsFactors = F)
  peakTable$Sample <- paste(peakTable$Genotype, peakTable$Time, peakTable$Tissue, sep = ".")
  peakTable$grp <- factor(paste(peakTable$Genotype, peakTable$Time, sep = "."))
  
  peakTable$Time <- factor(peakTable$Time, levels = levels(sampleSheet$Time))
  peakTable$grp <- relevel(peakTable$grp, grp_level)
  return(peakTable)
}


###
# Computing differences in Z-score between faire-datasets
###


compute_deltaZ <- function(peaks, sampleInfo, condition1, condition2, 
                           by = "grp", 
                           signal_file_colname = "BigwigZnorm", 
                           summary_type = "mean"){
  # computes change in z score (condtion1 - condition2)
  # returns `peaks` (GRanges) with deltaZ column appended to metadata in the form: dZ_<condition1>_<condition2>
  #
  # condition1/2 = character vector of sample id
  # by = colname that contains condition-level information 
  # (ie if c1/2 are in a column called "grp", by = "grp")
   
  sampleList <- sampleInfo %>% 
    split(., .[[by]])
  
  if (!(condition1 %in% names(sampleList))) {
    stop(paste0("Error: condition1 (", condition1,") is not contained in the column: ", by))
  }
  
  if (!(condition2 %in% names(sampleList))) {
    stop(paste0("Error: condition2 (", condition2,") is not contained in the column: ", by)) 
  } 
    
  #c1_bw_path <- sampleList[[condition1]] %>% .[[signal_file_colname]]
    
  #c1_bw <- rtracklayer::BigWigFile(c1_bw_path)
  
  c1_bw <- sampleList[[condition1]][[signal_file_colname]] %>% 
   rtracklayer::BigWigFile(.)
  
  c2_bw <- sampleList[[condition2]][[signal_file_colname]] %>% 
   rtracklayer::BigWigFile(.)
  
  c1_signal <- summary(c1_bw, peaks, type = summary_type) %>% 
    unlist
  c2_signal <- summary(c2_bw, peaks, type = summary_type) %>% 
    unlist
  
  dZ_col.name <- paste0("dZ_", condition1, "_", condition2)
  
  mcols(peaks)[[dZ_col.name]] <- mcols(c1_signal)$score - mcols(c2_signal)$score
  return(peaks)
}

#out <- compute_deltaZ(test_peaks, peakTable, "E93.GOF", "OR.3LW")

pairwise_deltaZ <- function(comparisons, sampleInfo, peaks, signal_file_colname = "BigwigZnorm", summary_type = "mean"){
  # returns input peaks with delta_Z appended to metadata for each pairwise comparison
  # comparisons = output from pairwise_comparisons() function
  # sampleInfo = sampleSheet/metafile containing sample-specific information
  # signal_file_colname = column in sampleInfo that contains a path to the Z-normalized bigwig (MUST BE BIGWIG) 
  # summary_type = any of rtracklayer::BigWigFile summary() types. Untested for values other than "mean".
  
  # TODO: could expand this fuction for uses other than dZ by adding "prefix" variable instead of "dZ_" to each colname
  # nothing about it actually prevents use for signal comparison elsewhere...
  
  # WARNING: 
  # THERE IS NO TYPE CHECKING YET
  
  # BETTER WAY TO DO THIS:
  # Compute z-scores for each dataset once into a data frame, then do the pairwise subtractions and add those columns to metadata
  # would be significantly faster. current version does too much I/O.
  
  dZ <- lapply(comparisons, function(comparison){
    # returns list of colname & deltaZ
    # where colname is the deltaZ comparison column name to be added to the peak metadata
    # deltaZ is in the same order as the peaks in peaklist
    filter_str_c1 <- paste0(comparison$contrast[1], " == \"", comparison$condition1, "\"")
    filter_str_c2 <- paste0(comparison$contrast[1], " == \"", comparison$condition2, "\"")
    
    # make bigwig file connection, then summarize Zscores within peaks (default mean-Z) 
    c1_bw <- sampleInfo %>% 
      dplyr::filter_(filter_str_c1) %>% 
      .[[signal_file_colname]] %>% 
      rtracklayer::BigWigFile(.)
    c1_score <- summary(c1_bw, peaks, type = summary_type) %>% 
      unlist %>% 
      mcols(.) %>% 
      .$score
    
    c2_bw <- sampleInfo %>% 
      dplyr::filter_(filter_str_c2) %>% 
      .[[signal_file_colname]] %>% 
      rtracklayer::BigWigFile(.)
    c2_score <- summary(c2_bw, peaks, type = summary_type) %>% 
      unlist %>% 
      mcols(.) %>% 
      .$score
    
    
    delta_z <- c1_score - c2_score
    
    # create new column to hold deltaZ values:
    # parse from contrast name:
    # dZ_<condition1>_<condition2>
    # where c1 - c2 = dZ
    dZ_col.name <- paste0("dZ_", comparison$condition1, "_", comparison$condition2)
    dZ_df <- data.frame(delta_z) 
    names(dZ_df) <- dZ_col.name
    return(dZ_df)
  })
  
  dZ_cols <- do.call("cbind", dZ)
  # append deltaZ columns to metadata
  mcols(peaks) <- c(mcols(peaks), dZ_cols)
  return(peaks)
}


overlap_pairwise_deltaZ <- function(overlappingPeaks, comparisons, sampleInfo, 
                                    signal_file_colname = "BigwigZnorm",
                                    summary_type = "mean"){
  
  # function takes pairwise deltaZ of chippeakanno overlap object
  # input is overlappingPeaks object from ChIPpeakAnno::findOverlapsOfPeaks()
  # comparisons = output of pairwise_comparisons()
  # sampleInfo = metafile describing datasets (should be same file used to generate comparisons)
  # signal_file_colname = column containing path to bigwig file to use
  # summary_type = type of bigwig summary to take within each range (default = mean, untested with other values)
  
  if (class(overlappingPeaks) != "overlappingPeaks"){ stop("Error: overlappingPeaks must be of class 'overlappingPeaks'")} 
    
  allPeak_overlaps_dz <- lapply(names(overlappingPeaks$peaklist), function(name){
    peaks <- overlappingPeaks$peaklist[[name]]
    #peaks <- compute_deltaZ(peaks, peakTable, "E93GOF.3LW", "OR.3LW", by = "grp")
    #peaks <- compute_deltaZ(peaks, peakTable, "OR.24APF", "OR.3LW", by = "grp")
    peaks <- pairwise_deltaZ(comparisons, 
                             sampleInfo = sampleInfo, 
                             peaks = peaks, 
                             summary_type = summary_type, 
                             signal_file_colname = signal_file_colname)
    mcols(peaks)$peak_type <- name
    return(data.frame(peaks))
  }) %>% 
    do.call("rbind", .)
  return(allPeak_overlaps_dz)
}

drop_na <- function(df, colName = "peak_binding_description", naVal = "na"){
  # for dropping the peaks shared between wt and gof to avoid doubleCounting
  filter_str = paste0(colName, " != \"", naVal, "\"")
  df %<>% 
    dplyr::filter_(filter_str)
  return(df)
}
