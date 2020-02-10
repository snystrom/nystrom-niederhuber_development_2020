source("R_Pipeline/pipeline_methods.R")

main <- function(sample_sheet_path = "SampleSheet.tsv", peak_table_path = "PeakTable.tsv", 
                 n_threads = 4, faire_qValueCutoff = 40){
  setup_packages(n_threads = n_threads) 
  sampleSheet <- prepare_sample_sheet(path = sample_sheet_path) %>% 
    dplyr::mutate(Assay = "FAIRE")
  peakTable <- prepare_peak_table(path = peak_table_path, sampleSheet) %>% 
    dplyr::mutate(Assay = "FAIRE")
 
  if (length(peakTable$Sample) != length(unique(peakTable$Sample))) {
     stop("peakTable SampleID's are non-unique!")
  }
  # TODO: 
  # add path environment variable for writing out these data, append & save filepath, etc.
  # for now will just writeout to R_Pipeline directory 
  # (check dir exist, make if not, etc. also needs add w/ this feature)
  # End.
  save(sampleSheet, peakTable, file = "R_Pipeline/output_data/metadata.rda")
  
  # import all peaks for each replicate and pooled peaks 
  allPeaks_by_rep <- getPeakData(sampleSheet, by = "Sample", narrowPeak_colname = "Peakfile") %>% 
    GRanges()
  
  allPeaks_pooled <- getPeakData(peakTable, by = "Sample", narrowPeak_colname = "Peakfile") %>% 
    GRanges()
   
  # write out peak lists: 
  # TODO: better file output parsing
  writeGRanges(allPeaks_by_rep, "R_Pipeline/output_data/allPeaks_by_rep.GRanges", bed = F)
  writeGRanges(allPeaks_pooled, "R_Pipeline/output_data/allPeaks_pooled.GRanges", bed = F)
  save(allPeaks_by_rep, file = "R_Pipeline/output_data/allPeaks_by_rep.rda")
  save(allPeaks_pooled, file = "R_Pipeline/output_data/allPeaks_pooled.rda")
  
  ###
  # Filter peaks by qValue cutoff ----------
  ###
  allPeaks_pooled_grpSplit <- allPeaks_pooled %>% split(., mcols(.)$grp)
  allPeaks_pooled_qCutoff <- lapply(allPeaks_pooled_grpSplit, function(peaks){
     peaks.df <- data.frame(peaks) 
     peaks_qCutoff <- peaks.df %>% 
       dplyr::filter(qValue >= faire_qValueCutoff) %>% 
       GRanges()
     
     return(peaks_qCutoff)
  })
  
  nPeak_cutoff <- sapply(allPeaks_pooled_qCutoff, length) %>% min
  print(paste0("Using top ", nPeak_cutoff, " Peaks"))
  
  allFairePeaks <- lapply(allPeaks_pooled_qCutoff, function(peaks){
    peaks_cutoff <- data.frame(peaks) %>% 
      dplyr::arrange(desc(qValue)) %>% 
      .[1:nPeak_cutoff,] %>% 
      GRanges()
  }) %>% 
    GRangesList
  
  ###
  # Make union peak list ---------- 
  ###
 
  unionFairePeaks <- allFairePeaks %>% 
    unlist %>% 
    reduce 
  
  save(allFairePeaks, unionFairePeaks, file = "R_Pipeline/output_data/final_faire_peaks.rda")
  
  ### 
  # ChIP ----------
  ###
  
  chip_sampleSheet <- prepare_sample_sheet("ChIP_SampleSheet.tsv") %>% 
    dplyr::mutate(Assay = "ChIP")
  chip_peakTable <- prepare_peak_table("ChIP_PeakTable.tsv", chip_sampleSheet, grp_level = "OR.24APF") %>% 
    dplyr::mutate(Assay = "ChIP")
  
  chip_allPeaks_by_rep <- getPeakData(chip_sampleSheet, by = "Sample", narrowPeak_colname = "Peakfile") %>% 
    GRanges()
  
  chip_allPeaks_pooled <- getPeakData(chip_peakTable, by = "Sample", narrowPeak_colname = "Peakfile") %>% 
    GRanges()
  
  save(chip_sampleSheet, chip_peakTable, file = "R_Pipeline/output_data/chip_metadata.rda")
  save(chip_allPeaks_by_rep, chip_allPeaks_pooled, file = "R_Pipeline/output_data/chip_allPeaks.rda")
  
  chip_rep_ol <- peak_overlaps_by_condition(chip_allPeaks_by_rep %>% split(., mcols(.)$grp)) 
  
  # output only merged peaks 
  allChipPeaks_repMerge <- lapply(chip_rep_ol, function(x) {x$mergedPeaks}) 
  
  # return peaks from pooled dataset that overlap RepMerged peaks 
  chip_allPeaks_pooled_grpSplit <- chip_allPeaks_pooled %>%
    split(., mcols(.)$grp)
  
  allChipPeaks_repMergeOvPool <- lapply(names(allChipPeaks_repMerge), function(name){
    allPoolPeaks <- chip_allPeaks_pooled_grpSplit[[name]]
    repMergePeaks <- allChipPeaks_repMerge[[name]]
    
    pool_ov_peaks <- subsetByOverlaps(allPoolPeaks, repMergePeaks)
    
    return(pool_ov_peaks)
  })
  names(allChipPeaks_repMergeOvPool) <- names(allChipPeaks_repMerge)
  
  # output pooled peak calls
  allChipPeaksList <- chip_allPeaks_pooled_grpSplit
  save(allChipPeaksList, allChipPeaks_repMergeOvPool, chip_allPeaks_pooled_grpSplit, file = "R_Pipeline/output_data/final_chip_peaks.rda")
}

# FAIRE
main(faire_qValueCutoff = 20)

