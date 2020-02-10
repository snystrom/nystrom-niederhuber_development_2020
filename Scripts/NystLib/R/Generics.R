# Accessor & Setter Functions ----------  
sampleInfo <- function(x) UseMethod("sampleInfo")

peakInfo <- function(x) UseMethod("peakInfo")

'peakInfo<-' <- function(x, value) UseMethod('peakInfo<-', x)

meta <- function(x) UseMethod("meta") 

'meta<-' <- function(x, value) UseMethod('meta<-', x) 

metaRanges <- function(x) UseMethod("metaRanges") 

#'metaRanges<-' <- function(x) UseMethod('metaRanges<-', x) 

sampleRanges <- function(x) UseMethod("sampleRanges") 

dds <- function(x) UseMethod("dds") 

res <- function(x) UseMethod("res") 

'res<-' <- function(x, value) UseMethod('res<-', x) 

setInfo <- function(x) UseMethod("setInfo")

misc <- function(x) UseMethod("misc")

#'misc<-' <- function(x, value) UseMethod('misc<-', x)

summits <- function(x) UseMethod("summits")

# Functions ----------  

annotate_results_subsets <- function(peakExperimentSet) UseMethod("annotate_results_subsets")

get_sample_summits_of_affected_peaks <- function(peakExperimentSet, ...) UseMethod("get_sample_summits_of_affected_peaks")

computeContrastSummits <- function(peakExperimentSet, ...) UseMethod("computeContrastSummits")

computeSummits <- function(peaks, ...) UseMethod("computeSummits")

generate_DESeq2_contrast_strings <- function(peakExperimentSet, ...) {
  # TODO:
  # add method for dataframe & RNAExperimentSet
  UseMethod("generate_DESeq2_contrast_strings")
}

generate_DESeq2_contrasts <- function(peakExperimentSet, ...) {
  # TODO:
  # add method for dataframe & RNAExperimentSet
  UseMethod("generate_DESeq2_contrasts")
}

subset_by_condition <- function(results, ...){
  UseMethod("subset_by_condition")
}

topQPeakNames <- function(peaks, ...){
  UseMethod("topQPeakNames")
}

readNarrowPeak <- function(path, ...){
  UseMethod("readNarrowPeak")
}


