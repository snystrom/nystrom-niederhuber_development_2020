peakExperimentSet <- function(DESeqDataSet, allResults = list(), sampleFeatures, metaFeatures, peakInfo = data.frame){
  # misc slot = list can store whatever. Is never checked. User can store any type of data here to
      # to keep it flat.
  # allResults = list of contrasted GRanges results objects
  if (!class(DESeqDataSet) == "DESeqDataSet") {
    if (!is.null(class(DESeqDataSet))){
      stop(paste0("ERROR: DESeqDataSet is not of class = DESeqDataSet. Currently: ", class(DESeqDataSet),
                  ". If you want this slot to be empty, use DESeqDataSet = NULL"))
    }
  } 
  
  if (!class(sampleFeatures) == "GRangesList") {
    stop(paste0("ERROR: sampleFeatures is not of class = GRangesList Currently: ", class(sampleFeatures)))
  } 
  
  if (!class(metaFeatures) == "GRanges") {
    stop(paste0("ERROR: metaFeatures is not of class = GRanges Currently: ", class(metaFeatures)))
  } 
  
  if (!class(allResults) == "list") {
    stop(paste0("ERROR: allResults is not of class = list Currently: ", class(allResults)))
  } 
  
  if (!class(peakInfo) == "data.frame") {
    stop(paste0("ERROR: peakInfo is not of class = list Currently: ", class(allResults)))
  } 


  setInfo = list(contrast_strings = list())
  structure(list(DESeqDataSet = DESeqDataSet,
                 allResults = allResults,
                 sampleFeatures = sampleFeatures,
                 metaFeatures = metaFeatures,
                 setInfo = setInfo,
                 peakInfo = peakInfo,
                 misc = list()), 
            class = "peakExperimentSet")
}

# TODO:
# is.peakData method, etc.
#
