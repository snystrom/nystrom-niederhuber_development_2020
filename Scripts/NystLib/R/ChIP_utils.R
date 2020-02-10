annotate_bound_peaks <- function(results, resPeaks, chipPeaks, colName, idcol = "id"){
  # Input: results (dataframe), resPeaks (GRanges), chipPeaks (GRanges), colname (string)
  # REQUIRES: results and resPeaks to have an 'id' column that allows joining peaks.
  # Typical usage: resPeaks will be unionPeaklist. 
  # Function will subset the unionpeaks to those with matching ids in results dataframe.
  # OUTPUT: new column in results (named `colName`) that reports 'bound' vs. 'unbound' based on simple GRanges overlap
  
  peaks <- resPeaks[mcols(resPeaks)[[idcol]] %in% results[[idcol]]]
  ov <- peaks[peaks %over% chipPeaks]
  
  results[[colName]] <- "unbound"
  
  results[results[[idcol]] %in% mcols(ov)[[idcol]],][[colName]] <- "bound"
  return(results)
}

annotate_peak_features <- function(peaks, min = 0, max = 0, 
                       txDb = TxDb.Dmelanogaster.UCSC.dm3.ensGene::TxDb.Dmelanogaster.UCSC.dm3.ensGene,
                       orgDb = org.Dm.eg.db::org.Dm.eg.db, annoDb = "org.Dm.eg.db"){
  # input: peaks (GRanges), which will be annotated to nearest gene, feature, etc. 
  # returns all chipSeeker annotations along with another column: `annotation_simple` which strips gene-level information
  # from "Intron/Exon (geneid ...)" in `annotation` column for easy grouping & plotting
  
  # output: dataframe with all metadata
  # usage: meta(faire_data)  <- annotate_peak_features(metaRanges(faire_data))
  
  GRcolnames <- c("seqnames", "start", "end", "width", "strand")
  anno_colnames <- c("annotation", "annotation_simple", "geneChr", "geneStart", "geneEnd", "geneLength", "geneStrand", 
                     "geneId", "transcriptId", "distanceToTSS", "ENTREZID", "SYMBOL", "GENENAME") 
  
  annoPeaks <- ChIPseeker::annotatePeak(peaks, tssRegion = c(min, max), TxDb = txDb, annoDb = annoDb) %>% 
    ChIPseeker::as.GRanges(.)
 
  # drop GRanges position data 
  annotation_data <- annoPeaks %>% 
    data.frame %>%
    dplyr::mutate(annotation_simple = gsub(" \\(.*\\)", "", annotation)) %>% 
    #dplyr::select_(one_of(paste0("-", GRcolnames)))
    .[-c(1:5)] 
   
  return(annotation_data)
}
