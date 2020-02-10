# GOutils

getGenesByGO <- function(go, term, by = "geneID"){
  # go = enrichGO object
  # term = go term description
  # by = output
  genes <- go@result %>% 
    data.frame %>%
    .[.$Description == term,] %>% 
    .[[by]] %>% 
    stringr::str_split(., "/") %>% 
    unlist
  return(genes)
}

goBarPlot <- function(go_out, n = 30, name = NULL){
  # barplot of go results from clusterProfiler enrichGO output
  go_out %>% 
    head(n) %>% 
    ggplot(aes(reorder(Description, -log10(p.adjust)), -log10(p.adjust))) +
      geom_col() + 
      coord_flip() +
      xlab(NULL) +
      ggtitle(name) +
      theme_cowplot()
}
