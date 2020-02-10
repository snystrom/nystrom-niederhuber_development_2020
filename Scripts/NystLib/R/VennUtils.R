# Venn Diagram Utilities
# Pipeline:
# make_3_category_ov_list_from_chipPeakAnno %>% 
#   venn.diagram(., filename = NULL, <configure params>) %>% 
#   get_venn_parts() %>% 
#   plot_3color_venn()


make_3category_ov_list_from_chipPeakAnno <- function(ov_anno){
  # take chippeak anno output and convert to overlap list amenable to plotting with VennDiagram
  shared <- ov_anno$peaklist[[3]]$peakNames %>% 
    paste(., collapse = "_") %>% 
    as.character()
  list1_name <- names(ov_anno$peaklist)[[1]]
  list2_name <- names(ov_anno$peaklist)[[2]]
  
  A <- ov_anno$peaklist[[1]]$peakNames %>% 
    as.character %>% 
    c(., shared)
  B <- ov_anno$peaklist[[2]]$peakNames %>% 
    as.character %>% 
    c(., shared)
  ov_list <- list(list1_name = A, list2_name = B)
  
  return(ov_list)
}

get_venn_parts <- function(venn){
  #library(sp)
  #library(rgeos)
  # grab x- and y-values from first circle
  venn1_x <- venn[[3]][["x"]]
  venn1_y <- venn[[3]][["y"]]
  
  # grab x- and y-values from second circle
  venn2_x <- venn[[4]][["x"]]
  venn2_y <- venn[[4]][["y"]]
  
  venn1 <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(cbind(venn1_x, venn1_y))), ID = 1))) 
  venn2 <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(cbind(venn2_x, venn2_y))), ID = 2))) 
  
  intersect <- rgeos::gIntersection(venn1, venn2) 
  
  ix <- sapply(venn, function(x) grepl("text", x$name, fixed = TRUE))
  labs <- do.call(rbind.data.frame, lapply(venn[ix], `[`, c("x", "y", "label")))
  
  list("venn1" = venn1, 
       "venn2" = venn2, 
       "venn1_x" = venn1_x, 
       "venn1_y" = venn1_y, 
       "venn2_x" = venn2_x, 
       "venn2_y" = venn2_y, 
       "intersect" = intersect,
       "labs" = labs)
}

plot_3color_venn <- function(vennParts, venn_color = c(NULL, NULL, NULL), ...){
  plot(vennParts$venn1, 
       xlim = range(c(vennParts$venn1_x, vennParts$venn2_x)),
       ylim = range(c(vennParts$venn1_y, vennParts$venn2_y)), col = venn_color[1], ...)
  plot(vennParts$venn2, add = T, col = venn_color[3], ...)
  plot(vennParts$intersect, add = T, col = venn_color[2], ...)
  # todo: ifelse labs = T add labs
  #text(x = labs$x, y = labs$y, labels = labs$label)
}

# TEsting code:
#shared <- wt_ov_gof$peaklist$`OR.24APF///E93GOF.3LW`$peakNames %>% 
#  paste(., collapse = "_") %>% 
#  as.character()
#A <- wt_ov_gof$peaklist$OR.24APF$peakNames %>% as.character() %>% c(., shared)
#B <- wt_ov_gof$peaklist$E93GOF.3LW$peakNames %>% as.character() %>% c(., shared)
#
#ov_list <- list("A" = A, "B" = B)
#venn <- venn.diagram(ov_list, filename = NULL, 
#                           category = c("Endogenous\nE93 24APF ChIP", "Exogenous\nE93 3LW ChIP"),
#                           fontfamily = "sans",
#                           cat.pos = c(200,180), 
#                           cat.dist = c(.08, 0.02),
#                           cat.fontfamily = "sans",
#                           cat.face = "bold", 
#                           ext.length = 0.9, 
#                           ext.line.lwd = 2) 

#vennParts <- get_venn_parts(venn)
#plot_3color_venn(vennParts, c("red", "green", "blue"))
#
#grid.newpage()
#grid.draw(venn)