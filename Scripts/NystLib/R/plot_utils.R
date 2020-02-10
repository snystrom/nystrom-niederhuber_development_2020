plot_savr <- function(plot, directory, fileName, 
                      fileTypes = c(".png", ".svg", ".eps"), 
                      width = unit(10, "in"), height  = unit(7.5, "in")){
    # TODO:
    # dont plot if exists unless force = T
    # End.
    dir.create(directory)
    for (fType in fileTypes) {
        outFile = paste0(directory, fileName, fType)
        ggsave(file = outFile, plot, width = width, height = height)
    }
}



# set break midpoint to 0
# See complicated example below:
# Relevant piece of code is the definition of `myBreaks`
# from: https://stackoverflow.com/questions/31677923/set-0-point-for-pheatmap-in-r

breaks_about_zero <- function(data, palette){
        paletteLength <- length(palette)
        breaks <- c(seq(min(heatmap), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(heatmap)/paletteLength, max(heatmap), length.out=floor(paletteLength/2)))
        return(breaks)
}


cowplot_title <- function(plot, title, rel_heights = c(0.1, 1), ...){
  
  title_plot <- ggdraw() +
    draw_text(title, ...)
  
  cowplot::plot_grid(title_plot, plot, nrow = 2, rel_heights = rel_heights)
}



## ##########
## my_axis allows adding images to axis text in ggplot2
## From: https://jcarroll.com.au/2016/06/03/images-as-x-axis-labels-updated/
## ##########
# user-level interface to the element grob

#my_axis = function(img) {
#  structure(
#    list(img=img),
#    class = c("element_custom","element_blank", "element") # inheritance test workaround
#  )
#}
## returns a gTree with two children: the text label, and a rasterGrob below
#element_grob.element_custom <- function(element, x,...)  {
#  stopifnot(length(x) == length(element$img))
#  tag <- names(element$img)
#  # add vertical padding to leave space
#  g1 <- textGrob(paste0(tag, "\n\n\n\n\n"), x=x, vjust=0.6)
#  g2 <- mapply(rasterGrob, x=x, image=element$img[tag], 
#               MoreArgs=list(vjust=0.7, interpolate=FALSE,
#                             height=unit(3,"lines")),
#               SIMPLIFY=FALSE)
#  
#  gTree(children=do.call(gList, c(g2, list(g1))), cl="custom_axis")
#}
## gTrees don't know their size and ggplot would squash it, so give it room
#grobHeight.custom_axis = heightDetails.custom_axis = function(x, ...)
#  unit(6, "lines") 
### not sure if above has side effects...


## ##########
## END
## ##########
