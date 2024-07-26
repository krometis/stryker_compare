DOTE_theme = theme(panel.background = element_rect(fill = "white", color = 'black', linewidth = 2),
                   panel.grid.major = element_line(color = "gray", linetype = "dashed"),
                   aspect.ratio = 1)

DOTE_hist <- function(data, xlab=NULL, ylab=NULL, title=NULL, breaks=NULL, fill = 'black', alpha = 1){
  p = ggplot() +
    geom_histogram(aes(x=data, y=after_stat(density)), fill = fill,
                   color = 'black', breaks = breaks, alpha = alpha) +
    xlab(xlab) + ylab(ylab) + ggtitle(paste(title))
  
  return(p)
}
