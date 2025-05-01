global_median_norm_batch <- function(m) {
  
  sample.bias <- apply(m, 2, median, na.rm = T)
  m <- sweep(m, 2, sample.bias, "-")
  
  batch.bias <- apply(m, 1, median, na.rm = T)
  m <- sweep(m, 1, batch.bias, "-")
  
  return(m)
}



global_median_norm <- function(m) {
  
  sample.bias <- apply(m, 2, median, na.rm = T)
  m <- sweep(m, 2, sample.bias, "-")
  
  return(m)
}


plotClusters <- function(vec, nPixels, ncol, nrow, print.matrix = FALSE,
                         title = NULL, byrow = TRUE,
                         #col = colorRampPalette(brewer.pal(12, "Set3"))(n.col),
                         col = colorRampPalette(brewer.pal(9, "Set1"))(n.col), na.col = "white",
                         ...) {
  
  if (all(is.na(vec))) {
    next
  }
  
  # this is the number of colors
  n.col = length(unique(as.numeric(vec)[!is.na(as.numeric(vec))]))
  
  
  r <- raster(matrix(vec, nrow = nrow, ncol = ncol, byrow = F))
  
  if (print.matrix == TRUE) {
    print(as.matrix(r))
  }
  
  plot(as.matrix(r), xlab = '', ylab = '', xaxt = 'n', na.col = na.col,
       col = col, main = title, ...)
  #col = colorRampPalette(c("gray", "black"))(n.col), main = title, ...)
}
