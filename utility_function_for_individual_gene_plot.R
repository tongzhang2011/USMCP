
# create a meta data: location of each sample

samples <- 1:493
nrow <- 17
ncol <- 29



# matrix(1:72, nrow = nrow, ncol = ncol, byrow = TRUE)

spatial_protein_plot_log2_scale <- function(my_gene) {
  
  
  dt <- d1[my_gene,]
  
  
  # dt <- 2^dt
  # max <- max(dt, na.rm = T)
  # dt <- 100* dt / max
  
  exp_data <- as.numeric(dt)
  
  exp_data <- exp_data - median(exp_data, na.rm = T)
  
  # put the data into a matrix
  pm <- matrix(exp_data, nrow = nrow, ncol = ncol)
  
  
  data <- pm
  
  
  # Generate heatmap using pheatmap with handling of NA values
  p <-   pheatmap(data, scale = "none", 
                  # breaks = breaksList,
                  # color = my_palette_with_na, 
                  cluster_rows = F,
                  cluster_cols = F,
                  #display_numbers = mat2,
                  legend = T,
                  main = gsub("_MOUSE", "", my_gene))
  return(p)
  
}

spatial_protein_plot_log2_scale_silent <- function(my_gene) {
  
  
  dt <- d1[my_gene,]
  
  
  # dt <- 2^dt
  # max <- max(dt, na.rm = T)
  # dt <- 100* dt / max
  
  exp_data <- as.numeric(dt)
  
  exp_data <- exp_data - median(exp_data, na.rm = T)
  
  # put the data into a matrix
  pm <- matrix(exp_data, nrow = nrow, ncol = ncol)
  
  
  data <- pm
  
  
  # Generate heatmap using pheatmap with handling of NA values
  p <-   pheatmap(data, scale = "none", 
                  # breaks = breaksList,
                  # color = my_palette_with_na, 
                  cluster_rows = F,
                  cluster_cols = F,
                  #display_numbers = mat2,
                  legend = T,
                  main = gsub("_MOUSE", "", my_gene),
                  silent = TRUE)
  return(p)
  
}




################ for grid. arrange ###########

library(pheatmap)
library(gridExtra)
library(grid)

spatial_protein_plot_log2_scale_grid <- function(my_gene) {
  
  dt <- d1[my_gene,]
  
  
  # dt <- 2^dt
  # max <- max(dt, na.rm = T)
  # dt <- 100* dt / max
  
  exp_data <- as.numeric(dt)
  
  exp_data <- exp_data - median(exp_data, na.rm = T)
  
  # put the data into a matrix
  pm <- matrix(exp_data, nrow = nrow, ncol = ncol)
  
  
  data <- pm
  
  # Generate heatmap using pheatmap and store the output in a variable
  p <- pheatmap(data, scale = "none", 
                cluster_rows = FALSE, 
                cluster_cols = FALSE,
                legend = TRUE,
                main = gsub("_MOUSE", "", my_gene), 
                silent = TRUE)  # 'silent = TRUE' suppresses direct plot output
  
  # Capture the plot as a grob
  plot_grob <- p$gtable
  
  return(plot_grob)
}

