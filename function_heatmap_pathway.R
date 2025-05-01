library(dplyr)
library(stringr)
library(gplots)
library(pheatmap)
library(gridExtra)
library(grid)

# update function with customerized color coding ############################################
heatmap_of_pathway <- function(bpi, fold_range = 3) {
  
  bp <- bp_file %>% 
    filter(Description == bpi) %>% 
    arrange(desc(Count)) # take the one the large number of genes
  
  # extract proteins in the BP of interest
  gene_string <- as.character(bp[1,10])
  genes <- str_split(gene_string, "/")[[1]]
  
  
  
  # select rows (genes)
  d2 <- d1[genes,]
  
  pro.heat <- as.matrix(d2)
  
  
  color_range <- colorRampPalette(c("blue", "white", "red"))(100)
  breaks <- seq(-fold_range, fold_range, length.out = 101)  # Adjust range of color mapping here
  
  # to annote the columns
  column_annotation <- data.frame(
    Cluster = gsub(".*_Cluster_(\\d+)$", "\\1", colnames(pro.heat))
  )
  rownames(column_annotation) <- colnames(pro.heat)
  
  # Define custom colors for the clusters
  custom_colors <- c(
    "1" = "#E41A1C",  # Red
    "2" = "#377EB8",  # Blue
    "3" = "#4DAF4A",  # Green
    "4" = "#FF7F00",  # Orange
    "5" = "#F781BF",  # Pink
    "6" = "#A65628",  # Brown
    "7" = "#984EA3"   # Purple
  )
  
  
  p <- pheatmap(pro.heat, 
                scale = "row",    # Scales the rows (genes) to have mean = 0 and sd = 1
                clustering_distance_rows = "euclidean",  # Method for distance calculation in rows
                cluster_cols = FALSE,  # Do not cluster columns
                show_rownames = TRUE,  # Show row names (genes)
                show_colnames = FALSE, # Do not show column names (samples)
                color = color_range,  # Color gradient for heatmap
                breaks = breaks,  # Define the color scale range
                main = bpi,  # Title of the heatmap
                annotation_col = column_annotation,  # Add cluster annotations to columns
                annotation_colors = list(Cluster = custom_colors),  # Apply custom colors to the clusters
                silent = TRUE
  )
  
  return(p)
}
# update function with customerized color coding ############################################

