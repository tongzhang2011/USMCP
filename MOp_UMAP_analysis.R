
library(umap)
library(RColorBrewer)
library(gplots)
library(ggplot2)


############## final UMAP ################
d1 <- read.csv("MOP_with_clustering.csv",
               row.names = 1, header = T, stringsAsFactors = F)

colnames(d1)
d1 <- d1[, 8: dim(d1)[2]]

dim(d1)

d2 <- d1[complete.cases(d1),] # there are only 31 proteins in this case

d3 <- d1
d3[is.na(d3)] <- 0

# to perform U-MAP, each sample should be in a row

d3 <- t(d3)

umap_result <- umap(d3)

umap_df <- data.frame(UMAP1 = umap_result$layout[, 1],
                      UMAP2 = umap_result$layout[, 2],
                      Cluster = as.factor(  gsub(".*Cluster_", "", 
                                                 row.names(umap_result$layout) )))



custom_colors <- c(
  "1" = "#E41A1C",  # Red
  "2" = "#377EB8",  # Blue
  "3" = "#4DAF4A",  # Green
  "4" = "#FF7F00",  # Orange
  "5" = "#F781BF",  # Pink
  "6" = "#A65628",  # Brown
  "7" = "#984EA3"   # Purple
)

dev.off()

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 2) +
  labs(title = "",
       x = "UMAP 1", y = "UMAP 2") +
  scale_color_manual(values = custom_colors) +  # Apply the custom colors
  theme_minimal()+
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect( colour = "black",
                                 fill = NA,
                                 size = 1)
  )


################# ################## ######################   ######################   
############################# UMAP without cluter 2 and 6 ##################    
umap_df <- data.frame(UMAP1 = umap_result$layout[, 1],
                      UMAP2 = umap_result$layout[, 2],
                      Cluster = as.factor(  gsub(".*Cluster_", "", 
                                                 row.names(umap_result$layout) )))

library(dplyr)

umap_df <- filter(umap_df, Cluster != "2")
umap_df <- filter(umap_df, Cluster != "6")

table( umap_df$Cluster)

custom_colors <- c(
  "1" = "#E41A1C",  # Red
  
  "3" = "#4DAF4A",  # Green
  "4" = "#FF7F00",  # Orange
  "5" = "#F781BF",  # Pink
  
  "7" = "#984EA3"   # Purple
)



ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
  geom_point(size = 2) +
  labs(title = "",
       x = "UMAP 1", y = "UMAP 2") +
  scale_color_manual(values = custom_colors) +  # Apply the custom colors
  theme_minimal()+
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect( colour = "black",
                                 fill = NA,
                                 size = 1)
  )


################# ################## ######################   ######################   
############################# UMAP without cluter 2 and 6 ################## 

