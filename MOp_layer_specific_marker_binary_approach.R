
d <- read.csv("MOp_Raw_intensities.csv", header = T, stringsAsFactors = F)

row.names(d) <- paste0(d$PG.Genes, "_", d$PG.ProteinGroups)


d1 <- d[, 8:dim(d)[2]]
colnames(d1)

# select voxels in the core of layers
L2_3 <- c(52:60, 69:77, 86:94, 103:111, 120:127, 140:144)
L5_6 <- c(171:176, 188:195, 205:212, 222:229, 239:245, 256:262, 273:279)
L6 <- c(307:308, 324:328, 341:347, 358:363, 375:380, 392:397, 409:414, 429:431)

L_all <- c(L2_3, L5_6, L6)

colnames(d1)
d2 <- d1[, 1:493 %in% L_all]
colnames(d2)

col_names_to_num <- function(x) {
  sample_num <- gsub("_.*", "",x)
  sample_num <- gsub("rR", "", sample_num)
  sample_num <- gsub("R", "", sample_num)
  sample_num <- gsub("\\.2", "", sample_num)
  
  return(sample_num)
}

my_samples <- col_names_to_num( colnames(d2))

d2a <- d2[,my_samples %in% L2_3]
d2b <- d2[,my_samples %in% L5_6]
d2c <- d2[,my_samples %in% L6]

colnames(d2a) <- paste0(colnames(d2a), "_L2")
colnames(d2b) <- paste0(colnames(d2b), "_L5")
colnames(d2c) <- paste0(colnames(d2c), "_L6")

d2 <- cbind(d2a, d2b, d2c)
colnames(d2)

a <- dim(d2a)[2]
b <- dim(d2b)[2]
c <- dim(d2c)[2]

# binary expression
# protein only expressed in layer L2

d3 <- d2
d3$valid_cnt_L2 <- apply(dplyr::select(d3, ends_with("L2")), 1, function(x) sum(!is.na(x)))
d3$valid_cnt_L5 <- apply(dplyr::select(d3, ends_with("L5")), 1, function(x) sum(!is.na(x)))
d3$valid_cnt_L6 <- apply(dplyr::select(d3, ends_with("L6")), 1, function(x) sum(!is.na(x)))

cut1 <- 0.75; cut2 <- 0.25
cut1 <- 0.60; cut2 <- 0.40



d3a <- filter(d3, valid_cnt_L2 > cut1*a, valid_cnt_L5 < cut2*b, valid_cnt_L6 < cut2*c)
d3b <- filter(d3, valid_cnt_L5 > cut1*b, valid_cnt_L2 < cut2*a, valid_cnt_L6 < cut2*c)
d3c <- filter(d3, valid_cnt_L6 > cut1*c, valid_cnt_L2 < cut2*a, valid_cnt_L5 < cut2*b)

write.csv(d3a, "Binary_expression_L2.csv")
write.csv(d3b, "Binary_expression_L5.csv")
write.csv(d3c, "Binary_expression_L6.csv")


################### the following code was not run #############
# working on dysrugulated proteins
LOI <- "L2"
d4 <- d2

higher_marker <- function(LOI) {
  
  d41 <- select(d4, ends_with(LOI))
  d42 <- select(d4, !ends_with(LOI))
  
  n1 <- dim(d41)[2]
  n2 <- dim(d42)[2]
  
  d4t <- cbind(d41, d42)
  group2_start <- n1 +1
  group2_end <- n1+n2
  
  d4t <- d4t[rowSums(!is.na(d4t)) > 50,]
  
  d4t$log2FC <- apply(d4t, 1, function(x) mean(x[1:n1], na.rm = T) - mean(x[group2_start : group2_end], na.rm = T))
  d4t$p.value <- apply(d4t, 1, function(x) t.test(x[1:n1], x[group2_start : group2_end], na.rm = T)$p.value)
  
  
  d4t <- filter(d4t, log2FC > 1, p.value  < 0.01)
  
  return(d4t)
}


d4a <- higher_marker("L2")
d4b <- higher_marker("L5")
d4c <- higher_marker("L6")


dev.off()

spatial_protein_plot_log2_scale("NFL_MOUSE")

# plot for d3a: ###############################
dim(d3a)

for(k in 1:1) {
  png(paste0("MoP_mapping_protein_Diff_set_", k, ".png" ),
      width = 4500, height = 3000, res = 300)
  
  
  start <- 1
  end <- 10
  
  
  heatmap_list <- vector("list", 10)
  
  for(i in start:end) {
    plot_tmp <- spatial_protein_plot_log2_scale(row.names(d3a)[i])$gtable
    
    # turn start number to 1, 2, 3, to 20
    heat_map_index <- i %% 20 # the last one is 0, need to turn to 20
    if( heat_map_index == 0) {heat_map_index <- 20}
    
    heatmap_list[[heat_map_index]] <- plot_tmp
    
  }
  
  grid.arrange(grobs = heatmap_list, nrow = 4, ncol = 5)
  
  dev.off()
}
# plot for d3a: ###############################


######### plot for d3b
for(k in 1:1) {
  png(paste0("MoP_mapping_protein_Diff_set_d3b_", k, ".png" ),
      width = 4500, height = 3000, res = 300)
  
  
  start <- 1
  end <- 4
  
  
  heatmap_list <- vector("list", 4)
  
  for(i in start:end) {
    plot_tmp <- spatial_protein_plot_log2_scale(row.names(d3b)[i])$gtable
    
    # turn start number to 1, 2, 3, to 20
    heat_map_index <- i %% 20 # the last one is 0, need to turn to 20
    if( heat_map_index == 0) {heat_map_index <- 20}
    
    heatmap_list[[heat_map_index]] <- plot_tmp
    
  }
  
  grid.arrange(grobs = heatmap_list, nrow = 4, ncol = 5)
  
  dev.off()
}
####################

dim(d3c)


for(k in 1:1) {
  png(paste0("MoP_mapping_protein_Diff_set_d3c_", k, ".png" ),
      width = 4500, height = 3000, res = 300)
  
  start <- k*20 - 19
  end <- k*20
  
  heatmap_list <- vector("list", 20)
  
  for(i in start:end) {
    plot_tmp <- spatial_protein_plot_log2_scale(row.names(d3c)[i])$gtable
    
    # turn start number to 1, 2, 3, to 20
    heat_map_index <- i %% 20 # the last one is 0, need to turn to 20
    if( heat_map_index == 0) {heat_map_index <- 20}
    
    heatmap_list[[heat_map_index]] <- plot_tmp
    
  }
  
  grid.arrange(grobs = heatmap_list, nrow = 4, ncol = 5)
  
  dev.off()
}




for(k in 1:1) {
  png(paste0("MoP_mapping_protein_Diff_set_d3c_2", k, ".png" ),
      width = 4500, height = 3000, res = 300)
  
  start <- 21
  end <- 33
  
  heatmap_list <- vector("list", 13)
  
  for(i in start:end) {
    plot_tmp <- spatial_protein_plot_log2_scale(row.names(d3c)[i])$gtable
    
    # turn start number to 1, 2, 3, to 20
    heat_map_index <- i %% 20 # the last one is 0, need to turn to 20
    if( heat_map_index == 0) {heat_map_index <- 20}
    
    heatmap_list[[heat_map_index]] <- plot_tmp
    
  }
  
  grid.arrange(grobs = heatmap_list, nrow = 4, ncol = 5)
  
  dev.off()
}





row.names(d3a)

for(k in 1:1) {
  png(paste0("MoP_mapping_protein_set_", k, ".png" ),
      width = 4500, height = 3000, res = 300)
  
  
  start <- k*20 - 19
  end <- k*20
  
  
  heatmap_list <- vector("list", 20)
  
  for(i in start:end) {
    plot_tmp <- spatial_protein_plot_log2_scale(row.names(d1)[i])$gtable
    
    # turn start number to 1, 2, 3, to 20
    heat_map_index <- i %% 20 # the last one is 0, need to turn to 20
    if( heat_map_index == 0) {heat_map_index <- 20}
    
    heatmap_list[[heat_map_index]] <- plot_tmp
    
  }
  
  grid.arrange(grobs = heatmap_list, nrow = 4, ncol = 5)
  
  dev.off()
}

