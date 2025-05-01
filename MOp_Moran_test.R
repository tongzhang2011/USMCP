rm(list = ls())
library(dplyr)
library(pheatmap)
library(pacman)
p_load(knitr, dplyr, reshape2, spdep, raster, doParallel, scales, viridis, plot.matrix,
       data.table, dynamicTreeCut, RColorBrewer, ggplot2)
source("function_moran_test.R")

d1 <- read.csv("MOp_dataset.csv",
               row.names = 1, header = T, stringsAsFactors = F)
colnames(d1)
d1$Gene.names <- row.names(d1)

dim(d1)
hist(apply(d1, 1, function(x) sum(!is.na(x))))
d1 <- d1[apply(d1, 1, function(x) sum(!is.na(x))) > 493*0.75,]

which(row.names(d1) %in% "Nefl")


nrow <- 17
ncol <- 29
nPixels <- nrow * ncol

cores = 4


17*29

impute = FALSE
plot = FALSE
nsim = 999
lag.max = (nrow*ncol)/2
cutoff = 8

data.moran <- moran.Parallel_tz(d1,
                                d1$Gene.names,
                                nPixels = 493,
                                multiPlate = FALSE,
                                nrow = 17,
                                ncol = 29,
                                cores = 4
)

# write.csv(data.moran, "Manual_selected_voxel_Moran_analysis.csv")
write.csv(data.moran, "Moran_analysis_at_least_75pct.csv")

data.moran <- read.csv("Moran_analysis_at_least_75pct.csv", header = T, stringsAsFactors = F)
colnames(data.moran)
data.moran <- data.moran[, -1] 

data.moran$adj.p <- p.adjust(data.moran$p.value, method = "BH")
hist(data.moran$moran)
range(data.moran$moran)

dm <- filter(data.moran, adj.p < 0.01, moran > 0.4)

which(dm$gene %in% "Nefl")

dm <- dm[order(dm$moran, decreasing = TRUE),]
dim(dm)

mygene <- dm$gene[1:18]


library(dplyr)
library(stringr)
library(VennDiagram)
library(pheatmap)
require(RColorBrewer)

library(pheatmap)
library(gplots)
library(ggplot2)
library(gridExtra)

source("utility_function_for_individual_gene_plot.R")

d1 <- read.csv("MOp_dataset.csv",
               row.names = 1, header = T, stringsAsFactors = F)

grid.arrange(spatial_protein_plot_log2_scale_grid(mygene[1]),
             spatial_protein_plot_log2_scale_grid(mygene[2]),
             spatial_protein_plot_log2_scale_grid(mygene[3]),
             spatial_protein_plot_log2_scale_grid(mygene[4]),
             spatial_protein_plot_log2_scale_grid(mygene[5]),
             spatial_protein_plot_log2_scale_grid(mygene[6]),
             spatial_protein_plot_log2_scale_grid(mygene[7]),
             spatial_protein_plot_log2_scale_grid(mygene[8]),
             spatial_protein_plot_log2_scale_grid(mygene[9]),
             spatial_protein_plot_log2_scale_grid(mygene[10]),
             spatial_protein_plot_log2_scale_grid(mygene[11]),
             spatial_protein_plot_log2_scale_grid(mygene[12]),
             spatial_protein_plot_log2_scale_grid(mygene[13]),
             spatial_protein_plot_log2_scale_grid(mygene[14]),
             spatial_protein_plot_log2_scale_grid(mygene[15]),
             spatial_protein_plot_log2_scale_grid(mygene[16]),
             spatial_protein_plot_log2_scale_grid(mygene[17]),
             spatial_protein_plot_log2_scale_grid(mygene[18]),
             ncol = 6)


