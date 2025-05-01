library(dplyr)
library(stringr)
library(gplots)
library(pheatmap)
library(gridExtra)
library(grid)

source("function_heatmap_pathway.R")



# input 1: the enrichment result file 

# file after Re-normalization, using p = 0.01, FC = 0.5 as a cutoff
bp_file <- read.csv("Heatmap_pathway/Enrichment_Reactome_FCp5_both_up_and_down_protein_Re_norm.csv",
                    header = TRUE,
                    stringsAsFactors = FALSE)


bps <- unique(bp_file$Description)
length(bps) # in total of 247 



# input 2: the specific pathway to be plot bpi
# example pathway; "Neurotransmitter release cycle"
# bpi = "Neurotransmitter release cycle"


# input 3: the mapping data be to plot
d1 <- read.csv("MOP_with_clustering.csv",
               row.names = 1, header = T, stringsAsFactors = F)

colnames(d1)

d1 <- d1[!duplicated(d1$PG.Genes), ]
dim(d1)
row.names(d1) <- d1$PG.Genes

d1 <- d1[, 8: dim(d1)[2]]


colnames(d1)

d1_c1 <- d1 %>% dplyr::select(contains("Cluster_1"))
d1_c2 <- d1 %>% dplyr::select(contains("Cluster_2"))
d1_c3 <- d1 %>% dplyr::select(contains("Cluster_3"))
d1_c4 <- d1 %>% dplyr::select(contains("Cluster_4"))
d1_c5 <- d1 %>% dplyr::select(contains("Cluster_5"))
d1_c6 <- d1 %>% dplyr::select(contains("Cluster_6"))
d1_c7 <- d1 %>% dplyr::select(contains("Cluster_7"))

d1 <- cbind(d1_c1, d1_c2, d1_c3, d1_c4, d1_c5, d1_c6, d1_c7)
colnames(d1)


heatmap_of_pathway("GABA synthesis, release, reuptake and degradation")
heatmap_of_pathway("Neurotransmitter receptors and postsynaptic signal transmission")
dev.off()
heatmap_of_pathway("Neurotransmitter receptors and postsynaptic signal transmission",2)

heatmap_of_pathway("GABA synthesis, release, reuptake and degradation")
heatmap_of_pathway("The citric acid (TCA) cycle and respiratory electron transport")



