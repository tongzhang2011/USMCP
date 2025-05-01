rm(list = ls())
library(dplyr)
library(pheatmap)
library(pacman)
p_load(knitr, dplyr, reshape2, spdep, raster, doParallel, scales, viridis, plot.matrix,
       data.table, dynamicTreeCut, RColorBrewer, ggplot2)
source("Utility_functions.R")

d <- read.csv("MOp_norm_for_Supplemental_Data.csv", row.names = 1, header = T, stringsAsFactors = F)

row.names(d) <- d$PG.ProteinGroups
colnames(d)
d1 <- d[, 8:dim(d)[2]]


d1$Gene.names <- row.names(d1)
colnames(d1)

nrow <- 17
ncol <- 29
nPixels <- nrow * ncol


## Assign the voxel type

voxel.type <- rep('Tissue',493)
empty_index <- c(8:17, 30:34)
voxel.type[empty_index] <- "Empty"
table(voxel.type)

#### plot the voxel assignment table
voxel.type.m = matrix(voxel.type, nrow = 17, ncol = 29,
                      byrow = F)
par(mar = c(2,0,2,6))

plot(voxel.type.m, axis.col=NULL, axis.row=NULL, xlab='', ylab='', main = '', asp = T,
     col = c("#FFFFCC","#225EA8"), breaks = c("Tissue","Empty"))

dev.off()
###### voxel plot 

data <- d1
data.norm <- d1

tissue <- names(data.norm[,seq(1:493)[voxel.type == 'Tissue']])
tissue.index <- seq(1:493)[voxel.type == 'Tissue']
tissue.empty <- names(data.norm[,seq(1:493)[voxel.type == 'Empty']])
tissue.empty.index <- seq(1:493)[voxel.type == 'Empty']

data.tissue <- data[, voxel.type == "Tissue"]
dim(data.tissue)

data.tissue.t <- transpose(data.tissue)
rownames(data.tissue.t) <- colnames(data.tissue)
colnames(data.tissue.t) <- data.tissue.t["Gene.names", ]


# calculate the distance matrix 
dist.tissue <- dist(data.tissue.t[1:478,])  # selecte the top 478 rows after the transpose
hclust.tissue <- hclust(dist.tissue, method = "ward.D")
  
    dynamicClusters <- cutree(hclust.tissue, k = 7)
    
    temp <- c(rep(NaN,478))
    temp[tissue.index] <- dynamicClusters

# set up the plotting parameters    
    par(mar = c(2,2,2,2))
    
    tiff_filename <-  "MOp_Clustering.tiff"
    
    tiff(tiff_filename,
         width = 10, height = 8, units = "in", res = 100)
    
    plotClusters(temp,
                 title = "",
                 axis.col=NULL, axis.row=NULL,
                 na.col = "gray70", na.print = FALSE,
                 nrow = 17,
                 ncol = 29,
                 fmt.cell='%.0f', cex = 1.25,
                 key = NULL,
                 asp = T)
    
    dev.off()
##################################



####### output the clustering results for grid plotting 

clustering_result <- data.frame(sample_name = names(dynamicClusters),
                                Cluster_id = as.numeric(dynamicClusters))

write.csv(clustering_result, "MOp_clustering_Results_Re_norm.csv")


table(clustering_result$Cluster_id)


### re-extract the data without the gene annotations
d1 <- d[, 8:dim(d)[2]]

## note that clustering results is only for 478 tissue voxels
table(colnames(d1) %in% clustering_result$sample_name)
d1 <- d1[, colnames(d1) %in% clustering_result$sample_name]  # select tissue samples that have a cluster assigned
colnames(d1)

colnames(d1) <- paste0(colnames(d1), "_Cluster_", clustering_result$Cluster_id ) # adding the clustering ID into the sample names


# write a dataframe with the cluster ID in the colnames
x <- cbind(d[,1:7], d1)

write.csv(x, "MOP_with_clustering.csv",
          row.names = TRUE)


############ create another grid view of the clustering results

voxel_num <- colnames(d1)
voxel_num <- gsub("_.*", "", voxel_num)
voxel_num <- gsub("R", "", voxel_num)
voxel_num <- as.numeric(voxel_num)

cluster_num <- gsub(".*Cluster_", "", colnames(d1) )

# assign cluster 8 to all the empty voxels
empty_voxel <- setdiff(1:493, voxel_num)
empty_voxel_cluster <- rep(8, length(empty_voxel))



# build a dataframe for plotting
df <- data.frame(id = as.numeric(c(voxel_num, empty_voxel)),
                 cluster = as.factor(c(cluster_num, empty_voxel_cluster)))

df <- df[order(df$id), ]

df$row = rep(1:17, times = 29)
df$col = rep(1:29, each = 17)



custom_colors <- c(
  "1" = "#E41A1C",  # Red
  "2" = "#377EB8",  # Blue
  "3" = "#4DAF4A",  # Green
  "4" = "#FF7F00",  # Orange
  "5" = "#F781BF",  # Pink
  "6" = "#A65628",  # Brown
  "7" = "#984EA3",   # Purple
  "8" = "gray"
)

head(df, 20)

ggplot(df, aes(x = col, y = row, fill = factor(cluster))) +
  geom_tile(color = "white") +  # Add white borders between cells
  scale_fill_manual(values = custom_colors) +  # Apply custom colors
  labs(title = "Cluster Grid Plot", fill = "Cluster") +  # Set title and legend label
  theme_minimal() +  # Minimal theme
  theme(
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_blank(),  # Remove y-axis title
    axis.text.x = element_blank(),   # Remove x-axis labels
    axis.text.y = element_blank(),   # Remove y-axis labels
    axis.ticks = element_blank(),    # Remove axis ticks
    panel.grid = element_blank()     # Remove grid lines
  ) +
  coord_fixed(ratio = 1) +  # Ensure square tiles
  scale_y_reverse()  # Reverse the y-axis to correct the orientation


























# do not need to run this block again
### Start replacing the row.name with gene name 
########################################################################################################
library(dplyr)

d1 <- read.csv("Manual_Voxel_MOP_Mapping_SN_selected_493_samples_clean_norm.csv", row.names = 1, header = T, stringsAsFactors = F)
colnames(d1)

d1 <- global_median_norm(d1)  # in "utility_functions.R"


# note that the row.name is not gene name
# need to import original dataset and replace VPS29_MOUSE with acutal gene name
# od: original data

od <- read.table("20240802_221630_MOP_Mapping_666xSNE_ReportPG-nonorm.tsv", 
                 header = T, 
                 stringsAsFactors = F, 
                 sep = "\t",
                 quote = "") # this line is necessary to read all data; otherwise half rows will be dropped.
od <- dplyr::select(od, PG.Genes, PG.ProteinNames)
od$PG.ProteinNames <- gsub(";.*", "", od$PG.ProteinNames)
head(od)

table(row.names(d1) %in% od$PG.ProteinNames)
d1$PG.ProteinNames <- row.names(d1)
colnames(d1)
merged <- merge(d1, od, 
                by = "PG.ProteinNames", 
                all.x = TRUE)

merged <- filter(merged, PG.Genes != "")
merged$PG.Genes <- gsub(";.*", "", merged$PG.Genes)
merged <- merged[!duplicated(merged$PG.Genes), ]
row.names(merged) <- merged$PG.Genes
merged <- dplyr::select(merged, -c(PG.Genes,PG.ProteinNames))

colnames(merged)
write.csv(merged, 
          "Manual_Voxel_MOP_Mapping_SN_selected_493_samples_clean_norm_gene_name_Re_norm.csv",
          row.names = TRUE)

### Finished replacing the row.name with gene name 
########################################################################################################



############## differential analysis ################
d1 <- read.csv("Manual_Voxel_MOP_Mapping_SN_selected_493_samples_clean_norm_gene_name_Re_norm.csv",
               row.names = 1, header = T, stringsAsFactors = F)
colnames(d1)

clustering_result <- read.csv("Manual_Selected_Protein_ward.D2_Fixed_7_clustering_Results_Re_norm.csv",
                              header = T, stringsAsFactors = F)
table(clustering_result$Cluster_id)

table(colnames(d1) %in% clustering_result$sample_name)
d1 <- d1[, colnames(d1) %in% clustering_result$sample_name]  # select tissue samples that have a cluster assigned
colnames(d1)

colnames(d1) <- paste0(colnames(d1), "_Cluster_", clustering_result$Cluster_id ) # adding the clustering ID into the sample names

dim(d1)
write.csv(d1, "Manual_Voxel_MOP_Mapping_SN_selected_493_samples_clean_norm_gene_name_with_clustering_Re_norm.csv",
          row.names = TRUE)

source("utility_functions_for_stats.R")
# change the cutoff as needed 

dev.off()
# one cluster vs all others
par(mfrow=c(2,3), tcl=-0.5, family="serif", mai=c(0.5,0.5,0.5,0.5))

par(mfrow=c(1,6), mar = c(5.1, 3, 4.1, 2))

vplot_one_vs_all_other(d1, "Cluster_1")
vplot_one_vs_all_other(d1, "Cluster_2")
vplot_one_vs_all_other(d1, "Cluster_3")
vplot_one_vs_all_other(d1, "Cluster_4")
vplot_one_vs_all_other(d1, "Cluster_5")
vplot_one_vs_all_other(d1, "Cluster_7")

dev.off()



# output significant proteins
sig_list_one_vs_all_others(d1, "Cluster_1")
sig_list_one_vs_all_others(d1, "Cluster_2")
sig_list_one_vs_all_others(d1, "Cluster_3")
sig_list_one_vs_all_others(d1, "Cluster_4")
sig_list_one_vs_all_others(d1, "Cluster_5")
sig_list_one_vs_all_others(d1, "Cluster_7")



# move all the files into a folder "Sig_list_Re_norm_p001_FC1"
# move all the files into a folder "Sig_list_Re_norm_p001_FCp5"

######## enrichment analysis 
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(purrr)
library(stringr)
library(dplyr)
library(ReactomePA)


# Set the subfolder path (relative to the working directory)
subfolder <- "Sig_One_vs_Others_FC_1"  # Replace with your subfolder name
subfolder <- "Sig_One_vs_Others_FC_1_UP_protein_only"  # Replace with your subfolder name
subfolder <- "Sig_list_Re_norm_p001_FC1"  # Replace with your subfolder name
subfolder <- "Sig_list_Re_norm_p001_FCp5"  # Replace with your subfolder name


# List all .csv files in the subfolder
csv_files <- list.files(path = subfolder, pattern = "\\.csv$", full.names = TRUE)

# Read each .csv file into a separate dataframe and store in a list

csv_data_list <- map(csv_files, ~read.csv(.x, row.names = 1)) # make sure row.names setting


# Extract the filenames without extension and assign as names to the list
names(csv_data_list) <- gsub(".csv$", "", basename(csv_files))
names(csv_data_list) <- gsub("_vs_Others.*", "", names(csv_data_list))


# a function to turn to ENTREZID
# the input is dataframe
# ENTREZID is Gene ID used in NCBI. it is a number.
turn_to_ENTREZID <- function(data) {
  gene_symbols <- gsub(";.*", "", row.names(data))
  entrez_ids <- AnnotationDbi::select(org.Mm.eg.db, keys = gene_symbols, keytype = "SYMBOL", columns = "ENTREZID")
  entrez_ids <- entrez_ids$ENTREZID
  
  return(entrez_ids)
}

list_length <- length(csv_data_list)

gene_clusters <- vector("list", length = list_length)
for (i in 1:list_length) {
  gene_clusters[[i]] <- turn_to_ENTREZID(csv_data_list[[i]])
}



names(gene_clusters) <- names(csv_data_list)
sapply(gene_clusters, length) # note the numbers matches with the actual number



gene_clusters <- lapply(gene_clusters, function(x) x[!is.na(x)])

#################### the above is data data prep ######################


################ Gene ontology ##################### ##################
ck <- compareCluster(geneCluster = gene_clusters, fun = enrichGO,
                     OrgDb = org.Mm.eg.db, keyType="ENTREZID",
                     ont = "BP")

ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
head(ck) 

str(ck)
xx <- ck@compareClusterResult

# manually select item for plotting
dotplot(ck)


# filter by GeneRatio

enriched_data <- ck@compareClusterResult  # extract enrichment results

enriched_data <- enriched_data %>%
  mutate(
    GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)),  # Extract the numerator
    GeneRatio_den = as.numeric(sub(".*/", "", GeneRatio)),  # Extract the denominator
    GeneRatio_value = GeneRatio_num / GeneRatio_den         # Calculate the ratio
  )

# Now you can filter based on GeneRatio_value
enriched_data <- enriched_data %>%
  filter(GeneRatio_value > 0.1) %>%
  select(-c(GeneRatio_num,
            GeneRatio_den,
            GeneRatio_value, ))

ck@compareClusterResult <- enriched_data
dotplot(ck, font = 9)

############################ done with GO ###################################



################ KEGG pathway enrichment ##################### ##################
# perform KEGG pathway enrichment
ck <- compareCluster(geneCluster = gene_clusters, fun = enrichKEGG,
                     organism     = 'mmu')

enriched_data <- ck@compareClusterResult  # extract enrichment results
enriched_data$Description <- gsub(" - Mus musculus \\(house mouse\\)", "",
                                  enriched_data$Description)

ck@compareClusterResult <- enriched_data
dotplot(ck, font = 9)

# filter by GeneRatio

enriched_data <- ck@compareClusterResult  # extract enrichment results

enriched_data <- enriched_data %>%
  mutate(
    GeneRatio_num = as.numeric(sub("/.*", "", GeneRatio)),  # Extract the numerator
    GeneRatio_den = as.numeric(sub(".*/", "", GeneRatio)),  # Extract the denominator
    GeneRatio_value = GeneRatio_num / GeneRatio_den         # Calculate the ratio
  )

# Now you can filter based on GeneRatio_value
enriched_data <- enriched_data %>%
  filter(GeneRatio_value > 0.1) %>%
  select(-c(GeneRatio_num,
            GeneRatio_den,
            GeneRatio_value, ))

ck@compareClusterResult <- enriched_data
dotplot(ck, font = 9)
############################ done with KEGG ###################################



########################## Reactome analysis ################

######### use gene name, not number id in the result table 
ck1 <- compareCluster(geneCluster = gene_clusters, 
                      fun = enrichPathway,
                      organism     = 'mouse',
                      readable = TRUE)
enriched_data <- ck1@compareClusterResult  # extract enrichment results


write.csv(enriched_data, "Enrichment_Reactome_FCp5_both_up_and_down_protein_Re_norm.csv")
dotplot(ck1, font = 9)



# selected pathway for plotting
selected_pathways <- c(
  "The citric acid (TCA) cycle and respiratory electron transport",
  "Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC)",
  "Glutamate binding, activation of AMPA receptors and synaptic plasticity",
  "Eukaryotic Translation Initiation",
  "L13a−mediated translational silencing of Ceruloplasmin expression",
  "MAPK family signaling cascades",
  "Neurotransmitter receptors and postsynaptic signal transmission",
  "GABA synthesis, release, reuptake and degradation",
  "Dopamine Neurotransmitter Release Cycle",
  "Glutamate Neurotransmitter Release Cycle",
  "Ion transport by P−type ATPases",
  "Axon guidance",
  "Cellular responses to stress",
  "Cellular responses to stimuli",
  "Signaling by WNT",
  "RHO GTPases Activate WASPs and WAVEs",
  "Fcgamma receptor (FCGR) dependent phagocytosis",
  "Retrograde neurotrophin signalling")
length(selected_pathways)

selection_index <- enriched_data$Description %in% selected_pathways
table(selection_index)

enriched_data <- enriched_data[selection_index,]

# re-assign to the plotting object
ck1@compareClusterResult <- enriched_data

dotplot(ck1, font = 9, showCategory = 18)
?dotplot

############################ done with Reactome ###################################







############### old way ################

######## enrichment analysis 
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)

s1 <- read.csv("Cluster_1_vs_Others_sig_proteins_FC_1_Mapping.csv", header = TRUE, row.names = 1)
s2 <- read.csv("Cluster_2_vs_Others_sig_proteins_FC_1_Mapping.csv", header = TRUE, row.names = 1)
s3 <- read.csv("Cluster_3_vs_Others_sig_proteins_FC_1_Mapping.csv", header = TRUE, row.names = 1)
s4 <- read.csv("Cluster_4_vs_Others_sig_proteins_FC_1_Mapping.csv", header = TRUE, row.names = 1)
s5 <- read.csv("Cluster_5_vs_Others_sig_proteins_FC_1_Mapping.csv", header = TRUE, row.names = 1)
s7 <- read.csv("Cluster_7_vs_Others_sig_proteins_FC_1_Mapping.csv", header = TRUE, row.names = 1)

row.names(s1)

# convert to ENTREZID

turn_to_ENTREZID <- function(data) {
  gene_symbols <- gsub(";.*", "", row.names(data))
  entrez_ids <- select(org.Mm.eg.db, keys = gene_symbols, keytype = "SYMBOL", columns = "ENTREZID")
  entrez_ids <- entrez_ids$ENTREZID
  
  return(entrez_ids)
}

gene_clusters <- vector("list", length = 6)
gene_clusters[[1]] <- turn_to_ENTREZID(s1)
gene_clusters[[2]] <- turn_to_ENTREZID(s2)
gene_clusters[[3]] <- turn_to_ENTREZID(s3)
gene_clusters[[4]] <- turn_to_ENTREZID(s4)
gene_clusters[[5]] <- turn_to_ENTREZID(s5)
gene_clusters[[6]] <- turn_to_ENTREZID(s7)

gene_clusters <- as.list(gene_clusters)
names(gene_clusters) <- paste0("Cluser", c(1:5,7))

gene_clusters <- lapply(gene_clusters, function(x) x[!is.na(x)])

ck <- compareCluster(geneCluster = gene_clusters, fun = enrichGO,
                     OrgDb = org.Mm.eg.db, keyType="ENTREZID",
                     ont = "BP")

ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType="ENTREZID")

save(ck, file = "Mapping_comparing_cluster_BP.RData")

str(ck)
xx <- ck@compareClusterResult

# manually select item for plotting
dotplot(ck)


# do MF enrichment #####################################
ck <- compareCluster(geneCluster = gene_clusters, fun = enrichGO,
                     OrgDb = org.Mm.eg.db, keyType="ENTREZID",
                     ont = "MF")

ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
head(ck) 

str(ck)
xx <- ck@compareClusterResult

# xx$Description[xx$Description == "ligand-gated monoatomic ion channel activity involved in regulation of presynaptic membrane potential"] <- "ligand-gated ion channel"

ck@compareClusterResult <- xx
# manually select item for plotting
dotplot(ck)
# do MF enrichment #####################################



















############## use a for loop for plotting ##############
ROIs <- paste0("Cluster_", 1:7)

ROIs <- paste0("Cluster_", c(1:5,7)) # remove 6

par(mfrow=c(3,5), tcl=-0.5, family="serif", mai=c(0.25,0.25,0.25,0.25))

for (i in 1:6) {
  start_number <- i+1
  for (j in start_number:7) {
    vplot(d1, ROIs[i], ROIs[j])
  }
}

# p_value plot
for (i in 1:7) {
  start_number <- i+1
  for (j in start_number:7) {
    p_value_plot(d1, ROIs[i], ROIs[j])
  }
}
############## use a for loop for plotting ##############

dev.off()



# output significant proteins
for (i in 1:7) {
  start_number <- i+1
  for (j in start_number:7) {
    sig_list(d1, ROIs[i], ROIs[j])
  }
}





region1 <- "Cluster_1"
region2 <- "Cluster_2"
vplot(d1, "Cluster_1", "Cluster_2")


vplot <- function(d1,region1, region2) {
  ROI1 <- region1
  ROI2 <- region2
  
  my_title <- paste(ROI2, ROI1, sep = " / ")
  
  df1 <-  dplyr::select(d1, ends_with(ROI1))
  df2 <-  dplyr::select(d1, ends_with(ROI2))
  
  df <- cbind(df1, df2)
  
  df$p_value <- t.test_v2_df(df1, df2)
  df$FC <- FC(df1, df2)
  
  
  valid_p <- as.numeric(table(is.na(df$p_value))[1])
  my_title <- paste(my_title,  
                    paste0("(n = ", valid_p, ")"))
  
  p_cut <- 0.01
  FC_cut <- 1
  
  #  p_cut <- 0.01
  # FC_cut <- 0.5
  
  
  
  d_up <- filter(df, p_value < p_cut, FC > FC_cut)
  d_dn <- filter(df, p_value < p_cut, FC < -FC_cut)
  
  Increase <- dim(d_up)[1]
  Decrease <- dim(d_dn)[1]
  
  point_size <- 1
  
  # create a plot
  with(df,plot(FC, -log(p_value,10),pch = 20,
               xlab = "",
               ylab = "",
               cex.axis = 1, font.axis = 1,
               cex.lab = 0.8, font.lab = 1,
               cex = point_size))
  
  xmin <- par("usr")[1]
  xmax <- par("usr")[2]
  ymin <- par("usr")[3]
  ymax <- par("usr")[4]
  
  xrange <- xmax - xmin
  yrange <- ymax - ymin
  xlab_pos_1 <- xmin + 0.2*xrange
  xlab_pos_2 <- xmin + 0.8*xrange
  ylab_pos <- ymin + 0.9*yrange
  
  
  with(subset(df,p_value > p_cut),points(FC, -log(p_value,10),pch = 20, col = "gray",
                                         cex = point_size))
  with(subset(df,p_value < p_cut  & FC > FC_cut ),points(FC, -log(p_value,10),pch = 20, col = "red",
                                                         cex = point_size))
  with(subset(df,p_value < p_cut  & FC < -FC_cut ),points(FC, -log(p_value,10),pch = 20, col = "blue",
                                                          cex = point_size))
  abline(h = -log(p_cut,10),lty =2 ,col = "black")
  abline(v = FC_cut, lty =2 ,col = "red")
  abline(v = -FC_cut, lty =2 ,col = "blue")
  box(lwd = 2)
  text(xlab_pos_1, ylab_pos, paste("Decrease: ", Decrease, sep = ""),
       cex = 1.2, col = "blue")
  text(xlab_pos_2, ylab_pos, paste("Increase: ", Increase, sep = ""),
       cex = 1.2, col = "red")
  title(main = my_title)
  
}




p_value_plot <- function(d1,region1, region2) {
  ROI1 <- region1
  ROI2 <- region2
  
  my_title <- paste(ROI2, ROI1, sep = " / ")
  
  df1 <-  dplyr::select(d1, ends_with(ROI1))
  df2 <-  dplyr::select(d1, ends_with(ROI2))
  
  df <- cbind(df1, df2)
  
  df$p_value <- t.test_v2_df(df1, df2)
  # df$FC <- FC(df1, df2)
  
  
  valid_p <- as.numeric(table(is.na(df$p_value))[1])
  my_title <- paste(my_title,  
                    paste0("(n = ", valid_p, ")"))
  
  hist(df$p_value, 
       breaks = seq(0, 1, length.out = 100),
       col = c("blue", rep("gray", 99)),
       xlab = "",
       ylab = "",
       main = my_title)
  
}

dev.off()




sig_list <- function(d1, region1, region2) {
  ROI1 <- region1
  ROI2 <- region2
  
  my_title <- paste(ROI2, ROI1, sep = "_vs_")
  
  df1 <-  dplyr::select(d1, ends_with(ROI1))
  df2 <-  dplyr::select(d1, ends_with(ROI2))
  
  df <- cbind(df1, df2)
  
  df$p_value <- t.test_v2_df(df1, df2)
  df$FC <- FC(df1, df2)
  
  
  valid_p <- as.numeric(table(is.na(df$p_value))[1])
  my_title <- paste(my_title,  
                    paste0("(n = ", valid_p, ")"))
  
  p_cut <- 0.01
  FC_cut <- 1
  
  df_sig <- filter(df,p_value < p_cut)
  df_sig <- filter(df_sig,abs(FC) > FC_cut)
  
  write.csv(df_sig, paste(my_title, "sig_proteins_FC_1.csv", sep = "_"))
  
}

####################################################################################calcluate p and FC
# step 1: define t.test and FC for this case
###########################################
#a t.test function that requires at least 2 valid numbers from both groups
at_least_number <- 5

t.test_at_least_2_v2 <- function(num1, num2) {
  if( sum(!is.na(num1)) >= at_least_number &  sum(!is.na(num2)) >= at_least_number ) {
    my_p.value <- t.test(num1,num2,na.rm=TRUE,var.equal=TRUE)$p.value
  } else {  my_p.value <- NA
  }
  return( my_p.value)
}

#apply this funtion to two dataframe.
#df1 contains the first group;
#df2 contains the second group;
t.test_v2_df <- function(df1,df2) {
  n <- dim(df1)[1]
  my_p_value <- vector(length = n)
  
  for (i in 1:n) {
    my_p_value[i] <- t.test_at_least_2_v2(as.numeric(df1[i,]), as.numeric(df2[i,]))
  }
  return(my_p_value)
}
#end of the function
##########################################

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#similarly, FC also requires at least two 
#valid values from each comparing group
fc_at_least_2  <- function(num1, num2) {
  if( sum(!is.na(num1)) >= at_least_number &  sum(!is.na(num2)) >= at_least_number ) {
    fc <-  mean(num2, na.rm=TRUE) - mean(num1, na.rm=TRUE)
  } else { fc <- NA}
  return( fc)
}

#define a function use two dataframes as input
#i.e: do the calculation at the datafrmae level
FC <- function(df1, df2) {
  n <- dim(df1)[1]
  my_fc <- vector(length = n)
  
  for (i in 1:n) {
    my_fc[i] <- fc_at_least_2(as.numeric(df1[i,]), as.numeric(df2[i,]))
  }
  return(my_fc)
}
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

























