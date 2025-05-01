
source("utility_functions_for_stats.R")
# change the cutoff as needed 


d1 <- read.csv("MOP_with_clustering.csv",
               row.names = 1, header = T, stringsAsFactors = F)

colnames(d1)

d1 <- d1[!duplicated(d1$PG.Genes), ]
dim(d1)
row.names(d1) <- d1$PG.Genes

d1 <- d1[, 8: dim(d1)[2]]



# one cluster vs all others
par(mfrow=c(2,3), tcl=-0.5, family="serif", mai=c(0.5,0.5,0.5,0.5))


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





