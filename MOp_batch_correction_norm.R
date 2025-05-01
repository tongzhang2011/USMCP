library(dplyr)
source("Utility_functions.R")

## the "MOp_Raw_intensities.csv" file is from the SN searching
d <- read.csv("MOp_Raw_intensities.csv", header = T, stringsAsFactors = F)
colnames(d)
# colnames(d) <- gsub("_.*", "", colnames(d))
# write.csv(d, "MOp_Raw_for_Supplemental_Data.csv", row.names = F)


d1 <- d[, 8:dim(d)[2]]

###### log2 transformation 
d1 <- log(d1,2)

############ batch correction ###############
table(grepl("_C62$", colnames(d1))) # 77
table(grepl("_C63$", colnames(d1))) # 234
table(grepl("_C64$", colnames(d1))) # 182

x1 <- select(d1, ends_with("C62"))
x1 <- global_median_norm_batch(x1)

x2 <- select(d1, ends_with("C63"))
x2 <- global_median_norm_batch(x2)


x3 <- select(d1, ends_with("C64"))
x3 <- global_median_norm_batch(x3)

d1 <- cbind(x1, x2, x3)
dim(d1)
############ batch correction ###############



###### normalization after batch correction
d1 <- global_median_norm(d1)
# boxplot(d1)
apply(d1, 2, function(x) median(x, na.rm = T))


# ensure that samples are in 1:493 order
sample_num <- gsub("_.*", "", colnames(d1))
sample_num <- gsub("R", "", sample_num)
sample_num <- gsub("\\.2", "", sample_num)

d1 <- d1[,match(1:493, sample_num)]  # reorder the samples
colnames(d1)

colnames(d1) <- gsub("_.*", "", colnames(d1))

write.csv(cbind(d[,1:7],d1), "MOp_norm_for_Supplemental_Data.csv", row.names = T)