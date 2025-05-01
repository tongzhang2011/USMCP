library(dplyr)
source("Utility_functions.R")

d1 <- read.csv("QC_mouse_liver_dataset.csv", header = T, stringsAsFactors = F, row.names = 1)
colnames(d1)


at_least_cutoff <- 0.5
d1 <- d1[apply(d1, 1,function(x) sum(!is.na(x))) >= 20*at_least_cutoff, ]

# alternatively
d1 <- d1[complete.cases(d1),]


d1 <- log(d1, 2)
# boxplot(d1)


# data normalization
d3 <- d1
d3 <- global_median_norm(d3)
# boxplot(d3)  
d3 <- 2^d3  


############# cv distribution


cvs <- apply(d3, 1, function(x) cv(x))

table(is.na(cvs)) # note some NAs

cvs <- cvs[!is.na(cvs)]
# hist(cvs)


cvs <- cvs[ cvs < 500]

hist(cvs, 
     breaks = seq(0, 500, by = 5),
     xlab = expression(bold("Coefficient of variation (%)")),
     ylab = expression(bold("Frequency")),
     col = "gold",
     main = "")

abline(v = median(cvs), lty = 2, lwd =2 , col = "red")

text(200, 150, 
     paste0("Median = ", round(median(cvs), 2), "%"),
     col = "red")



######################## finish the plot #################


