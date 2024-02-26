setwd("~/Downloads/STOTEN Revision Data and Code/DEseq")
library(dplyr)
library(vegan)

FPM = read.csv("FPMcounts_alldata.csv")
design_table = read.csv("design_table.csv")
str(FPM)

FPM.mat = as.data.frame(x = t(FPM), stringsAsFactors = FALSE)

colnames(FPM.mat) <- unlist(FPM.mat[1, ])
FPM.mat <- FPM.mat[-1, ]
FPM.mat = as.data.frame(FPM.mat)

str(FPM.mat)
View(FPM.mat)

FPM.mat[] <- lapply(FPM.mat, as.numeric) 


adonis2(FPM.mat~age+orig_depth+age*orig_depth, 
        data=design_table, 
        permutations = 999, 
        method = "euclidean")

adonis2(FPM.mat~age*orig_depth*trans_depth, 
        data=design_table, 
        permutations = 999, 
        method = "euclidean")

