#!/usr/bin/env Rscript

print("load libraries")
#Load libraries
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(stringi)
library(doParallel)
library(strict)
library(tools)

data <- fread("data") # needs to be more specific -- BED, BAM, CRAM? 
# e.g. suppose we have a wrapper function getFileType(bedOrBamOrCram) 
# further suppose that we call another function, say, getRangesFromFile(ibid)
# then we don't have to give a damn whether the input is bed, bam, or cram
# all we have to know is the number of counts for each contig and its sequence

# for example see CRAMtest.R against SCRAP runs

# in a 10X BAM this would have the edit-corrected UMI in a tag called "UB"
# the raw UMI barcode is in a tag "UR", its quality in a tag "UY", 
# and the same is true for the cell barcode (CB, CR, CY, as above)
#
# there's no particular reason not to adopt this strategy given its ubiquity; 
# we are even starting to use it for plate-based scRNA/scATAC/SCRAP merging.
# but in the case of the spike-ins, all of the below are available either from
# the BED/BAM/CRAM itself, or from annotations of the contigs in the reference
#
data_melt <- melt(data,
                  id.vars = c("GC", "UMI", "fragment_len", "CpG", "pos", "chr"),
                  measure.vars = "read_count")
data_melt$value <- as.numeric(as.character(data_melt$value))

data_melt$CpG_3 <- (data_melt$CpG) ^ (1 / 3)

data_melt$GC <- as.numeric(as.character(data_melt$GC))
data_melt$fragment_len <- as.numeric(data_melt$fragment_len)

##Source the GLM from the control GLM code
load("~/Projects/2020_Control_Project/2020_BatchAnalysis/Batch3_DT/2020_Gaussian_batch3.rda")

print("start fmol prediction")
##Predict fmol values
pred_conc <- as.data.frame(predict(gaussian_glm, data_melt))
names(pred_conc) <- "value"
##Replace value with conc
data_melt$value <- NULL
data_melt <- cbind(data_melt, pred_conc)
data_melt$variable <- NULL
colnames(data_melt) <- c("GC","UMI", "fragment_len", "CpG", "pos", "chr", "CpG3", "fmol")
data_melt <- data_melt[,c("chr", "pos", "fragment_len", "UMI", "CpG", "CpG3", "GC", "fmol")]

print(head(data_melt))

write.csv(data_melt, file = paste0("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch3_DT/", file_path_sans_ext(basename(filename)), "_fmol.csv"))
