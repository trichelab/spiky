#!/usr/bin/env R

#Load libraries
library(rmarkdown)
library(tidyverse)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(boot)
library(BlandAltmanLeh)
library(scales)
library(stringr)
library(forcats)
library(varhandle)
library(purrr)
library(stringi)

data <- fread("data") # CSV currently, should become BED ? 

##calculate methylation specificty here to report in the glm summary
### methylated reads / total reads

##Remove methylated reads

##Add in molar amount information
##This will need to be hard set for "fragment_len" column not 10 column
## The molar amount info is in picomole
for (i in 1:nrow(data)){
  data[i,10] <- ifelse (data[i,2] == "160", 0.002, ifelse(data[i,2] == "80", 0.004, 0.001))
}

names(data)[10] <- "conc"

##Adjusting for the 0.01ng dilution
data$conc <- data$conc*0.9

##Cube root CpG number to return to normal distribution
data$CpG <- as.integer(data$CpG)
data$CpG_3 <- (data$CpG) ^ (1 / 3)

##Gaussian model
gaussian_glm <- glm(formula = conc ~ value + fragment_len + GC + CpG_3, data = data_melt, family = gaussian)

##Calculating Gaussian r2
r2_gaussian= 1 - (gaussian_glm$deviance / gaussian_glm$null.deviance)
r2_gaussian #  reprot in summary

summary(gaussian_glm) ##save summary
save(gaussian_glm, file = "2020_Gaussian_batch1.rda")

##Calculating Gaussian r2
r2_gaussian= 1-(gaussian_glm$deviance/gaussian_glm$null.deviance)
r2_gaussian #  reprot in summary


# refactor into separate function
##Bland-Altman plot
#
png("2020_BA_gaussian.png", height = 6, width = 6, units = "in", res = 300)
BA <- bland.altman.plot(gaussian_glm$data[,1], gaussian_glm$resid, conf.int = .95, col.points = gaussian_glm$data[,2], pch = 19, graph.sys = "ggplot2")

BA + aes(color = gaussian_glm$data[,8]) + 
  theme_bw()+
  scale_color_manual(values = c("black", "darkgrey", "navy"))+
  theme(legend.title = element_blank(), legend.position = "right", legend.text = element_text(size = 14), 
        axis.title.y = element_text(size = 20), axis.text.x  = element_text(size = 16),
        axis.text.y = element_text(vjust = 0.5, size = 16), axis.title.x = element_text(size = 20)) +
  xlab("Mean of measurements") +
  ylab("Difference")
dev.off()

