#! /usr/bin/Rscript
getwd()
setwd("/home/shey/Escritorio/mat_project/new_db")
library(ggplot2)
data <- read.csv("existence.csv")
data2 <- read.csv("freq_abs.csv")
data3 <- read.csv("abs_and_rel.csv")

#Set variables
dataList = c("data$Category","data2$Test","data3$Test")
for(item in dataList){
  item <- factor(item, levels = c("Disp", "No Disp", "Head", "Paleo"))
}

#Boxplot of existence
ggplot(data, aes(x=Category, y=Intra, color=Category)) +
  geom_boxplot(notch=TRUE)
ggsave("boxplot_existence.png", width = 5, height = 4)

#Boxplot of absolute frequencies
ylim1 = boxplot.stats(data2$Intra)$stats[c(1, 5)]
ggplot(data2, aes(x=Test, y=Intra, color=Test)) +
        geom_boxplot(outlier.shape = NA) + 
        coord_cartesian(ylim = ylim1*1.1)
ggsave("boxplot_abs_freq.png", width = 5, height = 4)


