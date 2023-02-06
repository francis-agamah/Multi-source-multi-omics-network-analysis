'''
statistical analysis of coexpression networks
'''
library(dplyr)
library(corrplot)
library("PerformanceAnalytics")

wkdir<-setwd("D:/PhD/Corona_Project")
prepped_datas<-file.path(wkdir,"datasets/prepared_data/")
trans1_all<-read.table(file.path(prepped_datas,"trans1_Over_coexp_network_all_healthy-mild-moderate-severe.txt"), header = T)
prot1_all<-file.path(prepped_datas,"prot1_Overmyer_coexp_network_all_healthy-mild-moderate-severe")


#using performance analytics
chart.Correlation(trans1_all[,3:ncol(trans1_all)], histogram=TRUE, pch=19)
plot(trans1_all)
