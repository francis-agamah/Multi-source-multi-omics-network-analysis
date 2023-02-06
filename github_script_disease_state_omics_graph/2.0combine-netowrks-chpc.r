library(dplyr)
library(tidyverse)

#this script is used to filter the combined network into various disease states

wkdir<-setwd("D:/PhD/Corona_Project")
prepped_datas<-file.path(wkdir,"datasets/prepared_data/")
unprepped_datas1<-file.path(prepped_datas,"unprepped_data/Overmyer/")
unprepped_datas2<-file.path(prepped_datas,"unprepped_data/Su/")

load("Trans1_coexpression.RData")
load("Trans2_coexpression.RData")
load("Prot1_coexpression.RData")
load("Prot2_coexpression.RData")
load("metabo1_coexpression.RData")
load("metabo2_coexpression.RData")
load("lipid1_coexpression.RData")

#all_trans_and_prot_disease_state<-rbind(pcc_trans1_disease_state,pcc_trans2_disease_state,pcc_prot1_disease_state,pcc_prot2_disease_state)
#write.table(all_trans_and_prot_disease_state,paste(prepped_datas,"All_prot_trans_and_metabo_all_healthy-mild-moderate-severe.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

#all_metabo_disease_state<- rbind(pcc_metabo1_disease_state,..)



################################################################################
#combing  interactions for disease state as per omics data type
################################################################################
combined_prot_network<-as.data.frame(rbind(pcc_prot1_disease_state,pcc_prot2_disease_state)) #for all cases
combined_prot_network$feature_type="protein-protein"
combined_prot_network$color = "blue"
write.table(combined_prot_network,paste(prepped_datas,"All_prot_all_healthy-mild-moderate-severe.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


combined_trans_network<-as.data.frame(rbind(pcc_trans1_disease_state,pcc_trans2_disease_state)) #for all cases
combined_trans_network$feature_type="transcript-transcript"
combined_trans_network$color="Yellow"
write.table(combined_trans_network,paste(prepped_datas,"All_trans_all_healthy-mild-moderate-severe.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

combined_metabolite_network<-as.data.frame(rbind(pcc_metabo1_disease_state,pcc_metabo2_disease_state)) #for all cases
combined_metabolite_network$feature_type="metabolite-metabolite"
combined_metabolite_network$color="orange"
write.table(combined_metabolite_network,paste(prepped_datas,"All_metabo_all_healthy-mild-moderate-severe.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

combined_lipid_network<-as.data.frame(pcc_lipid1_disease_state)
combined_lipid_network$feature_type="lipid-lipid"
combined_lipid_network$color="pink"
write.table(combined_lipid_network,paste(prepped_datas,"All_lipid_all_healthy-mild-moderate-severe.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

#for the purpouse of this rbind all networks
colnames(combined_metabolite_network)=c("GeneID1","GeneID2","PCC_healthy","p.value_healthy","pcc_direction_healthy",
                                        "PCC_mild","p.value_mild","pcc_direction_mild",
                                        "PCC_moderate","p.value_moderate","pcc_direction_moderate",
                                        "PCC_severe","p.value_severe","pcc_direction_severe",
                                        "PCC_case","p.value_case","pcc_direction_case","feature_type", "color")

#for the purpouse of this rbind all networks
colnames(combined_lipid_network)=c("GeneID1","GeneID2","PCC_healthy","p.value_healthy","pcc_direction_healthy",
                                        "PCC_mild","p.value_mild","pcc_direction_mild",
                                        "PCC_moderate","p.value_moderate","pcc_direction_moderate",
                                        "PCC_severe","p.value_severe","pcc_direction_severe",
                                        "PCC_case","p.value_case","pcc_direction_case","feature_type", "color")


all_trans_prot_lipid_and_metabo_disease_state<-rbind(combined_prot_network, combined_trans_network, combined_metabolite_network,combined_lipid_network)
write.table(all_trans_prot_lipid_and_metabo_disease_state,paste(prepped_datas,"All_prot_trans_lipid_and_metabo_all_healthy-mild-moderate-severe.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

#all_metabo_disease_state<- rbind(pcc_metabo1_disease_state,..
'''
generate plots to show the relationships between correlation coefficients
'''
par(mfrow = c(3, 3))
network<-read.table(file.path(prepped_datas,"All_prot_trans_lipid_and_metabo_all_healthy-mild-moderate-severe.txt"),sep='\t', header = T, row.names = NULL, stringsAsFactors = F,na.strings = T)
plot(network$PCC_mild, network$PCC_moderate, col='red', pch=10, cex=0.2,
     xlab='mild', ylab='moderate', main='mild vs moderate')

plot(network$PCC_mild, network$PCC_severe, col='red', pch=10, cex=0.2,
     xlab='mild', ylab='severe', main='mild vs severe')

plot(network$PCC_mild, network$PCC_healthy, col='red', pch=10, cex=0.2,
     xlab='mild', ylab='healthy', main='mild vs healthy')


plot(network$PCC_moderate, network$PCC_severe, col='red', pch=10, cex=0.2,
     xlab='moderate', ylab='severe', main='moderate vs severe')

plot(network$PCC_moderate, network$PCC_healthy, col='red', pch=10, cex=0.2,
     xlab='moderate', ylab='healthy', main='moderate vs healthy')

plot(network$PCC_severe, network$PCC_healthy, col='red', pch=10, cex=0.2,
     xlab='severe', ylab='healthy', main='severe vs healthy')

#plot(network$PCC_case, network$PCC_healthy, col='red', pch=10, cex=0.2,
#     xlab='case', ylab='healthy', main='healthy vs case')

