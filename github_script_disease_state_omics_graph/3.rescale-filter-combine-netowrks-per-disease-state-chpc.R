library(dplyr)
library(tidyverse)

#this script is used to filter the combined network into various disease states
wkdir<-setwd("D:/PhD/Corona_Project")
prepped_datas<-file.path(wkdir,"datasets/prepared_data/")
unprepped_datas1<-file.path(prepped_datas,"unprepped_data/Overmyer/")
unprepped_datas2<-file.path(prepped_datas,"unprepped_data/Su/")

NetworkFile<-read.table(file.path(prepped_datas,"All_prot_trans_lipid_and_metabo_all_healthy-mild-moderate-severe.txt"),header = T,
                        sep = "\t", row.names = NULL,na.strings = T, stringsAsFactors = F, quote="", fill = F)

#all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% filter(abs(as.numeric(PCC_mild))>=0.2 & feature_type=="lipid-lipid")
#Take the abs PCC score, select id1-id1-pcc score-feature_type
###################################################################################
#PROTEIN-PROTEIN INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & feature_type=="protein-protein")
write.table(all_sig_mild,paste(prepped_datas,"Proteomics-all_sig_mild_absolute_PCC.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & feature_type=="protein-protein")
write.table(all_sig_moderate,paste(prepped_datas,"Proteomics-all_sig_moderate_absolute_PCC.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & feature_type=="protein-protein")
write.table(all_sig_severe,paste(prepped_datas,"Proteomics-all_sig_severe_absolute_PCC.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & feature_type=="protein-protein")
write.table(all_sig_healthy,paste(prepped_datas,"Proteomics-all_sig_healthy_absolute_PCC.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

#---------------------------------------------------------------------------------
###################################################################################
#TRANSCRIPT-TRANSCRIPT INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & feature_type=="transcript-transcript")
write.table(all_sig_mild,paste(prepped_datas,"Transcriptomics-all_sig_mild_absolute_PCC.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & feature_type=="transcript-transcript")
write.table(all_sig_moderate,paste(prepped_datas,"Transcriptomics-all_sig_moderate_absolute_PCC.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & feature_type=="transcript-transcript")
write.table(all_sig_severe,paste(prepped_datas,"Transcriptomics-all_sig_severe_absolute_PCC.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & feature_type=="transcript-transcript")
write.table(all_sig_healthy,paste(prepped_datas,"Transcriptomics-all_sig_healthy_absolute_PCC.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


###################################################################################
#METABOLITE-METABOLITE INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & feature_type=="metabolite-metabolite")
write.table(all_sig_mild,paste(prepped_datas,"Metabolomics-all_sig_mild_absolute_PCC.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & feature_type=="metabolite-metabolite")
write.table(all_sig_moderate,paste(prepped_datas,"Metabolomics-all_sig_moderate_absolute_PCC.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & feature_type=="metabolite-metabolite")
write.table(all_sig_severe,paste(prepped_datas,"Metabolomics-all_sig_severe_absolute_PCC.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & feature_type=="metabolite-metabolite")
write.table(all_sig_healthy,paste(prepped_datas,"Metabolomics-all_sig_healthy_absolute_PCC.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)



###################################################################################
#LIPID-LIPID INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & feature_type=="lipid-lipid")
write.table(all_sig_mild,paste(prepped_datas,"Lipidomics-all_sig_mild_absolute_PCC.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & feature_type=="lipid-lipid")
write.table(all_sig_moderate,paste(prepped_datas,"Lipidomics-all_sig_moderate_absolute_PCC.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & feature_type=="lipid-lipid")
write.table(all_sig_severe,paste(prepped_datas,"Lipidomics-all_sig_severe_absolute_PCC.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & feature_type=="lipid-lipid")
write.table(all_sig_healthy,paste(prepped_datas,"Lipidomics-all_sig_healthy_absolute_PCC.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###################################################################################
#PROTEIN-PROTEIN INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.1) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.1 & feature_type=="protein-protein")
write.table(all_sig_mild,paste(prepped_datas,"Proteomics-all_sig_mild_absolute_PCC_0.1.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.1 & feature_type=="protein-protein")
write.table(all_sig_moderate,paste(prepped_datas,"Proteomics-all_sig_moderate_absolute_PCC_0.1.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.1 & feature_type=="protein-protein")
write.table(all_sig_severe,paste(prepped_datas,"Proteomics-all_sig_severe_absolute_PCC_0.1.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.1 & feature_type=="protein-protein")
write.table(all_sig_healthy,paste(prepped_datas,"Proteomics-all_sig_healthy_absolute_PCC_0.1.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


###################################################################################
#TRANSCRIPT-TRANSCRIPT INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.1) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.1 & feature_type=="transcript-transcript")
write.table(all_sig_mild,paste(prepped_datas,"Transcriptomics-all_sig_mild_absolute_PCC_0.1.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.1 & feature_type=="transcript-transcript")
write.table(all_sig_moderate,paste(prepped_datas,"Transcriptomics-all_sig_moderate_absolute_PCC_0.1.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.1 & feature_type=="transcript-transcript")
write.table(all_sig_severe,paste(prepped_datas,"Transcriptomics-all_sig_severe_absolute_PCC_0.1.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.1 & feature_type=="transcript-transcript")
write.table(all_sig_healthy,paste(prepped_datas,"Transcriptomics-all_sig_healthy_absolute_PCC_0.1.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


###################################################################################
#METABOLITE-METABOLITE INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.1) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.1 & feature_type=="metabolite-metabolite")
write.table(all_sig_mild,paste(prepped_datas,"Metabolomics-all_sig_mild_absolute_PCC_0.1.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.1 & feature_type=="metabolite-metabolite")
write.table(all_sig_moderate,paste(prepped_datas,"Metabolomics-all_sig_moderate_absolute_PCC_0.1.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.1 & feature_type=="metabolite-metabolite")
write.table(all_sig_severe,paste(prepped_datas,"Metabolomics-all_sig_severe_absolute_PCC_0.1.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.1 & feature_type=="metabolite-metabolite")
write.table(all_sig_healthy,paste(prepped_datas,"Metabolomics-all_sig_healthy_absolute_PCC_0.1.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)



###################################################################################
#LIPID-LIPID INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.1 ) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.1 & feature_type=="lipid-lipid")
write.table(all_sig_mild,paste(prepped_datas,"Lipidomics-all_sig_mild_absolute_PCC_0.1.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.1 & feature_type=="lipid-lipid")
write.table(all_sig_moderate,paste(prepped_datas,"Lipidomics-all_sig_moderate_absolute_PCC_0.1.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.1 & feature_type=="lipid-lipid")
write.table(all_sig_severe,paste(prepped_datas,"Lipidomics-all_sig_severe_absolute_PCC_0.1.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.1 & feature_type=="lipid-lipid")
write.table(all_sig_healthy,paste(prepped_datas,"Lipidomics-all_sig_healthy_absolute_PCC_0.1.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)



#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###################################################################################
#PROTEIN-PROTEIN INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.2) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.2 & feature_type=="protein-protein")
write.table(all_sig_mild,paste(prepped_datas,"Proteomics-all_sig_mild_absolute_PCC_0.2.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.2 & feature_type=="protein-protein")
write.table(all_sig_moderate,paste(prepped_datas,"Proteomics-all_sig_moderate_absolute_PCC_0.2.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.2 & feature_type=="protein-protein")
write.table(all_sig_severe,paste(prepped_datas,"Proteomics-all_sig_severe_absolute_PCC_0.2.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.2 & feature_type=="protein-protein")
write.table(all_sig_healthy,paste(prepped_datas,"Proteomics-all_sig_healthy_absolute_PCC_0.2.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


###################################################################################
#TRANSCRIPT-TRANSCRIPT INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.2) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.2 & feature_type=="transcript-transcript")
write.table(all_sig_mild,paste(prepped_datas,"Transcriptomics-all_sig_mild_absolute_PCC_0.2.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.2 & feature_type=="transcript-transcript")
write.table(all_sig_moderate,paste(prepped_datas,"Transcriptomics-all_sig_moderate_absolute_PCC_0.2.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.2 & feature_type=="transcript-transcript")
write.table(all_sig_severe,paste(prepped_datas,"Transcriptomics-all_sig_severe_absolute_PCC_0.2.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.2 & feature_type=="transcript-transcript")
write.table(all_sig_healthy,paste(prepped_datas,"Transcriptomics-all_sig_healthy_absolute_PCC_0.2.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


###################################################################################
#METABOLITE-METABOLITE INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.2) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.2 & feature_type=="metabolite-metabolite")
write.table(all_sig_mild,paste(prepped_datas,"Metabolomics-all_sig_mild_absolute_PCC_0.2.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.2 & feature_type=="metabolite-metabolite")
write.table(all_sig_moderate,paste(prepped_datas,"Metabolomics-all_sig_moderate_absolute_PCC_0.2.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.2 & feature_type=="metabolite-metabolite")
write.table(all_sig_severe,paste(prepped_datas,"Metabolomics-all_sig_severe_absolute_PCC_0.2.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.2 & feature_type=="metabolite-metabolite")
write.table(all_sig_healthy,paste(prepped_datas,"Metabolomics-all_sig_healthy_absolute_PCC_0.2.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)



###################################################################################
#LIPID-LIPID INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.2 ) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.2 & feature_type=="lipid-lipid")
write.table(all_sig_mild,paste(prepped_datas,"Lipidomics-all_sig_mild_absolute_PCC_0.2.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.2 & feature_type=="lipid-lipid")
write.table(all_sig_moderate,paste(prepped_datas,"Lipidomics-all_sig_moderate_absolute_PCC_0.2.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.2 & feature_type=="lipid-lipid")
write.table(all_sig_severe,paste(prepped_datas,"Lipidomics-all_sig_severe_absolute_PCC_0.2.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.2 & feature_type=="lipid-lipid")
write.table(all_sig_healthy,paste(prepped_datas,"Lipidomics-all_sig_healthy_absolute_PCC_0.2.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)



#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###################################################################################
#PROTEIN-PROTEIN INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.3) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.3 & feature_type=="protein-protein")
write.table(all_sig_mild,paste(prepped_datas,"Proteomics-all_sig_mild_absolute_PCC_0.3.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.3 & feature_type=="protein-protein")
write.table(all_sig_moderate,paste(prepped_datas,"Proteomics-all_sig_moderate_absolute_PCC_0.3.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.3 & feature_type=="protein-protein")
write.table(all_sig_severe,paste(prepped_datas,"Proteomics-all_sig_severe_absolute_PCC_0.3.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.3 & feature_type=="protein-protein")
write.table(all_sig_healthy,paste(prepped_datas,"Proteomics-all_sig_healthy_absolute_PCC_0.3.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


###################################################################################
#TRANSCRIPT-TRANSCRIPT INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.3) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.3 & feature_type=="transcript-transcript")
write.table(all_sig_mild,paste(prepped_datas,"Transcriptomics-all_sig_mild_absolute_PCC_0.3.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.3 & feature_type=="transcript-transcript")
write.table(all_sig_moderate,paste(prepped_datas,"Transcriptomics-all_sig_moderate_absolute_PCC_0.3.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.3 & feature_type=="transcript-transcript")
write.table(all_sig_severe,paste(prepped_datas,"Transcriptomics-all_sig_severe_absolute_PCC_0.3.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.3 & feature_type=="transcript-transcript")
write.table(all_sig_healthy,paste(prepped_datas,"Transcriptomics-all_sig_healthy_absolute_PCC_0.3.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


###################################################################################
#METABOLITE-METABOLITE INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.3) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.3 & feature_type=="metabolite-metabolite")
write.table(all_sig_mild,paste(prepped_datas,"Metabolomics-all_sig_mild_absolute_PCC_0.3.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.3 & feature_type=="metabolite-metabolite")
write.table(all_sig_moderate,paste(prepped_datas,"Metabolomics-all_sig_moderate_absolute_PCC_0.3.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.3 & feature_type=="metabolite-metabolite")
write.table(all_sig_severe,paste(prepped_datas,"Metabolomics-all_sig_severe_absolute_PCC_0.3.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.3 & feature_type=="metabolite-metabolite")
write.table(all_sig_healthy,paste(prepped_datas,"Metabolomics-all_sig_healthy_absolute_PCC_0.3.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)



###################################################################################
#LIPID-LIPID INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.3 ) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.3 & feature_type=="lipid-lipid")
write.table(all_sig_mild,paste(prepped_datas,"Lipidomics-all_sig_mild_absolute_PCC_0.3.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.3 & feature_type=="lipid-lipid")
write.table(all_sig_moderate,paste(prepped_datas,"Lipidomics-all_sig_moderate_absolute_PCC_0.3.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.3 & feature_type=="lipid-lipid")
write.table(all_sig_severe,paste(prepped_datas,"Lipidomics-all_sig_severe_absolute_PCC_0.3.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.3 & feature_type=="lipid-lipid")
write.table(all_sig_healthy,paste(prepped_datas,"Lipidomics-all_sig_healthy_absolute_PCC_0.3.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)



#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###################################################################################
#PROTEIN-PROTEIN INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.4) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.4 & feature_type=="protein-protein")
write.table(all_sig_mild,paste(prepped_datas,"Proteomics-all_sig_mild_absolute_PCC_0.4.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.4 & feature_type=="protein-protein")
write.table(all_sig_moderate,paste(prepped_datas,"Proteomics-all_sig_moderate_absolute_PCC_0.4.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.4 & feature_type=="protein-protein")
write.table(all_sig_severe,paste(prepped_datas,"Proteomics-all_sig_severe_absolute_PCC_0.4.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.4 & feature_type=="protein-protein")
write.table(all_sig_healthy,paste(prepped_datas,"Proteomics-all_sig_healthy_absolute_PCC_0.4.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


###################################################################################
#TRANSCRIPT-TRANSCRIPT INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.4) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.4 & feature_type=="transcript-transcript")
write.table(all_sig_mild,paste(prepped_datas,"Transcriptomics-all_sig_mild_absolute_PCC_0.4.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.4 & feature_type=="transcript-transcript")
write.table(all_sig_moderate,paste(prepped_datas,"Transcriptomics-all_sig_moderate_absolute_PCC_0.4.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.4 & feature_type=="transcript-transcript")
write.table(all_sig_severe,paste(prepped_datas,"Transcriptomics-all_sig_severe_absolute_PCC_0.4.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.4 & feature_type=="transcript-transcript")
write.table(all_sig_healthy,paste(prepped_datas,"Transcriptomics-all_sig_healthy_absolute_PCC_0.4.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


###################################################################################
#METABOLITE-METABOLITE INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.4) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.4 & feature_type=="metabolite-metabolite")
write.table(all_sig_mild,paste(prepped_datas,"Metabolomics-all_sig_mild_absolute_PCC_0.4.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.4 & feature_type=="metabolite-metabolite")
write.table(all_sig_moderate,paste(prepped_datas,"Metabolomics-all_sig_moderate_absolute_PCC_0.4.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.4 & feature_type=="metabolite-metabolite")
write.table(all_sig_severe,paste(prepped_datas,"Metabolomics-all_sig_severe_absolute_PCC_0.4.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.4 & feature_type=="metabolite-metabolite")
write.table(all_sig_healthy,paste(prepped_datas,"Metabolomics-all_sig_healthy_absolute_PCC_0.4.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)



###################################################################################
#LIPID-LIPID INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.4 ) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.4 & feature_type=="lipid-lipid")
write.table(all_sig_mild,paste(prepped_datas,"Lipidomics-all_sig_mild_absolute_PCC_0.4.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.4 & feature_type=="lipid-lipid")
write.table(all_sig_moderate,paste(prepped_datas,"Lipidomics-all_sig_moderate_absolute_PCC_0.4.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.4 & feature_type=="lipid-lipid")
write.table(all_sig_severe,paste(prepped_datas,"Lipidomics-all_sig_severe_absolute_PCC_0.4.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.4 & feature_type=="lipid-lipid")
write.table(all_sig_healthy,paste(prepped_datas,"Lipidomics-all_sig_healthy_absolute_PCC_0.4.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###################################################################################
#PROTEIN-PROTEIN INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.5) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.5 & feature_type=="protein-protein")
write.table(all_sig_mild,paste(prepped_datas,"Proteomics-all_sig_mild_absolute_PCC_0.5.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.5 & feature_type=="protein-protein")
write.table(all_sig_moderate,paste(prepped_datas,"Proteomics-all_sig_moderate_absolute_PCC_0.5.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.5 & feature_type=="protein-protein")
write.table(all_sig_severe,paste(prepped_datas,"Proteomics-all_sig_severe_absolute_PCC_0.5.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.5 & feature_type=="protein-protein")
write.table(all_sig_healthy,paste(prepped_datas,"Proteomics-all_sig_healthy_absolute_PCC_0.5.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


###################################################################################
#TRANSCRIPT-TRANSCRIPT INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.5) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.5 & feature_type=="transcript-transcript")
write.table(all_sig_mild,paste(prepped_datas,"Transcriptomics-all_sig_mild_absolute_PCC_0.5.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.5 & feature_type=="transcript-transcript")
write.table(all_sig_moderate,paste(prepped_datas,"Transcriptomics-all_sig_moderate_absolute_PCC_0.5.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.5 & feature_type=="transcript-transcript")
write.table(all_sig_severe,paste(prepped_datas,"Transcriptomics-all_sig_severe_absolute_PCC_0.5.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.5 & feature_type=="transcript-transcript")
write.table(all_sig_healthy,paste(prepped_datas,"Transcriptomics-all_sig_healthy_absolute_PCC_0.5.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


###################################################################################
#METABOLITE-METABOLITE INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.5) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.5 & feature_type=="metabolite-metabolite")
write.table(all_sig_mild,paste(prepped_datas,"Metabolomics-all_sig_mild_absolute_PCC_0.5.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.5 & feature_type=="metabolite-metabolite")
write.table(all_sig_moderate,paste(prepped_datas,"Metabolomics-all_sig_moderate_absolute_PCC_0.5.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.5 & feature_type=="metabolite-metabolite")
write.table(all_sig_severe,paste(prepped_datas,"Metabolomics-all_sig_severe_absolute_PCC_0.5.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.5 & feature_type=="metabolite-metabolite")
write.table(all_sig_healthy,paste(prepped_datas,"Metabolomics-all_sig_healthy_absolute_PCC_0.5.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)



###################################################################################
#LIPID-LIPID INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.5 ) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.5 & feature_type=="lipid-lipid")
write.table(all_sig_mild,paste(prepped_datas,"Lipidomics-all_sig_mild_absolute_PCC_0.5.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.5 & feature_type=="lipid-lipid")
write.table(all_sig_moderate,paste(prepped_datas,"Lipidomics-all_sig_moderate_absolute_PCC_0.5.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.5 & feature_type=="lipid-lipid")
write.table(all_sig_severe,paste(prepped_datas,"Lipidomics-all_sig_severe_absolute_PCC_0.5.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.5 & feature_type=="lipid-lipid")
write.table(all_sig_healthy,paste(prepped_datas,"Lipidomics-all_sig_healthy_absolute_PCC_0.5.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###################################################################################
#PROTEIN-PROTEIN INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.6) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.6 & feature_type=="protein-protein")
write.table(all_sig_mild,paste(prepped_datas,"Proteomics-all_sig_mild_absolute_PCC_0.6.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.6 & feature_type=="protein-protein")
write.table(all_sig_moderate,paste(prepped_datas,"Proteomics-all_sig_moderate_absolute_PCC_0.6.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.6 & feature_type=="protein-protein")
write.table(all_sig_severe,paste(prepped_datas,"Proteomics-all_sig_severe_absolute_PCC_0.6.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.6 & feature_type=="protein-protein")
write.table(all_sig_healthy,paste(prepped_datas,"Proteomics-all_sig_healthy_absolute_PCC_0.6.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


###################################################################################
#TRANSCRIPT-TRANSCRIPT INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.6) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.6 & feature_type=="transcript-transcript")
write.table(all_sig_mild,paste(prepped_datas,"Transcriptomics-all_sig_mild_absolute_PCC_0.6.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.6 & feature_type=="transcript-transcript")
write.table(all_sig_moderate,paste(prepped_datas,"Transcriptomics-all_sig_moderate_absolute_PCC_0.6.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.6 & feature_type=="transcript-transcript")
write.table(all_sig_severe,paste(prepped_datas,"Transcriptomics-all_sig_severe_absolute_PCC_0.6.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.6 & feature_type=="transcript-transcript")
write.table(all_sig_healthy,paste(prepped_datas,"Transcriptomics-all_sig_healthy_absolute_PCC_0.6.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


###################################################################################
#METABOLITE-METABOLITE INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.6) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.6 & feature_type=="metabolite-metabolite")
write.table(all_sig_mild,paste(prepped_datas,"Metabolomics-all_sig_mild_absolute_PCC_0.6.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.6 & feature_type=="metabolite-metabolite")
write.table(all_sig_moderate,paste(prepped_datas,"Metabolomics-all_sig_moderate_absolute_PCC_0.6.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.6 & feature_type=="metabolite-metabolite")
write.table(all_sig_severe,paste(prepped_datas,"Metabolomics-all_sig_severe_absolute_PCC_0.6.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.6 & feature_type=="metabolite-metabolite")
write.table(all_sig_healthy,paste(prepped_datas,"Metabolomics-all_sig_healthy_absolute_PCC_0.6.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)



###################################################################################
#LIPID-LIPID INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.6 ) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.6 & feature_type=="lipid-lipid")
write.table(all_sig_mild,paste(prepped_datas,"Lipidomics-all_sig_mild_absolute_PCC_0.6.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.6 & feature_type=="lipid-lipid")
write.table(all_sig_moderate,paste(prepped_datas,"Lipidomics-all_sig_moderate_absolute_PCC_0.6.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.6 & feature_type=="lipid-lipid")
write.table(all_sig_severe,paste(prepped_datas,"Lipidomics-all_sig_severe_absolute_PCC_0.6.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.6 & feature_type=="lipid-lipid")
write.table(all_sig_healthy,paste(prepped_datas,"Lipidomics-all_sig_healthy_absolute_PCC_0.6.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###################################################################################
#PROTEIN-PROTEIN INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.7) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.7 & feature_type=="protein-protein")
write.table(all_sig_mild,paste(prepped_datas,"Proteomics-all_sig_mild_absolute_PCC_0.7.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.7 & feature_type=="protein-protein")
write.table(all_sig_moderate,paste(prepped_datas,"Proteomics-all_sig_moderate_absolute_PCC_0.7.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.7 & feature_type=="protein-protein")
write.table(all_sig_severe,paste(prepped_datas,"Proteomics-all_sig_severe_absolute_PCC_0.7.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.7 & feature_type=="protein-protein")
write.table(all_sig_healthy,paste(prepped_datas,"Proteomics-all_sig_healthy_absolute_PCC_0.7.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


###################################################################################
#TRANSCRIPT-TRANSCRIPT INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.7) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.7 & feature_type=="transcript-transcript")
write.table(all_sig_mild,paste(prepped_datas,"Transcriptomics-all_sig_mild_absolute_PCC_0.7.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.7 & feature_type=="transcript-transcript")
write.table(all_sig_moderate,paste(prepped_datas,"Transcriptomics-all_sig_moderate_absolute_PCC_0.7.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.7 & feature_type=="transcript-transcript")
write.table(all_sig_severe,paste(prepped_datas,"Transcriptomics-all_sig_severe_absolute_PCC_0.7.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.7 & feature_type=="transcript-transcript")
write.table(all_sig_healthy,paste(prepped_datas,"Transcriptomics-all_sig_healthy_absolute_PCC_0.7.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


###################################################################################
#METABOLITE-METABOLITE INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.7) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.7 & feature_type=="metabolite-metabolite")
write.table(all_sig_mild,paste(prepped_datas,"Metabolomics-all_sig_mild_absolute_PCC_0.7.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.7 & feature_type=="metabolite-metabolite")
write.table(all_sig_moderate,paste(prepped_datas,"Metabolomics-all_sig_moderate_absolute_PCC_0.7.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.7 & feature_type=="metabolite-metabolite")
write.table(all_sig_severe,paste(prepped_datas,"Metabolomics-all_sig_severe_absolute_PCC_0.7.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.7 & feature_type=="metabolite-metabolite")
write.table(all_sig_healthy,paste(prepped_datas,"Metabolomics-all_sig_healthy_absolute_PCC_0.7.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)



###################################################################################
#LIPID-LIPID INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.7 ) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.7 & feature_type=="lipid-lipid")
write.table(all_sig_mild,paste(prepped_datas,"Lipidomics-all_sig_mild_absolute_PCC_0.7.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.7 & feature_type=="lipid-lipid")
write.table(all_sig_moderate,paste(prepped_datas,"Lipidomics-all_sig_moderate_absolute_PCC_0.7.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.7 & feature_type=="lipid-lipid")
write.table(all_sig_severe,paste(prepped_datas,"Lipidomics-all_sig_severe_absolute_PCC_0.7.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.7 & feature_type=="lipid-lipid")
write.table(all_sig_healthy,paste(prepped_datas,"Lipidomics-all_sig_healthy_absolute_PCC_0.7.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###################################################################################
#PROTEIN-PROTEIN INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.9) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.9 & feature_type=="protein-protein")
write.table(all_sig_mild,paste(prepped_datas,"Proteomics-all_sig_mild_absolute_PCC_0.9.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.9 & feature_type=="protein-protein")
write.table(all_sig_moderate,paste(prepped_datas,"Proteomics-all_sig_moderate_absolute_PCC_0.9.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.9 & feature_type=="protein-protein")
write.table(all_sig_severe,paste(prepped_datas,"Proteomics-all_sig_severe_absolute_PCC_0.9.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.9 & feature_type=="protein-protein")
write.table(all_sig_healthy,paste(prepped_datas,"Proteomics-all_sig_healthy_absolute_PCC_0.9.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


###################################################################################
#TRANSCRIPT-TRANSCRIPT INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.8) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.8 & feature_type=="transcript-transcript")
write.table(all_sig_mild,paste(prepped_datas,"Transcriptomics-all_sig_mild_absolute_PCC_0.8.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.8 & feature_type=="transcript-transcript")
write.table(all_sig_moderate,paste(prepped_datas,"Transcriptomics-all_sig_moderate_absolute_PCC_0.8.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.8 & feature_type=="transcript-transcript")
write.table(all_sig_severe,paste(prepped_datas,"Transcriptomics-all_sig_severe_absolute_PCC_0.8.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.8 & feature_type=="transcript-transcript")
write.table(all_sig_healthy,paste(prepped_datas,"Transcriptomics-all_sig_healthy_absolute_PCC_0.8.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


###################################################################################
#METABOLITE-METABOLITE INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.8) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.8 & feature_type=="metabolite-metabolite")
write.table(all_sig_mild,paste(prepped_datas,"Metabolomics-all_sig_mild_absolute_PCC_0.8.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.8 & feature_type=="metabolite-metabolite")
write.table(all_sig_moderate,paste(prepped_datas,"Metabolomics-all_sig_moderate_absolute_PCC_0.8.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.8 & feature_type=="metabolite-metabolite")
write.table(all_sig_severe,paste(prepped_datas,"Metabolomics-all_sig_severe_absolute_PCC_0.8.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.8 & feature_type=="metabolite-metabolite")
write.table(all_sig_healthy,paste(prepped_datas,"Metabolomics-all_sig_healthy_absolute_PCC_0.8.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)



###################################################################################
#LIPID-LIPID INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.8 ) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.8 & feature_type=="lipid-lipid")
write.table(all_sig_mild,paste(prepped_datas,"Lipidomics-all_sig_mild_absolute_PCC_0.8.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.8 & feature_type=="lipid-lipid")
write.table(all_sig_moderate,paste(prepped_datas,"Lipidomics-all_sig_moderate_absolute_PCC_0.8.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.8 & feature_type=="lipid-lipid")
write.table(all_sig_severe,paste(prepped_datas,"Lipidomics-all_sig_severe_absolute_PCC_0.8.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.8 & feature_type=="lipid-lipid")
write.table(all_sig_healthy,paste(prepped_datas,"Lipidomics-all_sig_healthy_absolute_PCC_0.8.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###################################################################################
#PROTEIN-PROTEIN INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.9) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.9 & feature_type=="protein-protein")
write.table(all_sig_mild,paste(prepped_datas,"Proteomics-all_sig_mild_absolute_PCC_0.9.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.9 & feature_type=="protein-protein")
write.table(all_sig_moderate,paste(prepped_datas,"Proteomics-all_sig_moderate_absolute_PCC_0.9.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.9 & feature_type=="protein-protein")
write.table(all_sig_severe,paste(prepped_datas,"Proteomics-all_sig_severe_absolute_PCC_0.9.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.9 & feature_type=="protein-protein")
write.table(all_sig_healthy,paste(prepped_datas,"Proteomics-all_sig_healthy_absolute_PCC_0.9.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


###################################################################################
#TRANSCRIPT-TRANSCRIPT INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.9) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.9 & feature_type=="transcript-transcript")
write.table(all_sig_mild,paste(prepped_datas,"Transcriptomics-all_sig_mild_absolute_PCC_0.9.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.9 & feature_type=="transcript-transcript")
write.table(all_sig_moderate,paste(prepped_datas,"Transcriptomics-all_sig_moderate_absolute_PCC_0.9.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.9 & feature_type=="transcript-transcript")
write.table(all_sig_severe,paste(prepped_datas,"Transcriptomics-all_sig_severe_absolute_PCC_0.9.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.9 & feature_type=="transcript-transcript")
write.table(all_sig_healthy,paste(prepped_datas,"Transcriptomics-all_sig_healthy_absolute_PCC_0.9.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


###################################################################################
#METABOLITE-METABOLITE INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.9) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.9 & feature_type=="metabolite-metabolite")
write.table(all_sig_mild,paste(prepped_datas,"Metabolomics-all_sig_mild_absolute_PCC_0.9.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.9 & feature_type=="metabolite-metabolite")
write.table(all_sig_moderate,paste(prepped_datas,"Metabolomics-all_sig_moderate_absolute_PCC_0.9.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.9 & feature_type=="metabolite-metabolite")
write.table(all_sig_severe,paste(prepped_datas,"Metabolomics-all_sig_severe_absolute_PCC_0.9.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.9 & feature_type=="metabolite-metabolite")
write.table(all_sig_healthy,paste(prepped_datas,"Metabolomics-all_sig_healthy_absolute_PCC_0.9.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)



###################################################################################
#LIPID-LIPID INTERACTION
#we select only the node pairs and filter against the absoulute PCC value (PCC value !=0 & PCC value >=0.9 ) of the disease state
####################################################################################

all_sig_mild<-select(NetworkFile,c(GeneID1, GeneID2,PCC_mild,feature_type)) %>% mutate(PCC_mild_abs = abs(as.numeric(PCC_mild))) %>% filter(PCC_mild_abs!=0 & PCC_mild_abs!= "FALSE" & PCC_mild_abs >= 0.9 & feature_type=="lipid-lipid")
write.table(all_sig_mild,paste(prepped_datas,"Lipidomics-all_sig_mild_absolute_PCC_0.9.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_moderate<-select(NetworkFile,c(GeneID1, GeneID2,PCC_moderate,feature_type)) %>% mutate(PCC_moderate_abs = abs(as.numeric(PCC_moderate))) %>% filter(PCC_moderate_abs!=0 & PCC_moderate_abs!="FALSE" & PCC_moderate_abs >= 0.9 & feature_type=="lipid-lipid")
write.table(all_sig_moderate,paste(prepped_datas,"Lipidomics-all_sig_moderate_absolute_PCC_0.9.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_severe<-select(NetworkFile,c(GeneID1, GeneID2,PCC_severe,feature_type))  %>% mutate(PCC_severe_abs = abs(as.numeric(PCC_severe))) %>% filter(PCC_severe_abs!=0 & PCC_severe_abs!= "FALSE" & PCC_severe_abs >= 0.9 & feature_type=="lipid-lipid")
write.table(all_sig_severe,paste(prepped_datas,"Lipidomics-all_sig_severe_absolute_PCC_0.9.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

all_sig_healthy<-select(NetworkFile,c(GeneID1, GeneID2,PCC_healthy,feature_type)) %>% mutate(PCC_healthy_abs = abs(as.numeric(PCC_healthy))) %>% filter(PCC_healthy_abs!=0 & PCC_healthy_abs != "FALSE" & PCC_healthy_abs >= 0.9 & feature_type=="lipid-lipid")
write.table(all_sig_healthy,paste(prepped_datas,"Lipidomics-all_sig_healthy_absolute_PCC_0.9.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

