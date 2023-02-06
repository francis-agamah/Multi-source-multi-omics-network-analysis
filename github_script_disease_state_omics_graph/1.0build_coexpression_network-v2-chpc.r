library(dplyr)

wkdir<-setwd("D:/PhD/Corona_Project")
prepped_datas<-file.path(wkdir,"datasets/prepared_data/")
unprepped_datas1<-file.path(prepped_datas,"unprepped_data/Overmyer/")
unprepped_datas2<-file.path(prepped_datas,"unprepped_data/Su/")
covidKG_lit<-file.path(prepped_datas,"unprepped_data/covid19kg/covid19kg/export/sif")
###############################################################################################################################################
#LOAD THE ALREADY PREPARED DATA
#Once loaded, skip the steps and continue from INTEGRATING THE OMICS DATA AND UNIFIED GRAPH TO BUILD COEXPRESSION NETWORKS FOR EACH OMICS TYPE
###############################################################################################################################################

load("study-data.Rdata") # load all data and objects, then begin from pcc.
load("Trans1_coexpression.RData") # contains pcc_trans1_disease_state, pcc, pcc_healthy, pcc_mild, pcc_severe for Overmyer
load("Trans2_coexpression.RData") # Trans2 coexpression
load("Prot1_coexpression.RData")
load("Prot2_coexpression.RData")
load("metabo1_coexpression.RData")
load("metabo2_coexpression.RData")
load("lipid1_coexpression.RData")

#########################################################################################################################################################################
#LOADING THE PRE-PROCESSED OMICS MEASUREMENT FILES/DATA
#load prepared omics files
#########################################################################################################################################################################
trans1<-read.table(file.path(prepped_datas,"trans_omics_count_norm_transposed_Overmyer.tsv"),sep='\t', header = T, row.names = NULL, stringsAsFactors = F,na.strings = T)
trans2<-read.table(file.path(prepped_datas,"trans_omics_count_2_norm_transposed_Su.tsv"),sep='\t', header = T, row.names = NULL, stringsAsFactors = F,na.strings = T)

prot1<-read.table(file.path(prepped_datas,"prot_omics_count_norm_transposed_mapped_Overmyer.tsv"),sep='\t', header = T, row.names = NULL, stringsAsFactors = F,na.strings = T)
prot2<-read.table(file.path(prepped_datas,"prot_omics_count_2_norm_transposed_mapped_Su.tsv"),sep='\t', header = T, row.names = NULL)

metabo1<-read.table(file.path(prepped_datas,"metabo1_omics_count_rm_miss_norm.tsv"),sep='\t', header = T, row.names = NULL, stringsAsFactors = F,na.strings = T)
metabo2<-read.table(file.path(prepped_datas,"metabo2_omics_count_rm_miss_norm.tsv"),sep='\t', header = T, row.names = NULL, stringsAsFactors = F,na.strings = T)

lipido1<-read.table(file.path(prepped_datas,"lipid_omics_count_rm_miss_norm_Overmyer.tsv"),sep='\t', header = T, row.names = NULL, stringsAsFactors = F,na.strings = T)


####################################################################################################################################################################
#LOADING SAMPLE LABEL/METADATA TO SPLIT THE OMICS MEASUREMENT DATA ACCORDING TO THE VARIOUS DISEASE STATE
#load omics sample label
#header for Su sample label [condition	sample_id	Who_Ordinal_Scale	disease_state]
#Overmyer sample lable is the same for the other omics type [sample_id	Sample_label	condition	WHO_ordinal_at_day_28	Hospital_free_days_45	ICU_1	disease_state]
#####################################################################################################################################################################
trans1_sample_label<-read.table(file.path(prepped_datas,"Overmyer_transcriptomics_sample_label.txt"), header = T, sep = '\t')
prot1_sample_label<-trans1_sample_label
prot2_sample_label<-read.table(file.path(prepped_datas,"Su_proteomics_sample_label.txt"), header = T, sep = '\t') #sample label for Su sample differ. the difference is due to the additional healthy samples.
metabo2_sample_label<-read.table(file.path(prepped_datas,"Su_metabolomics_sample_label.txt"), header = T, sep = '\t')
#trans2_sample_label<-read.table(file.path(prepped_datas,"Su_transcriptomics_sample_label.txt"), header = T, sep = '\t')
trans2_sample_label<-read.table(file.path(prepped_datas,"Su_transcritomics_sample_description_plus_WOS.txt"), header = T, sep = '\t')

#extract samples label per disease state
disease_state_trans1<-factor(trans1_sample_label$disease_state, levels=c("healthy","mild","moderate","severe"))
trans1_sample_label_case<-trans1_sample_label %>% filter(condition=="COVID19")
trans1_sample_label_healthy<-trans1_sample_label %>% filter(condition!="COVID19")

disease_state_trans2<-factor(trans2_sample_label$disease_state, levels=c("healthy","mild","moderate","severe"))
trans2_sample_label_case<-trans2_sample_label %>% filter(condition=="COVID-19")
trans2_sample_label_healthy<-trans2_sample_label %>% filter(condition!="COVID-19")

disease_state_prot1<-factor(trans1_sample_label$disease_state, levels=c("healthy","mild","moderate","severe"))
prot1_sample_label_case<-trans1_sample_label %>% filter(condition=="COVID19")
prot1_sample_label_healthy<-trans1_sample_label %>% filter(condition!="COVID19")

disease_state_metabo1<-factor(trans1_sample_label$disease_state, levels=c("healthy","mild","moderate","severe"))
metabo1_sample_label_case<-trans1_sample_label %>% filter(condition=="COVID19")
metabo1_sample_label_healthy<-trans1_sample_label %>% filter(condition!="COVID19")

prot2_sample_label_case<-prot2_sample_label %>% filter(condition=="COVID19")
prot2_sample_label_healthy<-prot2_sample_label %>% filter(condition!="COVID19")

disease_state_metabo2<-factor(metabo2_sample_label$disease_state,levels = c("healthy","mild","moderate","severe"))
metabo2_sample_label_case<-metabo2_sample_label %>% filter(condition=="COVID19")
metabo2_sample_label_healthy<-metabo2_sample_label %>% filter(condition!="COVID19")

lipido1_sample_label_case<-trans1_sample_label %>% filter(condition=="COVID19")
lipido1_sample_label_healthy<-trans1_sample_label %>% filter(condition!="COVID19")

#extract omics data specific to samples for a condition (i.e, case or control)
`%notin%` <- Negate(`%in%`) # %notin% operator is not built-in and can be created by applying Negate function to %in%.

trans1_healthy<- cbind(trans1$Gene_ID, trans1[colnames(trans1)%in% filter(trans1_sample_label, disease_state =="healthy")$sample_id])
trans1_mild<- cbind(trans1$Gene_ID, trans1[colnames(trans1)%in% filter(trans1_sample_label, disease_state =="mild")$sample_id])
trans1_moderate<-cbind(trans1$Gene_ID, trans1[colnames(trans1)%in% filter(trans1_sample_label, disease_state =="moderate")$sample_id])
trans1_severe<-cbind(trans1$Gene_ID, trans1[colnames(trans1)%in% filter(trans1_sample_label, disease_state =="severe")$sample_id])
trans1_case<-cbind(trans1$Gene_ID ,trans1[(colnames(trans1)%in%trans1_sample_label_case$sample_id)])

prot1_healthy<- cbind(prot1$Gene_ID, prot1[colnames(prot1)%in% filter(trans1_sample_label, disease_state =="healthy")$sample_id])
prot1_mild<- cbind(prot1$Gene_ID, prot1[colnames(prot1)%in% filter(trans1_sample_label, disease_state =="mild")$sample_id])
prot1_moderate<-cbind(prot1$Gene_ID, prot1[colnames(prot1)%in% filter(trans1_sample_label, disease_state =="moderate")$sample_id])
prot1_severe<-cbind(prot1$Gene_ID, prot1[colnames(prot1)%in% filter(trans1_sample_label, disease_state =="severe")$sample_id])
prot1_case<-cbind(prot1$Gene_ID ,prot1[(colnames(prot1)%in%prot1_sample_label_case$sample_id)])

metabo1_healthy<- cbind(metabo1$metabolite_ID, metabo1[colnames(metabo1)%in% filter(trans1_sample_label, disease_state =="healthy")$sample_id])
metabo1_mild<- cbind(metabo1$metabolite_ID, metabo1[colnames(metabo1)%in% filter(trans1_sample_label, disease_state =="mild")$sample_id])
metabo1_moderate<-cbind(metabo1$metabolite_ID, metabo1[colnames(metabo1)%in% filter(trans1_sample_label, disease_state =="moderate")$sample_id])
metabo1_severe<-cbind(metabo1$metabolite_ID, metabo1[colnames(metabo1)%in% filter(trans1_sample_label, disease_state =="severe")$sample_id])
metabo1_case<-cbind(metabo1$metabolite_ID, metabo1[(colnames(metabo1)%in%metabo1_sample_label_case$sample_id)])

lipido1_healthy<- cbind(lipido1$lipid_ID, lipido1[colnames(lipido1)%in% filter(trans1_sample_label, disease_state =="healthy")$sample_id])
lipido1_mild<- cbind(lipido1$lipid_ID, lipido1[colnames(lipido1)%in% filter(trans1_sample_label, disease_state =="mild")$sample_id])
lipido1_moderate<-cbind(lipido1$lipid_ID, lipido1[colnames(lipido1)%in% filter(trans1_sample_label, disease_state =="moderate")$sample_id])
lipido1_severe<-cbind(lipido1$lipid_ID, lipido1[colnames(lipido1)%in% filter(trans1_sample_label, disease_state =="severe")$sample_id])
lipido1_case<-cbind(lipido1$lipid_ID ,lipido1[(colnames(lipido1)%in%trans1_sample_label_case$sample_id)])


#prot1_healthy<-prot1[(colnames(prot1)%notin%prot1_sample_label_case$sample_id)]
trans2_healthy<- cbind(trans2$Gene_ID, trans2[colnames(trans2)%in% filter(trans2_sample_label, disease_state =="healthy")$sample_id])
trans2_mild<- cbind(trans2$Gene_ID, trans2[colnames(trans2)%in% filter(trans2_sample_label, disease_state =="mild")$sample_id])
trans2_moderate<-cbind(trans2$Gene_ID, trans2[colnames(trans2)%in% filter(trans2_sample_label, disease_state =="moderate")$sample_id])
trans2_severe<-cbind(trans2$Gene_ID, trans2[colnames(trans2)%in% filter(trans2_sample_label, disease_state =="severe")$sample_id])
trans2_case<-cbind(trans2$Gene_ID ,trans2[(colnames(trans2)%in%trans2_sample_label_case$sample_id)])

prot2_healthy<- cbind(prot2$Gene_ID, prot2[colnames(prot2)%in% filter(prot2_sample_label, disease_state =="healthy")$sample_id])
prot2_mild<- cbind(prot2$Gene_ID, prot2[colnames(prot2)%in% filter(prot2_sample_label, disease_state =="mild")$sample_id])
prot2_moderate<-cbind(prot2$Gene_ID, prot2[colnames(prot2)%in% filter(prot2_sample_label, disease_state =="moderate")$sample_id])
prot2_severe<-cbind(prot2$Gene_ID, prot2[colnames(prot2)%in% filter(prot2_sample_label, disease_state =="severe")$sample_id])
prot2_case<-cbind(prot2$Gene_ID ,prot2[(colnames(prot2)%in%prot2_sample_label_case$sample_id)])

metabo2_healthy<- cbind(metabo2$metabolite_ID, metabo2[colnames(metabo2)%in% filter(metabo2_sample_label, disease_state =="healthy")$sample_id])
metabo2_mild<- cbind(metabo2$metabolite_ID, metabo2[colnames(metabo2)%in% filter(metabo2_sample_label, disease_state =="mild")$sample_id])
metabo2_moderate<-cbind(metabo2$metabolite_ID, metabo2[colnames(metabo2)%in% filter(metabo2_sample_label, disease_state =="moderate")$sample_id])
metabo2_severe<-cbind(metabo2$metabolite_ID, metabo2[colnames(metabo2)%in% filter(metabo2_sample_label, disease_state =="severe")$sample_id])
metabo2_case<-cbind(metabo2$metabolite_ID, metabo2[(colnames(metabo2)%in%metabo2_sample_label_case$sample_id)])


####################################################################################################################################################################
#BUILDING AN INTEGRATED KNOWLEDGE GRAPH
#covid knowledge graph
#####################################################################################################################################################################
covid_interactome<-read.table(file.path(prepped_datas,"genemania-interactions-1-1685.txt"),header = T, sep='\t')#may expand #put here plus covid knowlege graph
metabo_gene_interactome1<-read.table(file.path(prepped_datas,"Overmyer-metabolite-gene-network-KEGG.txt"),header = T, sep='\t')
metabo_gene_interactome2<-read.table(file.path(prepped_datas,"Su-metabolite-gene-network-KEGG.txt"),header = T, sep='\t')
metabo_gene_interactome<-rbind(metabo_gene_interactome1,metabo_gene_interactome2) #combining metabolite-gene interactome. the file has 4 columns |metabolite KEGG ID|Gene ID|Metabolite name as per the data|Gene name
covidKG<-read.table(file.path(covidKG_lit,"ExtractedCovidKnowledgeGraph.txt"),header = T, sep='\t',fill = T)
metabo1_metabo1_interactome<-read.table(file.path(prepped_datas,"metabo1_metabo1_interactome_all_healthy-mild-moderate-severe.txt"),header = T, sep='\t')
metabo1_metabo1_interactome<-filter(metabo1_metabo1_interactome, metabo1_metabo1_interactome$frequency>=0.7 & metabo2_metabo2_interactome$frequency!=1.00000000)
metabo2_metabo2_interactome<-read.table(file.path(prepped_datas,"metabo2_metabo2_interactome_all_healthy-mild-moderate-severe.txt"),header = T, sep='\t',fill = T)
metabo2_metabo2_interactome<-filter(metabo2_metabo2_interactome, metabo2_metabo2_interactome$frequency>=0.7 & metabo2_metabo2_interactome$frequency!=1.00000000)
lipido1_lipido1_interactome<-read.table(file.path(prepped_datas,"lipido1_lipido1_interactome_all_healthy-mild-moderate-severe.txt"),header = T, sep='\t',fill = T)
lipido1_lipido1_interactome<-filter(lipido1_lipido1_interactome, lipido1_lipido1_interactome$frequency>=0.7 & lipido1_lipido1_interactome$frequency!=1.00000000)


#intra-omics network
trans1_geneList<-trans1[,1]
trans2_geneList<-trans2[,1]
prot1_geneList<-prot1[,1]
prot2_geneList<-prot2[,1]
#combined_prot_trans_List<-cbind(trans1_geneList,trans2_geneList,prot1_geneList,prot2_geneList) #combining the protein and gene list. Easy to check metabolite-gene or metabolite-protein interaction
metabolite1<-read.table(file.path(prepped_datas,"Overmyer-metabolites.txt"), header = T, sep = '\t') #reading mapped metabolites. The file has 4 columns representing the metabolite in different IDs. |metabolite name|KEGG ID| PubChem ID| CHEBI
metabolite2<-read.table(file.path(prepped_datas,"Su-metabolites.txt"), header = T, sep = '\t') #reading mapped metabolites

#inter-omics network


covid_interactome_in_trans1 = covid_interactome[(covid_interactome[,1]%in%trans1_geneList) & (covid_interactome[,2]%in%trans1_geneList),  ] #filtering PPI to maintain only those in the interactome
#covid_interactome_in_trans1=covid_interactome_in_trans1[1:500, ] #for testing purposes
covid_interactome_in_trans2 = covid_interactome[(covid_interactome[,1]%in%trans2_geneList) & (covid_interactome[,2]%in%trans2_geneList),  ]
covid_interactome_in_prot1 = covid_interactome[(covid_interactome[,1]%in%prot1_geneList) & (covid_interactome[,2]%in%prot1_geneList),  ]
covid_interactome_in_prot2 = covid_interactome[(covid_interactome[,1]%in%prot2_geneList) & (covid_interactome[,2]%in%prot2_geneList),  ]
#covid_interactome_in_metabo1 = metabo_gene_interactome[(metabo_gene_interactome[,1]%in%metabolite1$KEGG.ID) & (metabo_gene_interactome[,4]%in%combined_prot_trans_List), ] #metabolite-gene or metabolite-protein interaction present across the metabolite and proteomics and transcriptomics data
#covid_interactome_in_metabo2 = metabo_gene_interactome[(metabo_gene_interactome[,1]%in%metabolite2$KEGG.ID) & (metabo_gene_interactome[,4]%in%combined_prot_trans_List), ]


colnames(trans1_case)[1]="Gene_ID"; colnames(trans1_mild)[1]="Gene_ID"; colnames(trans1_moderate)[1]="Gene_ID"; colnames(trans1_severe)[1]="Gene_ID"; colnames(trans1_healthy)[1]="Gene_ID";
colnames(prot1_case)[1]="Gene_ID"; colnames(prot1_mild)[1]="Gene_ID"; colnames(prot1_moderate)[1]="Gene_ID"; colnames(prot1_severe)[1]="Gene_ID"; colnames(prot1_healthy)[1]="Gene_ID";
colnames(metabo1_case)[1]="metabolite_ID"; colnames(metabo1_mild)[1]="metabolite_ID"; colnames(metabo1_moderate)[1]="metabolite_ID"; colnames(metabo1_severe)[1]="metabolite_ID"; colnames(metabo1_healthy)[1]="metabolite_ID";
colnames(trans2_case)[1]="Gene_ID"; colnames(trans2_mild)[1]="Gene_ID"; colnames(trans2_moderate)[1]="Gene_ID"; colnames(trans2_severe)[1]="Gene_ID"; colnames(trans2_healthy)[1]="Gene_ID";
colnames(prot2_case)[1]="Gene_ID"; colnames(prot2_mild)[1]="Gene_ID"; colnames(prot2_moderate)[1]="Gene_ID"; colnames(prot2_severe)[1]="Gene_ID"; colnames(prot2_healthy)[1]="Gene_ID";
colnames(metabo2_case)[1]="metabolite_ID"; colnames(metabo2_mild)[1]="metabolite_ID"; colnames(metabo2_moderate)[1]="metabolite_ID"; colnames(metabo2_severe)[1]="metabolite_ID"; colnames(metabo2_healthy)[1]="metabolite_ID";


save(trans1,trans2, prot1, prot2, trans1_sample_label, trans2_sample_label, prot2_sample_label, metabo2_sample_label, disease_state_trans1,trans1_sample_label_case,
     trans1_sample_label_healthy, disease_state_trans2,trans2_sample_label_case, trans2_sample_label_healthy, disease_state_prot1, prot1_sample_label_case,prot1_sample_label_healthy, disease_state_metabo1, metabo1_sample_label_case,
     metabo1_sample_label_healthy, prot2_sample_label_case, prot2_sample_label_healthy,disease_state_metabo2,metabo2_sample_label_case,
     metabo2_sample_label_healthy, trans1_healthy,trans1_mild, trans1_moderate,trans1_severe,trans1_case, prot1_healthy,
     prot1_mild,prot1_moderate,prot1_severe,prot1_case,metabo1_healthy,metabo1_mild,metabo1_moderate,metabo1_severe, metabo1, metabo2,
     metabo1_case,trans2_healthy, trans2_mild,trans2_moderate, trans2_severe,trans2_case,prot2_healthy,prot2_mild, prot2_moderate,prot2_severe,
     prot2_case, metabo2_healthy,metabo2_mild,metabo2_moderate,metabo2_severe, metabo2_case, covid_interactome,metabo_gene_interactome1, metabo_gene_interactome2,
     metabo_gene_interactome,trans1_geneList,trans2_geneList,prot1_geneList,prot2_geneList,combined_prot_trans_List, metabolite1,metabolite2, covid_interactome_in_trans1,covid_interactome_in_trans2,
     covid_interactome_in_prot1,covid_interactome_in_prot2,covid_interactome_in_metabo1, covid_interactome_in_metabo2,metabo1_metabo1_interactome,metabo2_metabo2_interactome,lipido1,
     covidKG, lipido1_sample_label_case, lipido1_sample_label_healthy,lipido1_healthy, lipido1_mild, lipido1_moderate,lipido1_case,lipido1_severe,lipido1_lipido1_interactome, file="study-data.Rdata")




##################################################################################################
#INTEGRATING THE OMICS DATA AND UNIFIED GRAPH TO BUILD COEXPRESSION NETWORKS FOR EACH OMICS TYPE
#################################################################################################
#for trans1
pcc = c() #for case
pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
for (i in 1:nrow(covid_interactome_in_trans1)){
  #print(i)
  temp_pcc = c() # this is for case (combined disease state except healthy)
  temp_pcc_healthy = c() #this is for healthy
  temp_pcc_mild = c(); temp_pcc_moderate = c(); temp_pcc_severe = c()
  
  temp_gene1 = covid_interactome_in_trans1[i,1]
  temp_gene1_row<-which(trans1$Gene_ID == temp_gene1)
  print(temp_gene1_row)

  temp_gene2 = covid_interactome_in_trans1[i,2]
  temp_gene2_row<-which(trans1$Gene_ID == temp_gene2)
  weight_temp_gene1<-covid_interactome_in_trans1[i,3]
  
  temp_exp1 = trans1_case[ temp_gene1_row ,2:ncol(trans1_case)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp2 = trans1_case[ temp_gene2_row ,2:ncol(trans1_case)]
  temp_pcc_test = cor.test( as.numeric(temp_exp1), as.numeric(temp_exp2),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_test
  
  temp_exp3 = trans1_healthy[ temp_gene1_row ,2:ncol(trans1_healthy)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp4 = trans1_healthy[ temp_gene2_row ,2:ncol(trans1_healthy)]
  temp_pcc_healthy_test = cor.test( as.numeric(temp_exp3), as.numeric(temp_exp4),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_healthy_test

  temp_exp5 = trans1_mild[ temp_gene1_row ,2:ncol(trans1_mild)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp6 = trans1_mild[ temp_gene2_row ,2:ncol(trans1_mild)]
  temp_pcc_mild_test = cor.test( as.numeric(temp_exp5), as.numeric(temp_exp6),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_mild_test
  
  temp_exp7 = trans1_moderate[ temp_gene1_row ,2:ncol(trans1_moderate)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp8 = trans1_moderate[ temp_gene2_row ,2:ncol(trans1_moderate)]
  temp_pcc_moderate_test = cor.test( as.numeric(temp_exp7), as.numeric(temp_exp8),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_moderate_test

  temp_exp9 = trans1_severe[ temp_gene1_row ,2:ncol(trans1_severe)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp10 = trans1_severe[ temp_gene2_row ,2:ncol(trans1_severe)]
  temp_pcc_severe_test = cor.test( as.numeric(temp_exp9), as.numeric(temp_exp10),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_severe_test
  
  temp_pcc[1] = temp_gene1
  temp_pcc[2] = temp_gene2
  temp_pcc[3] = temp_pcc_test$estimate
  temp_pcc[4] = temp_pcc_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if( temp_pcc[3]>0 ){
    temp_pcc[5] = 1 #PCC >0
  }else{
    temp_pcc[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc = as.data.frame(rbind(pcc,temp_pcc ))
  
  temp_pcc_healthy[1] = temp_gene1
  temp_pcc_healthy[2] = temp_gene2
  temp_pcc_healthy[3] = temp_pcc_healthy_test$estimate
  temp_pcc_healthy[4] = temp_pcc_healthy_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if( temp_pcc_healthy[3]>0 ){
    temp_pcc_healthy[5] = 1 #PCC >0
  }else{
    temp_pcc_healthy[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc_healthy = as.data.frame(rbind(pcc_healthy,temp_pcc_healthy ))
  
  temp_pcc_mild[1] = temp_gene1
  temp_pcc_mild[2] = temp_gene2
  temp_pcc_mild[3] = temp_pcc_mild_test$estimate
  temp_pcc_mild[4] = temp_pcc_mild_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if( temp_pcc_mild[3]>0 ){
    temp_pcc_mild[5] = 1 #PCC >0
  }else{
    temp_pcc_mild[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc_mild = as.data.frame(rbind(pcc_mild,temp_pcc_mild))
  
  temp_pcc_moderate[1] = temp_gene1
  temp_pcc_moderate[2] = temp_gene2
  temp_pcc_moderate[3] = temp_pcc_moderate_test$estimate
  temp_pcc_moderate[4] = temp_pcc_moderate_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if( temp_pcc_moderate[3]>0 ){
    temp_pcc_moderate[5] = 1 #PCC >0
  }else{
    temp_pcc_moderate[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc_moderate = as.data.frame(rbind(pcc_moderate,temp_pcc_moderate))
  
  temp_pcc_severe[1] = temp_gene1
  temp_pcc_severe[2] = temp_gene2
  temp_pcc_severe[3] = temp_pcc_severe_test$estimate
  temp_pcc_severe[4] = temp_pcc_severe_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if( temp_pcc_severe[3]>0 ){
    temp_pcc_severe[5] = 1 #PCC >0
  }else{
    temp_pcc_severe[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc_severe = as.data.frame(rbind(pcc_severe,temp_pcc_severe))
}


rm( list=ls(pat="temp_") )
rownames(pcc)=NULL
colnames(pcc)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

rownames(pcc_healthy)=NULL
colnames(pcc_healthy)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

rownames(pcc_mild)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc_mild)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

rownames(pcc_moderate)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc_moderate)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

rownames(pcc_severe)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc_severe)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

pcc_trans1_disease_state<-mutate(pcc_healthy, pcc_mild$PCC,pcc_mild$p.value,pcc_mild$pcc_direction,pcc_moderate$PCC,pcc_moderate$p.value,pcc_moderate$pcc_direction,
                                 pcc_severe$PCC, pcc_severe$p.value, pcc_severe$pcc_direction, pcc$PCC,pcc$p.value,pcc$pcc_direction) #pcc across the various disease state 

colnames(pcc_trans1_disease_state)=c("GeneID1","GeneID2","PCC_healthy","p.value_healthy","pcc_direction_healthy",
                                     "PCC_mild","p.value_mild","pcc_direction_mild",
                                     "PCC_moderate","p.value_moderate","pcc_direction_moderate",
                                     "PCC_severe","p.value_severe","pcc_direction_severe",
                                     "PCC_case","p.value_case","pcc_direction_case")

save(pcc_trans1_disease_state, pcc, pcc_healthy, pcc_mild, pcc_severe,file = "Trans1_coexpression.RData")
write.table(pcc_trans1_disease_state,paste(prepped_datas,"trans1_Over_coexp_network_all_healthy-mild-moderate-severe.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)
#pcc_sig=pcc[ pcc[,4]<=0.05    ,] #extracting significant interaction among case
#pcc_healthy_sig=pcc_healthy[ pcc_healthy[,4]<=0.05    ,] #extract significant interaction among healthy
'''
num_pcc_nodes<-length(unique(sort(cbind(pcc[,1], pcc[,2]))))
num_pcc_edges<-nrow(pcc)
num_pcc_sig_nodes<-length(unique(sort(cbind(pcc_sig[,1], pcc_sig[,2]))))
num_pcc_sig_edges<-nrow(pcc_sig)
num_pcc_healthy_nodes<-length(unique(sort(cbind(pcc_healthy[,1], pcc_healthy[,2]))))
num_pcc_healthy_edges<-nrow(pcc_healthy)
num_pcc_healthy_sig_nodes<-length(unique(sort(cbind(pcc_healthy_sig[,1], pcc_healthy_sig[,2]))))
num_pcc_healthy_sig_edges<-nrow(pcc_healthy_sig)
descp<-c("num_pcc_nodes","num_pcc_edges","num_pcc_sig_nodes","num_pcc_sig_edges","num_pcc_healthy_nodes","num_pcc_healthy_edges","num_pcc_healthy_sig_nodes", "num_pcc_healthy_sig_edges")
node_stats<-c(num_pcc_nodes,num_pcc_edges, num_pcc_sig_nodes,num_pcc_sig_edges,num_pcc_healthy_nodes,num_pcc_healthy_edges,num_pcc_healthy_sig_nodes, num_pcc_healthy_sig_edges)
print(data.frame(x=descp,y=node_stats))

write.table(pcc,paste(prepped_datas,"trans1_Over_coexp_network_all.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)
write.table(pcc_sig,paste(prepped_datas,"trans1_Over_coexp_network_significant.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)
write.table(pcc_healthy,paste(prepped_datas,"trans1_Over_coexp_network_all_healthy.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)
write.table(pcc_healthy_sig,paste(prepped_datas,"trans1_Over_coexp_network_significant_healthy.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)
'''


#for trans2
pcc1 = c()
pcc1_healthy=c(); pcc1_mild =c(); pcc1_moderate =c(); pcc1_severe =c();
for (i in 1:nrow(covid_interactome_in_trans2)){
  #print(i)
  temp_pcc = c() # this is for case (combined disease state except healthy)
  temp_pcc_healthy = c() #this is for healthy
  temp_pcc_mild = c(); temp_pcc_moderate = c(); temp_pcc_severe = c()
  
  temp_gene1 = covid_interactome_in_trans2[i,1]
  temp_gene1_row<-which(trans2$Gene_ID == temp_gene1)
  print(temp_gene1_row)
  
  temp_gene2 = covid_interactome_in_trans2[i,2]
  temp_gene2_row<-which(trans2$Gene_ID == temp_gene2)
  weight_temp_gene1<-covid_interactome_in_trans2[i,3]
  
  temp_exp1 = trans2_case[ temp_gene1_row ,2:ncol(trans2_case)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp2 = trans2_case[ temp_gene2_row ,2:ncol(trans2_case)]
  temp_pcc_test = cor.test( as.numeric(temp_exp1), as.numeric(temp_exp2),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_test
  
  temp_exp3 = trans2_healthy[ temp_gene1_row ,2:ncol(trans2_healthy)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp4 = trans2_healthy[ temp_gene2_row ,2:ncol(trans2_healthy)]
  temp_pcc_healthy_test = cor.test( as.numeric(temp_exp3), as.numeric(temp_exp4),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_healthy_test
  
  temp_exp5 = trans2_mild[ temp_gene1_row ,2:ncol(trans2_mild)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp6 = trans2_mild[ temp_gene2_row ,2:ncol(trans2_mild)]
  temp_pcc_mild_test = cor.test( as.numeric(temp_exp5), as.numeric(temp_exp6),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_mild_test
  
  temp_exp7 = trans2_moderate[ temp_gene1_row ,2:ncol(trans2_moderate)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp8 = trans2_moderate[ temp_gene2_row ,2:ncol(trans2_moderate)]
  temp_pcc_moderate_test = cor.test( as.numeric(temp_exp7), as.numeric(temp_exp8),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_moderate_test
  
  temp_exp9 = trans2_severe[ temp_gene1_row ,2:ncol(trans2_severe)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp10 = trans2_severe[ temp_gene2_row ,2:ncol(trans2_severe)]
  temp_pcc_severe_test = cor.test( as.numeric(temp_exp9), as.numeric(temp_exp10),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_severe_test
  
  temp_pcc[1] = temp_gene1
  temp_pcc[2] = temp_gene2
  temp_pcc[3] = temp_pcc_test$estimate
  temp_pcc[4] = temp_pcc_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc[3])){temp_pcc[3]=FALSE} 
  if( temp_pcc[3]>0 ){
    temp_pcc[5] = 1 #PCC >0
  }else{
    temp_pcc[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc1 = as.data.frame(rbind(pcc1,temp_pcc ))
  
  temp_pcc_healthy[1] = temp_gene1
  temp_pcc_healthy[2] = temp_gene2
  temp_pcc_healthy[3] = temp_pcc_healthy_test$estimate
  temp_pcc_healthy[4] = temp_pcc_healthy_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_healthy[3])){temp_pcc_healthy[3]=FALSE} 
  if( temp_pcc_healthy[3]>0 ){
    temp_pcc_healthy[5] = 1 #PCC >0
  }else{
    temp_pcc_healthy[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc1_healthy = as.data.frame(rbind(pcc1_healthy,temp_pcc_healthy ))
  
  temp_pcc_mild[1] = temp_gene1
  temp_pcc_mild[2] = temp_gene2
  temp_pcc_mild[3] = temp_pcc_mild_test$estimate
  temp_pcc_mild[4] = temp_pcc_mild_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_mild[3])){temp_pcc_mild[3]=FALSE} 
  if( temp_pcc_mild[3]>0 ){
    temp_pcc_mild[5] = 1 #PCC >0
  }else{
    temp_pcc_mild[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc1_mild = as.data.frame(rbind(pcc1_mild,temp_pcc_mild))
  
  temp_pcc_moderate[1] = temp_gene1
  temp_pcc_moderate[2] = temp_gene2
  temp_pcc_moderate[3] = temp_pcc_moderate_test$estimate
  temp_pcc_moderate[4] = temp_pcc_moderate_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_moderate[3])){temp_pcc_moderate[3]=FALSE} 
  if( temp_pcc_moderate[3]>0 ){
    temp_pcc_moderate[5] = 1 #PCC >0
  }else{
    temp_pcc_moderate[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc1_moderate = as.data.frame(rbind(pcc1_moderate,temp_pcc_moderate))
  
  temp_pcc_severe[1] = temp_gene1
  temp_pcc_severe[2] = temp_gene2
  temp_pcc_severe[3] = temp_pcc_severe_test$estimate
  temp_pcc_severe[4] = temp_pcc_severe_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_severe[3])){temp_pcc_severe[3]=FALSE} #fix and correct error
  if( temp_pcc_severe[3]>0 ){
    temp_pcc_severe[5] = 1 #PCC >0
  }else{
    temp_pcc_severe[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc1_severe = as.data.frame(rbind(pcc1_severe,temp_pcc_severe))
}


rm( list=ls(pat="temp_") )
rownames(pcc1)=NULL
colnames(pcc1)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

rownames(pcc1_healthy)=NULL
colnames(pcc1_healthy)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

rownames(pcc1_mild)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc1_mild)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

rownames(pcc1_moderate)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc1_moderate)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

rownames(pcc1_severe)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc1_severe)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

pcc_trans2_disease_state<-as.data.frame(mutate(pcc1_healthy, pcc1_mild$PCC,pcc1_mild$p.value,pcc1_mild$pcc_direction,pcc1_moderate$PCC,pcc1_moderate$p.value,pcc1_moderate$pcc_direction,
                                 pcc1_severe$PCC, pcc1_severe$p.value, pcc1_severe$pcc_direction, pcc1$PCC,pcc1$p.value,pcc1$pcc_direction)) #pcc across the various disease state 

colnames(pcc_trans2_disease_state)=c("GeneID1","GeneID2","PCC_healthy","p.value_healthy","pcc_direction_healthy",
                                     "PCC_mild","p.value_mild","pcc_direction_mild",
                                     "PCC_moderate","p.value_moderate","pcc_direction_moderate",
                                     "PCC_severe","p.value_severe","pcc_direction_severe",
                                     "PCC_case","p.value_case","pcc_direction_case")

save(pcc1_healthy, pcc1_mild, pcc1_moderate,pcc1_severe,pcc1,pcc_trans2_disease_state, file = "Trans2_coexpression.RData")
write.table(pcc_trans2_disease_state,paste(prepped_datas,"trans2_Su_coexp_network_all_healthy-mild-moderate-severe.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)
'''
#print number of nodes in each correlation network
num_pcc1_nodes<-length(unique(sort(cbind(pcc1[,1], pcc1[,2]))))
num_pcc1_edges<-nrow(pcc1)
num_pcc1_sig_nodes<-length(unique(sort(cbind(pcc1_sig[,1], pcc1_sig[,2]))))
num_pcc1_sig_edges<-nrow(pcc1_sig)
num_pcc1_healthy_nodes<-length(unique(sort(cbind(pcc1_healthy[,1], pcc1_healthy[,2]))))
num_pcc1_healthy_edges<-nrow(pcc1_healthy)
num_pcc1_healthy_sig_nodes<-length(unique(sort(cbind(pcc1_healthy_sig[,1], pcc1_healthy_sig[,2]))))
num_pcc1_healthy_sig_edges<-nrow(pcc1_healthy_sig)
descp1<-c("num_pcc1_nodes","num_pcc1_edges","num_pcc1_sig_nodes","num_pcc1_sig_edges","num_pcc1_healthy_nodes","num_pcc1_healthy_edges","num_pcc1_healthy_sig_nodes","num_pcc1_healthy_sig_edges")
node_stats1<-c(num_pcc1_nodes,num_pcc1_edges,num_pcc1_sig_nodes,num_pcc_sig_edges,num_pcc1_healthy_nodes,num_pcc_healthy_edges,num_pcc1_healthy_sig_nodes,num_pcc1_healthy_sig_edges )
print(data.frame(x=descp1,y=node_stats1))

write.table(pcc1,paste(prepped_datas,"trans2_Su_coexp_network_all.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)
write.table(pcc1_sig,paste(prepped_datas,"trans2_Su_coexp_network_significant.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)
write.table(pcc1_healthy,paste(prepped_datas,"trans2_Su_coexp_network_all_healthy.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)
write.table(pcc1_healthy_sig,paste(prepped_datas,"trans2_Su_coexp_network_significant_healthy.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)
'''

#for prot1
pcc2 = c()
pcc2_healthy=c(); pcc2_mild =c(); pcc2_moderate =c(); pcc2_severe =c();
for (i in 1:nrow(covid_interactome_in_prot1)){
  #print(i)
  temp_pcc = c() # this is for case (combined disease state except healthy)
  temp_pcc_healthy = c() #this is for healthy
  temp_pcc_mild = c(); temp_pcc_moderate = c(); temp_pcc_severe = c()
  
  temp_gene1 = covid_interactome_in_prot1[i,1]
  temp_gene1_row<-which(prot1$Gene_ID == temp_gene1)
  print(temp_gene1_row)
  
  temp_gene2 = covid_interactome_in_prot1[i,2]
  temp_gene2_row<-which(prot1$Gene_ID == temp_gene2)
  weight_temp_gene1<-covid_interactome_in_prot1[i,3]
  
  temp_exp1 = prot1_case[ temp_gene1_row ,2:ncol(prot1_case)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp2 = prot1_case[ temp_gene2_row ,2:ncol(prot1_case)]
  temp_pcc_test = cor.test( as.numeric(temp_exp1), as.numeric(temp_exp2),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_test
  
  temp_exp3 = prot1_healthy[ temp_gene1_row ,2:ncol(prot1_healthy)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp4 = prot1_healthy[ temp_gene2_row ,2:ncol(prot1_healthy)]
  temp_pcc_healthy_test = cor.test( as.numeric(temp_exp3), as.numeric(temp_exp4),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_healthy_test
  
  temp_exp5 = prot1_mild[ temp_gene1_row ,2:ncol(prot1_mild)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp6 = prot1_mild[ temp_gene2_row ,2:ncol(prot1_mild)]
  temp_pcc_mild_test = cor.test( as.numeric(temp_exp5), as.numeric(temp_exp6),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_mild_test
  
  temp_exp7 = prot1_moderate[ temp_gene1_row ,2:ncol(prot1_moderate)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp8 = prot1_moderate[ temp_gene2_row ,2:ncol(prot1_moderate)]
  temp_pcc_moderate_test = cor.test( as.numeric(temp_exp7), as.numeric(temp_exp8),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_moderate_test
  
  temp_exp9 = prot1_severe[ temp_gene1_row ,2:ncol(prot1_severe)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp10 = prot1_severe[ temp_gene2_row ,2:ncol(prot1_severe)]
  temp_pcc_severe_test = cor.test( as.numeric(temp_exp9), as.numeric(temp_exp10),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_severe_test
  
  temp_pcc[1] = temp_gene1
  temp_pcc[2] = temp_gene2
  temp_pcc[3] = temp_pcc_test$estimate
  temp_pcc[4] = temp_pcc_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc[3])){temp_pcc[3]=FALSE} 
  if( temp_pcc[3]>0 ){
    temp_pcc[5] = 1 #PCC >0
  }else{
    temp_pcc[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc2 = as.data.frame(rbind(pcc2,temp_pcc ))
  
  temp_pcc_healthy[1] = temp_gene1
  temp_pcc_healthy[2] = temp_gene2
  temp_pcc_healthy[3] = temp_pcc_healthy_test$estimate
  temp_pcc_healthy[4] = temp_pcc_healthy_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_healthy[3])){temp_pcc_healthy[3]=FALSE} 
  if( temp_pcc_healthy[3]>0 ){
    temp_pcc_healthy[5] = 1 #PCC >0
  }else{
    temp_pcc_healthy[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc2_healthy = as.data.frame(rbind(pcc2_healthy,temp_pcc_healthy ))
  
  temp_pcc_mild[1] = temp_gene1
  temp_pcc_mild[2] = temp_gene2
  temp_pcc_mild[3] = temp_pcc_mild_test$estimate
  temp_pcc_mild[4] = temp_pcc_mild_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_mild[3])){temp_pcc_mild[3]=FALSE} 
  if( temp_pcc_mild[3]>0 ){
    temp_pcc_mild[5] = 1 #PCC >0
  }else{
    temp_pcc_mild[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc2_mild = as.data.frame(rbind(pcc2_mild,temp_pcc_mild))
  
  temp_pcc_moderate[1] = temp_gene1
  temp_pcc_moderate[2] = temp_gene2
  temp_pcc_moderate[3] = temp_pcc_moderate_test$estimate
  temp_pcc_moderate[4] = temp_pcc_moderate_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_moderate[3])){temp_pcc_moderate[3]=FALSE} 
  if( temp_pcc_moderate[3]>0 ){
    temp_pcc_moderate[5] = 1 #PCC >0
  }else{
    temp_pcc_moderate[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc2_moderate = as.data.frame(rbind(pcc2_moderate,temp_pcc_moderate))
  
  temp_pcc_severe[1] = temp_gene1
  temp_pcc_severe[2] = temp_gene2
  temp_pcc_severe[3] = temp_pcc_severe_test$estimate
  temp_pcc_severe[4] = temp_pcc_severe_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_severe[3])){temp_pcc_severe[3]=FALSE} #fix and correct error
  if( temp_pcc_severe[3]>0 ){
    temp_pcc_severe[5] = 1 #PCC >0
  }else{
    temp_pcc_severe[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc2_severe = as.data.frame(rbind(pcc2_severe,temp_pcc_severe))
}


rm( list=ls(pat="temp_") )
rownames(pcc2)=NULL
colnames(pcc2)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

rownames(pcc2_healthy)=NULL
colnames(pcc2_healthy)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

rownames(pcc2_mild)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc2_mild)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

rownames(pcc2_moderate)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc2_moderate)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

rownames(pcc2_severe)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc2_severe)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

pcc_prot1_disease_state<-mutate(pcc2_healthy, pcc2_mild$PCC,pcc2_mild$p.value,pcc2_mild$pcc_direction,pcc2_moderate$PCC,pcc2_moderate$p.value,pcc2_moderate$pcc_direction,
                                 pcc2_severe$PCC, pcc2_severe$p.value, pcc2_severe$pcc_direction, pcc2$PCC,pcc2$p.value,pcc2$pcc_direction) #pcc across the various disease state 

colnames(pcc_prot1_disease_state)=c("GeneID1","GeneID2","PCC_healthy","p.value_healthy","pcc_direction_healthy",
                                     "PCC_mild","p.value_mild","pcc_direction_mild",
                                     "PCC_moderate","p.value_moderate","pcc_direction_moderate",
                                     "PCC_severe","p.value_severe","pcc_direction_severe",
                                     "PCC_case","p.value_case","pcc_direction_case")

save(pcc2_healthy, pcc2_mild, pcc2_moderate,pcc2_severe,pcc2, pcc_prot1_disease_state, file = "Prot1_coexpression.RData")

write.table(pcc_prot1_disease_state,paste(prepped_datas,"prot1_Overmyer_coexp_network_all_healthy-mild-moderate-severe.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)
'''
#print number of nodes in each correlation network
num_pcc2_nodes<-length(unique(sort(cbind(pcc2[,1], pcc2[,2]))))
num_pcc2_edges<-nrow(pcc2)
num_pcc2_sig_nodes<-length(unique(sort(cbind(pcc2_sig[,1], pcc2_sig[,2]))))
num_pcc2_sig_edges<-nrow(pcc2_sig)
num_pcc2_healthy_nodes<-length(unique(sort(cbind(pcc2_healthy[,1], pcc2_healthy[,2]))))
num_pcc2_healthy_edges<-nrow(pcc2_healthy)
num_pcc2_healthy_sig_nodes<-length(unique(sort(cbind(pcc2_healthy_sig[,1], pcc2_healthy_sig[,2]))))
num_pcc2_healthy_sig_edges<-nrow(pcc2_healthy_sig)
descp2<-c("num_pcc2_nodes","num_pcc2_edges","num_pcc2_sig_nodes","num_pcc2_sig_edges","num_pcc2_healthy_nodes","num_pcc2_healthy_edges","num_pcc2_healthy_sig_nodes","num_pcc2_healthy_sig_edges")
node_stats2<-c(num_pcc2_nodes,num_pcc2_edges,num_pcc2_sig_nodes,num_pcc2_sig_edges,num_pcc2_healthy_nodes,num_pcc2_healthy_edges,num_pcc2_healthy_sig_nodes,num_pcc2_healthy_sig_edges)
print(data.frame(x=descp2,y=node_stats2))

write.table(pcc2,paste(prepped_datas,"prot1_Over_coexp_network_all.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)
write.table(pcc2_sig,paste(prepped_datas,"prot1_Over_coexp_network_significant.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)
write.table(pcc2_healthy,paste(prepped_datas,"prot1_Over_coexp_network_all_healthy.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)
write.table(pcc2_healthy_sig,paste(prepped_datas,"prot1_Over_coexp_network_significant_healthy.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)
'''

#for prot2
pcc3 = c()
pcc3_healthy=c(); pcc3_mild =c(); pcc3_moderate =c(); pcc3_severe =c();
for (i in 1:nrow(covid_interactome_in_prot2)){
  #print(i)
  temp_pcc = c() # this is for case (combined disease state except healthy)
  temp_pcc_healthy = c() #this is for healthy
  temp_pcc_mild = c(); temp_pcc_moderate = c(); temp_pcc_severe = c()
  
  temp_gene1 = covid_interactome_in_prot2[i,1]
  temp_gene1_row<-which(prot2$Gene_ID == temp_gene1)
  print(temp_gene1_row)
  
  temp_gene2 = covid_interactome_in_prot2[i,2]
  temp_gene2_row<-which(prot2$Gene_ID == temp_gene2)
  weight_temp_gene1<-covid_interactome_in_prot2[i,3]
  
  temp_exp1 = prot2_case[ temp_gene1_row ,2:ncol(prot2_case)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp2 = prot2_case[ temp_gene2_row ,2:ncol(prot2_case)]
  temp_pcc_test = cor.test( as.numeric(temp_exp1), as.numeric(temp_exp2),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_test
  
  temp_exp3 = prot2_healthy[ temp_gene1_row ,2:ncol(prot2_healthy)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp4 = prot2_healthy[ temp_gene2_row ,2:ncol(prot2_healthy)]
  temp_pcc_healthy_test = cor.test( as.numeric(temp_exp3), as.numeric(temp_exp4),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_healthy_test
  
  temp_exp5 = prot2_mild[ temp_gene1_row ,2:ncol(prot2_mild)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp6 = prot2_mild[ temp_gene2_row ,2:ncol(prot2_mild)]
  temp_pcc_mild_test = cor.test( as.numeric(temp_exp5), as.numeric(temp_exp6),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_mild_test
  
  temp_exp7 = prot2_moderate[ temp_gene1_row ,2:ncol(prot2_moderate)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp8 = prot2_moderate[ temp_gene2_row ,2:ncol(prot2_moderate)]
  temp_pcc_moderate_test = cor.test( as.numeric(temp_exp7), as.numeric(temp_exp8),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_moderate_test
  
  temp_exp9 = prot2_severe[ temp_gene1_row ,2:ncol(prot2_severe)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp10 = prot2_severe[ temp_gene2_row ,2:ncol(prot2_severe)]
  temp_pcc_severe_test = cor.test( as.numeric(temp_exp9), as.numeric(temp_exp10),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_severe_test
  
  temp_pcc[1] = temp_gene1
  temp_pcc[2] = temp_gene2
  temp_pcc[3] = temp_pcc_test$estimate
  temp_pcc[4] = temp_pcc_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc[3])){temp_pcc[3]=FALSE} 
  if( temp_pcc[3]>0 ){
    temp_pcc[5] = 1 #PCC >0
  }else{
    temp_pcc[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc3 = as.data.frame(rbind(pcc3,temp_pcc ))
  
  temp_pcc_healthy[1] = temp_gene1
  temp_pcc_healthy[2] = temp_gene2
  temp_pcc_healthy[3] = temp_pcc_healthy_test$estimate
  temp_pcc_healthy[4] = temp_pcc_healthy_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_healthy[3])){temp_pcc_healthy[3]=FALSE} 
  if( temp_pcc_healthy[3]>0 ){
    temp_pcc_healthy[5] = 1 #PCC >0
  }else{
    temp_pcc_healthy[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc3_healthy = as.data.frame(rbind(pcc3_healthy,temp_pcc_healthy ))
  
  temp_pcc_mild[1] = temp_gene1
  temp_pcc_mild[2] = temp_gene2
  temp_pcc_mild[3] = temp_pcc_mild_test$estimate
  temp_pcc_mild[4] = temp_pcc_mild_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_mild[3])){temp_pcc_mild[3]=FALSE} 
  if( temp_pcc_mild[3]>0 ){
    temp_pcc_mild[5] = 1 #PCC >0
  }else{
    temp_pcc_mild[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc3_mild = as.data.frame(rbind(pcc3_mild,temp_pcc_mild))
  
  temp_pcc_moderate[1] = temp_gene1
  temp_pcc_moderate[2] = temp_gene2
  temp_pcc_moderate[3] = temp_pcc_moderate_test$estimate
  temp_pcc_moderate[4] = temp_pcc_moderate_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_moderate[3])){temp_pcc_moderate[3]=FALSE} 
  if( temp_pcc_moderate[3]>0 ){
    temp_pcc_moderate[5] = 1 #PCC >0
  }else{
    temp_pcc_moderate[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc3_moderate = as.data.frame(rbind(pcc3_moderate,temp_pcc_moderate))
  
  temp_pcc_severe[1] = temp_gene1
  temp_pcc_severe[2] = temp_gene2
  temp_pcc_severe[3] = temp_pcc_severe_test$estimate
  temp_pcc_severe[4] = temp_pcc_severe_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_severe[3])){temp_pcc_severe[3]=FALSE} #fix and correct error
  if( temp_pcc_severe[3]>0 ){
    temp_pcc_severe[5] = 1 #PCC >0
  }else{
    temp_pcc_severe[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc3_severe = as.data.frame(rbind(pcc3_severe,temp_pcc_severe))
}


rm( list=ls(pat="temp_") )
rownames(pcc3)=NULL
colnames(pcc3)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

rownames(pcc3_healthy)=NULL
colnames(pcc3_healthy)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

rownames(pcc3_mild)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc3_mild)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

rownames(pcc3_moderate)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc3_moderate)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

rownames(pcc3_severe)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc3_severe)=c("GeneID1","GeneID2","PCC","p.value","pcc_direction")

pcc_prot2_disease_state<-mutate(pcc3_healthy, pcc3_mild$PCC,pcc3_mild$p.value,pcc3_mild$pcc_direction,pcc3_moderate$PCC,pcc3_moderate$p.value,pcc3_moderate$pcc_direction,
                                pcc3_severe$PCC, pcc3_severe$p.value, pcc3_severe$pcc_direction, pcc3$PCC, pcc3$p.value, pcc3$pcc_direction) #pcc across the various disease state 

colnames(pcc_prot2_disease_state)=c("GeneID1","GeneID2","PCC_healthy","p.value_healthy","pcc_direction_healthy",
                                    "PCC_mild","p.value_mild","pcc_direction_mild",
                                    "PCC_moderate","p.value_moderate","pcc_direction_moderate",
                                    "PCC_severe","p.value_severe","pcc_direction_severe",
                                    "PCC_case","p.value_case","pcc_direction_case")

save(pcc3_healthy, pcc3_mild, pcc3_moderate,pcc3_severe,pcc3, pcc_prot2_disease_state, file = "Prot2_coexpression.RData")
write.table(pcc_prot2_disease_state,paste(prepped_datas,"prot2_Su_coexp_network_all_healthy-mild-moderate-severe.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

############################################################################################
#metabo1-metabo1 interaction
#metabo1_transpose<-t(metabo1)
#colnames(metabo1_transpose)<-metabo1_transpose[1,]
#metabo1_transpose <-as.data.frame(metabo1_transpose[2:nrow(metabo1_transpose),])
#
#metabo1_transpose[]<-lapply(metabo1_transpose, function(x) as.numeric(as.character(x)))
#metabo1_metabo1_interactome<-cor(metabo1_transpose)
#metabo1_metabo1_interactome_table<-as.data.frame(as.table(metabo1_metabo1_interactome))
#colnames(metabo1_metabo1_interactome_table)<-c("metaboliteID1","metaboliteID2","frequency")
#write.table(metabo1_metabo1_interactome_table,paste(prepped_datas,"metabo1_metabo1_interactome_all_healthy-mild-moderate-severe.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

pcc4 = c()
pcc4_healthy=c(); pcc4_mild =c(); pcc4_moderate =c(); pcc4_severe =c();
for (i in 1:nrow(metabo1_metabo1_interactome)){
  #print(i)
  temp_pcc = c() # this is for case (combined disease state except healthy)
  temp_pcc_healthy = c() #this is for healthy
  temp_pcc_mild = c(); temp_pcc_moderate = c(); temp_pcc_severe = c()
  
  #temp_gene1 = covid_interactome_in_metabo1[i,3]
  temp_gene1 = metabo1_metabo1_interactome[i,1]
  temp_gene1_row<-which(metabo1$metabolite_ID == temp_gene1)
  print(temp_gene1_row)
  
  temp_gene2 = metabo1_metabo1_interactome[i,2]
  temp_gene2_row<-which(metabo1$metabolite_ID == temp_gene2)
  #weight_temp_gene1<-covid_interactome_in_prot2[i,3]
  
  temp_exp1 = metabo1_case[ temp_gene1_row ,2:ncol(metabo1_case)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp2 = metabo1_case[ temp_gene2_row ,2:ncol(metabo1_case)]
  temp_pcc_test = cor.test( as.numeric(temp_exp1), as.numeric(temp_exp2),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_test
  
  temp_exp3 = metabo1_healthy[ temp_gene1_row ,2:ncol(metabo1_healthy)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp4 = metabo1_healthy[ temp_gene2_row ,2:ncol(metabo1_healthy)]
  temp_pcc_healthy_test = cor.test( as.numeric(temp_exp3), as.numeric(temp_exp4),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_healthy_test
  
  temp_exp5 = metabo1_mild[ temp_gene1_row ,2:ncol(metabo1_mild)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp6 = metabo1_mild[ temp_gene2_row ,2:ncol(metabo1_mild)]
  temp_pcc_mild_test = cor.test( as.numeric(temp_exp5), as.numeric(temp_exp6),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_mild_test
  
  temp_exp7 = metabo1_moderate[ temp_gene1_row ,2:ncol(metabo1_moderate)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp8 = metabo1_moderate[ temp_gene2_row ,2:ncol(metabo1_moderate)]
  temp_pcc_moderate_test = cor.test( as.numeric(temp_exp7), as.numeric(temp_exp8),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_moderate_test
  
  temp_exp9 = metabo1_severe[ temp_gene1_row ,2:ncol(metabo1_severe)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp10 = metabo1_severe[ temp_gene2_row ,2:ncol(metabo1_severe)]
  temp_pcc_severe_test = cor.test( as.numeric(temp_exp9), as.numeric(temp_exp10),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_severe_test
  
  temp_pcc[1] = temp_gene1
  temp_pcc[2] = temp_gene2
  temp_pcc[3] = temp_pcc_test$estimate
  temp_pcc[4] = temp_pcc_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc[3])){temp_pcc[3]=FALSE} 
  if( temp_pcc[3]>0 ){
    temp_pcc[5] = 1 #PCC >0
  }else{
    temp_pcc[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc4 = as.data.frame(rbind(pcc4,temp_pcc ))
  
  temp_pcc_healthy[1] = temp_gene1
  temp_pcc_healthy[2] = temp_gene2
  temp_pcc_healthy[3] = temp_pcc_healthy_test$estimate
  temp_pcc_healthy[4] = temp_pcc_healthy_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_healthy[3])){temp_pcc_healthy[3]=FALSE} 
  if( temp_pcc_healthy[3]>0 ){
    temp_pcc_healthy[5] = 1 #PCC >0
  }else{
    temp_pcc_healthy[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc4_healthy = as.data.frame(rbind(pcc4_healthy,temp_pcc_healthy ))
  
  temp_pcc_mild[1] = temp_gene1
  temp_pcc_mild[2] = temp_gene2
  temp_pcc_mild[3] = temp_pcc_mild_test$estimate
  temp_pcc_mild[4] = temp_pcc_mild_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_mild[3])){temp_pcc_mild[3]=FALSE} 
  if( temp_pcc_mild[3]>0 ){
    temp_pcc_mild[5] = 1 #PCC >0
  }else{
    temp_pcc_mild[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc4_mild = as.data.frame(rbind(pcc4_mild,temp_pcc_mild))
  
  temp_pcc_moderate[1] = temp_gene1
  temp_pcc_moderate[2] = temp_gene2
  temp_pcc_moderate[3] = temp_pcc_moderate_test$estimate
  temp_pcc_moderate[4] = temp_pcc_moderate_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_moderate[3])){temp_pcc_moderate[3]=FALSE} 
  if( temp_pcc_moderate[3]>0 ){
    temp_pcc_moderate[5] = 1 #PCC >0
  }else{
    temp_pcc_moderate[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc4_moderate = as.data.frame(rbind(pcc4_moderate,temp_pcc_moderate))
  
  temp_pcc_severe[1] = temp_gene1
  temp_pcc_severe[2] = temp_gene2
  temp_pcc_severe[3] = temp_pcc_severe_test$estimate
  temp_pcc_severe[4] = temp_pcc_severe_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_severe[3])){temp_pcc_severe[3]=FALSE} #fix and correct error
  if( temp_pcc_severe[3]>0 ){
    temp_pcc_severe[5] = 1 #PCC >0
  }else{
    temp_pcc_severe[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc4_severe = as.data.frame(rbind(pcc4_severe,temp_pcc_severe))
}


rm( list=ls(pat="temp_") )
rownames(pcc4)=NULL
colnames(pcc4)=c("metaboliteID1","metaboliteID2","PCC","p.value","pcc_direction")

rownames(pcc4_healthy)=NULL
colnames(pcc4_healthy)=c("metaboliteID1","metaboliteID2","PCC","p.value","pcc_direction")

rownames(pcc4_mild)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc4_mild)=c("metaboliteID1","metaboliteID2","PCC","p.value","pcc_direction")

rownames(pcc4_moderate)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc4_moderate)=c("metaboliteID1","metaboliteID2","PCC","p.value","pcc_direction")

rownames(pcc4_severe)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc4_severe)=c("metaboliteID1","metaboliteID2","PCC","p.value","pcc_direction")

pcc_metabo1_disease_state<-mutate(pcc4_healthy, pcc4_mild$PCC,pcc4_mild$p.value,pcc4_mild$pcc_direction,pcc4_moderate$PCC,pcc4_moderate$p.value,pcc4_moderate$pcc_direction,
                                pcc4_severe$PCC, pcc4_severe$p.value, pcc4_severe$pcc_direction, pcc4$PCC, pcc4$p.value, pcc4$pcc_direction) #pcc across the various disease state 

colnames(pcc_metabo1_disease_state)=c("metaboliteID1","metaboliteID2","PCC_healthy","p.value_healthy","pcc_direction_healthy",
                                    "PCC_mild","p.value_mild","pcc_direction_mild",
                                    "PCC_moderate","p.value_moderate","pcc_direction_moderate",
                                    "PCC_severe","p.value_severe","pcc_direction_severe",
                                    "PCC_case","p.value_case","pcc_direction_case")

save(pcc4_healthy, pcc4_mild, pcc4_moderate,pcc4_severe,pcc4, pcc_metabo1_disease_state, file = "metabo1_coexpression.RData")
write.table(pcc_metabo1_disease_state,paste(prepped_datas,"metabo1_Overmyer_coexp_network_all_healthy-mild-moderate-severe.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)






#metabo2-metabo2 interaction
#metabo2_transpose<-t(metabo2)
#colnames(metabo2_transpose)<-metabo2_transpose[1,]
#metabo2_transpose <-as.data.frame(metabo2_transpose[2:nrow(metabo2_transpose),])

#metabo2_transpose[]<-lapply(metabo2_transpose, function(x) as.numeric(as.character(x)))
#metabo2_metabo2_interactome<-cor(metabo2_transpose)
#metabo2_metabo2_interactome_table<-as.data.frame(as.table(metabo2_metabo2_interactome))
#colnames(metabo2_metabo2_interactome_table)<-c("metaboliteID1","metaboliteID2","frequency")
#write.table(metabo2_metabo2_interactome_table,paste(prepped_datas,"metabo2_metabo2_interactome_all_healthy-mild-moderate-severe.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)

pcc5 = c()
pcc5_healthy=c(); pcc5_mild =c(); pcc5_moderate =c(); pcc5_severe =c();
for (i in 1:nrow(metabo2_metabo2_interactome)){
  #print(i)
  temp_pcc = c() # this is for case (combined disease state except healthy)
  temp_pcc_healthy = c() #this is for healthy
  temp_pcc_mild = c(); temp_pcc_moderate = c(); temp_pcc_severe = c()
  
  #temp_gene1 = covid_interactome_in_metabo1[i,3]
  temp_gene1 = metabo2_metabo2_interactome[i,1]
  temp_gene1_row<-which(metabo2$metabolite_ID == temp_gene1)
  print(temp_gene1_row)
  
  temp_gene2 = metabo2_metabo2_interactome[i,2]
  temp_gene2_row<-which(metabo2$metabolite_ID == temp_gene2)
  #weight_temp_gene1<-covid_interactome_in_prot2[i,3]
  
  temp_exp1 = metabo2_case[ temp_gene1_row ,2:ncol(metabo2_case)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp2 = metabo2_case[ temp_gene2_row ,2:ncol(metabo2_case)]
  temp_pcc_test = cor.test( as.numeric(temp_exp1), as.numeric(temp_exp2),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_test
  
  temp_exp3 = metabo2_healthy[ temp_gene1_row ,2:ncol(metabo2_healthy)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp4 = metabo2_healthy[ temp_gene2_row ,2:ncol(metabo2_healthy)]
  temp_pcc_healthy_test = cor.test( as.numeric(temp_exp3), as.numeric(temp_exp4),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_healthy_test
  
  temp_exp5 = metabo2_mild[ temp_gene1_row ,2:ncol(metabo2_mild)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp6 = metabo2_mild[ temp_gene2_row ,2:ncol(metabo2_mild)]
  temp_pcc_mild_test = cor.test( as.numeric(temp_exp5), as.numeric(temp_exp6),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_mild_test
  
  temp_exp7 = metabo2_moderate[ temp_gene1_row ,2:ncol(metabo2_moderate)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp8 = metabo2_moderate[ temp_gene2_row ,2:ncol(metabo2_moderate)]
  temp_pcc_moderate_test = cor.test( as.numeric(temp_exp7), as.numeric(temp_exp8),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_moderate_test
  
  temp_exp9 = metabo2_severe[ temp_gene1_row ,2:ncol(metabo2_severe)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp10 = metabo2_severe[ temp_gene2_row ,2:ncol(metabo2_severe)]
  temp_pcc_severe_test = cor.test( as.numeric(temp_exp9), as.numeric(temp_exp10),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_severe_test
  
  temp_pcc[1] = temp_gene1
  temp_pcc[2] = temp_gene2
  temp_pcc[3] = temp_pcc_test$estimate
  temp_pcc[4] = temp_pcc_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc[3])){temp_pcc[3]=FALSE} 
  if( temp_pcc[3]>0 ){
    temp_pcc[5] = 1 #PCC >0
  }else{
    temp_pcc[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc5 = as.data.frame(rbind(pcc5,temp_pcc ))
  
  temp_pcc_healthy[1] = temp_gene1
  temp_pcc_healthy[2] = temp_gene2
  temp_pcc_healthy[3] = temp_pcc_healthy_test$estimate
  temp_pcc_healthy[4] = temp_pcc_healthy_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_healthy[3])){temp_pcc_healthy[3]=FALSE} 
  if( temp_pcc_healthy[3]>0 ){
    temp_pcc_healthy[5] = 1 #PCC >0
  }else{
    temp_pcc_healthy[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc5_healthy = as.data.frame(rbind(pcc5_healthy,temp_pcc_healthy ))
  
  temp_pcc_mild[1] = temp_gene1
  temp_pcc_mild[2] = temp_gene2
  temp_pcc_mild[3] = temp_pcc_mild_test$estimate
  temp_pcc_mild[4] = temp_pcc_mild_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_mild[3])){temp_pcc_mild[3]=FALSE} 
  if( temp_pcc_mild[3]>0 ){
    temp_pcc_mild[5] = 1 #PCC >0
  }else{
    temp_pcc_mild[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc5_mild = as.data.frame(rbind(pcc5_mild,temp_pcc_mild))
  
  temp_pcc_moderate[1] = temp_gene1
  temp_pcc_moderate[2] = temp_gene2
  temp_pcc_moderate[3] = temp_pcc_moderate_test$estimate
  temp_pcc_moderate[4] = temp_pcc_moderate_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_moderate[3])){temp_pcc_moderate[3]=FALSE} 
  if( temp_pcc_moderate[3]>0 ){
    temp_pcc_moderate[5] = 1 #PCC >0
  }else{
    temp_pcc_moderate[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc5_moderate = as.data.frame(rbind(pcc5_moderate,temp_pcc_moderate))
  
  temp_pcc_severe[1] = temp_gene1
  temp_pcc_severe[2] = temp_gene2
  temp_pcc_severe[3] = temp_pcc_severe_test$estimate
  temp_pcc_severe[4] = temp_pcc_severe_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_severe[3])){temp_pcc_severe[3]=FALSE} #fix and correct error
  if( temp_pcc_severe[3]>0 ){
    temp_pcc_severe[5] = 1 #PCC >0
  }else{
    temp_pcc_severe[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc5_severe = as.data.frame(rbind(pcc5_severe,temp_pcc_severe))
}


rm( list=ls(pat="temp_") )
rownames(pcc5)=NULL
colnames(pcc5)=c("metaboliteID1","metaboliteID2","PCC","p.value","pcc_direction")

rownames(pcc5_healthy)=NULL
colnames(pcc5_healthy)=c("metaboliteID1","metaboliteID2","PCC","p.value","pcc_direction")

rownames(pcc5_mild)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc5_mild)=c("metaboliteID1","metaboliteID2","PCC","p.value","pcc_direction")

rownames(pcc5_moderate)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc5_moderate)=c("metaboliteID1","metaboliteID2","PCC","p.value","pcc_direction")

rownames(pcc5_severe)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc5_severe)=c("metaboliteID1","metaboliteID2","PCC","p.value","pcc_direction")

pcc_metabo2_disease_state<-mutate(pcc5_healthy, pcc5_mild$PCC,pcc5_mild$p.value,pcc5_mild$pcc_direction,pcc5_moderate$PCC,pcc5_moderate$p.value,pcc5_moderate$pcc_direction,
                                  pcc5_severe$PCC, pcc5_severe$p.value, pcc5_severe$pcc_direction, pcc5$PCC, pcc5$p.value, pcc5$pcc_direction) #pcc across the various disease state 

colnames(pcc_metabo2_disease_state)=c("metaboliteID1","metaboliteID2","PCC_healthy","p.value_healthy","pcc_direction_healthy",
                                      "PCC_mild","p.value_mild","pcc_direction_mild",
                                      "PCC_moderate","p.value_moderate","pcc_direction_moderate",
                                      "PCC_severe","p.value_severe","pcc_direction_severe",
                                      "PCC_case","p.value_case","pcc_direction_case")

save(pcc5_healthy, pcc5_mild, pcc5_moderate,pcc5_severe,pcc5, pcc_metabo2_disease_state, file = "metabo2_coexpression.RData")
write.table(pcc_metabo1_disease_state,paste(prepped_datas,"metabo2_Su_coexp_network_all_healthy-mild-moderate-severe.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)





'''
LIPIDS
'''
#lipido1-lipido1 interaction
#lipido1_transpose<-t(lipido1)
#colnames(lipido1_transpose)<-lipido1_transpose[1,]
#lipido1_transpose <-as.data.frame(lipido1_transpose[2:nrow(lipido1_transpose),])

#lipido1_transpose[]<-lapply(lipido1_transpose, function(x) as.numeric(as.character(x)))
#lipido1_lipido1_interactome<-cor(lipido1_transpose)
#lipido1_lipido1_interactome_table<-as.data.frame(as.table(lipido1_lipido1_interactome))
#colnames(lipido1_lipido1_interactome_table)<-c("lipidID1","lipidID2","frequency")
#write.table(lipido1_lipido1_interactome_table,paste(prepped_datas,"lipido1_lipido1_interactome_all_healthy-mild-moderate-severe.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)
#
#splitt<-nrow(lipido1_lipido1_interactome)
#lipido1_lipido1_interactome<-lipido1_lipido1_interactome[1:splitt, ]
#
#lipido1_lipido1_interactome<-lipido1_lipido1_interactome[1:nrow(lipido1_lipido1_interactome)/2, ]
pcc6 = c()
pcc6_healthy=c(); pcc6_mild =c(); pcc6_moderate =c(); pcc6_severe =c();
for (i in 1:nrow(lipido1_lipido1_interactome)){
  #print(i)
  temp_pcc = c() # this is for case (combined disease state except healthy)
  temp_pcc_healthy = c() #this is for healthy
  temp_pcc_mild = c(); temp_pcc_moderate = c(); temp_pcc_severe = c()
  
  #temp_gene1 = covid_interactome_in_metabo1[i,3]
  temp_gene1 = lipido1_lipido1_interactome[i,1]
  temp_gene1_row<-which(lipido1$lipid_ID == temp_gene1)
  print(temp_gene1_row)
  
  temp_gene2 = lipido1_lipido1_interactome[i,2]
  temp_gene2_row<-which(lipido1$lipid_ID == temp_gene2)
  #weight_temp_gene1<-covid_interactome_in_prot2[i,3]
  
  temp_exp1 = lipido1_case[ temp_gene1_row ,2:ncol(lipido1_case)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp2 = lipido1_case[ temp_gene2_row ,2:ncol(lipido1_case)]
  temp_pcc_test = cor.test( as.numeric(temp_exp1), as.numeric(temp_exp2),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_test
  
  temp_exp3 = lipido1_healthy[ temp_gene1_row ,2:ncol(lipido1_healthy)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp4 = lipido1_healthy[ temp_gene2_row ,2:ncol(lipido1_healthy)]
  temp_pcc_healthy_test = cor.test( as.numeric(temp_exp3), as.numeric(temp_exp4),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_healthy_test
  
  temp_exp5 = lipido1_mild[ temp_gene1_row ,2:ncol(lipido1_mild)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp6 = lipido1_mild[ temp_gene2_row ,2:ncol(lipido1_mild)]
  temp_pcc_mild_test = cor.test( as.numeric(temp_exp5), as.numeric(temp_exp6),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_mild_test
  
  temp_exp7 = lipido1_moderate[ temp_gene1_row ,2:ncol(lipido1_moderate)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp8 = lipido1_moderate[ temp_gene2_row ,2:ncol(lipido1_moderate)]
  temp_pcc_moderate_test = cor.test( as.numeric(temp_exp7), as.numeric(temp_exp8),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_moderate_test
  
  temp_exp9 = lipido1_severe[ temp_gene1_row ,2:ncol(lipido1_severe)] #extracting expression of genes present in the interactome   #this section can be grouped based on cases and healthy
  temp_exp10 = lipido1_severe[ temp_gene2_row ,2:ncol(lipido1_severe)]
  temp_pcc_severe_test = cor.test( as.numeric(temp_exp9), as.numeric(temp_exp10),na.action = na.omit  ) ### Calculate intra-omics pairwise correlation  #repeat for case and healthy
  temp_pcc_severe_test
  
  temp_pcc[1] = temp_gene1
  temp_pcc[2] = temp_gene2
  temp_pcc[3] = temp_pcc_test$estimate
  temp_pcc[4] = temp_pcc_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc[3])){temp_pcc[3]=FALSE} 
  if( temp_pcc[3]>0 ){
    temp_pcc[5] = 1 #PCC >0
  }else{
    temp_pcc[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc6 = as.data.frame(rbind(pcc6,temp_pcc ))
  
  temp_pcc_healthy[1] = temp_gene1
  temp_pcc_healthy[2] = temp_gene2
  temp_pcc_healthy[3] = temp_pcc_healthy_test$estimate
  temp_pcc_healthy[4] = temp_pcc_healthy_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_healthy[3])){temp_pcc_healthy[3]=FALSE} 
  if( temp_pcc_healthy[3]>0 ){
    temp_pcc_healthy[5] = 1 #PCC >0
  }else{
    temp_pcc_healthy[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc6_healthy = as.data.frame(rbind(pcc6_healthy,temp_pcc_healthy ))
  
  temp_pcc_mild[1] = temp_gene1
  temp_pcc_mild[2] = temp_gene2
  temp_pcc_mild[3] = temp_pcc_mild_test$estimate
  temp_pcc_mild[4] = temp_pcc_mild_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_mild[3])){temp_pcc_mild[3]=FALSE} 
  if( temp_pcc_mild[3]>0 ){
    temp_pcc_mild[5] = 1 #PCC >0
  }else{
    temp_pcc_mild[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc6_mild = as.data.frame(rbind(pcc6_mild,temp_pcc_mild))
  
  temp_pcc_moderate[1] = temp_gene1
  temp_pcc_moderate[2] = temp_gene2
  temp_pcc_moderate[3] = temp_pcc_moderate_test$estimate
  temp_pcc_moderate[4] = temp_pcc_moderate_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_moderate[3])){temp_pcc_moderate[3]=FALSE} 
  if( temp_pcc_moderate[3]>0 ){
    temp_pcc_moderate[5] = 1 #PCC >0
  }else{
    temp_pcc_moderate[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc6_moderate = as.data.frame(rbind(pcc6_moderate,temp_pcc_moderate))
  
  temp_pcc_severe[1] = temp_gene1
  temp_pcc_severe[2] = temp_gene2
  temp_pcc_severe[3] = temp_pcc_severe_test$estimate
  temp_pcc_severe[4] = temp_pcc_severe_test$p.value 
  #temp_pcc[5] = weight_temp_gene1  #genemania weights
  if (is.na(temp_pcc_severe[3])){temp_pcc_severe[3]=FALSE} #fix and correct error
  if( temp_pcc_severe[3]>0 ){
    temp_pcc_severe[5] = 1 #PCC >0
  }else{
    temp_pcc_severe[5] = 0 #PCC <0
  }
  #pcc = rbind(pcc,abs(temp_pcc[2:]) ) 
  pcc6_severe = as.data.frame(rbind(pcc6_severe,temp_pcc_severe))
}


rm( list=ls(pat="temp_") )
rownames(pcc6)=NULL
colnames(pcc6)=c("lipidID1","lipidID2","PCC","p.value","pcc_direction")

rownames(pcc6_healthy)=NULL
colnames(pcc6_healthy)=c("lipidID1","lipidID2","PCC","p.value","pcc_direction")

rownames(pcc6_mild)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc6_mild)=c("lipidID1","lipidID2","PCC","p.value","pcc_direction")

rownames(pcc6_moderate)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc6_moderate)=c("lipidID1","lipidID2","PCC","p.value","pcc_direction")

rownames(pcc6_severe)=NULL # pcc_healthy=c(); pcc_mild =c(); pcc_moderate =c(); pcc_severe =c();
colnames(pcc6_severe)=c("lipidID1","lipidID2","PCC","p.value","pcc_direction")

pcc_lipid1_disease_state<-mutate(pcc6_healthy, pcc6_mild$PCC,pcc6_mild$p.value,pcc6_mild$pcc_direction,pcc6_moderate$PCC,pcc6_moderate$p.value,pcc6_moderate$pcc_direction,
                                  pcc6_severe$PCC, pcc6_severe$p.value, pcc6_severe$pcc_direction, pcc6$PCC, pcc6$p.value, pcc6$pcc_direction) #pcc across the various disease state 

colnames(pcc_lipid1_disease_state)=c("lipidID1","lipidID2","PCC_healthy","p.value_healthy","pcc_direction_healthy",
                                      "PCC_mild","p.value_mild","pcc_direction_mild",
                                      "PCC_moderate","p.value_moderate","pcc_direction_moderate",
                                      "PCC_severe","p.value_severe","pcc_direction_severe",
                                      "PCC_case","p.value_case","pcc_direction_case")

save(pcc6_healthy, pcc6_mild, pcc6_moderate,pcc6_severe,pcc6, pcc_lipid1_disease_state, file = "lipid1_coexpression.RData")
write.table(pcc_lipid1_disease_state,paste(prepped_datas,"lipido1_lipido1_interactome_all_healthy-mild-moderate-severe.txt",sep="/"),sep="\t",quote=F,col.names=T,row.names=F)


