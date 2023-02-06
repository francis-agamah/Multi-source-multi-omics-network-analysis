wkdir<-setwd("D:/PhD/Corona_Project")
prepped_datas<-file.path(wkdir,"datasets/prepared_data/")
unprepped_datas1<-file.path(prepped_datas,"unprepped_data/Overmyer/")
unprepped_datas2<-file.path(prepped_datas,"unprepped_data/Su/")

#load the data
lipid_omics_count<-read.table(file.path(unprepped_datas1,"Overmyer_lipidomics_measurements_transposed.txt"),sep='\t', header = T, row.names = 1, stringsAsFactors = F,na.strings = T)
trans_omics_count<-read.table(file.path(unprepped_datas2,"trans_omics_count_2_norm_transposed_Su.txt"), comment.char='',skip=0, header=TRUE, fill=TRUE,check.names=FALSE,row.names=1)

prot_omics_count<-read.table(file.path(unprepped_datas2,"Su_proteomics_measurements_transposed_mapped-not-normalized.txt"), comment.char='',skip=0, header=TRUE, fill=TRUE,check.names=FALSE,row.names=1)

prot_omics_count1<-read.table(file.path(unprepped_datas1,"Overmyer_proteomics_transposed_mapped.txt"), comment.char='',skip=0, header=TRUE, fill=TRUE,check.names=FALSE,row.names=1)






##remove features with more than 10% missing values
lipid_omics_count_rm_miss<-lipid_omics_count[rowSums(is.na(lipid_omics_count)) <=
                                               round(dim(lipid_omics_count)[2]*0.1),]

trans_omics_count_rm_miss<-trans_omics_count[rowSums(is.na(trans_omics_count)) <=
                                               round(dim(trans_omics_count)[2]*0.1),]


prot_omics_count_rm_miss<-prot_omics_count[rowSums(is.na(prot_omics_count)) <=
                                             round(dim(prot_omics_count)[2]*0.1),]

prot_omics_count1_rm_miss<-prot_omics_count1[rowSums(is.na(prot_omics_count1)) <=
                                             round(dim(prot_omics_count1)[2]*0.1),]

#normalize function
normalize.matrix<-function(data.matrix){
  num=data.matrix-rowMeans(data.matrix, na.rm=TRUE)
  should.keep=(apply(num, 1, function(x) sd(x, na.rm=TRUE))!=0)
  return((num/apply(num,1,function(x) sd(x,na.rm=TRUE)))[should.keep,])
}

#converting all columns to numeric
trans_omics_count_rm_miss[] <- lapply(trans_omics_count_rm_miss, function(x) as.numeric(as.character(x))) #[,2:ncol(trans_omics_count_rm_miss)]
sapply(trans_omics_count_rm_miss, class)

prot_omics_count_rm_miss[] <- lapply(prot_omics_count_rm_miss, function(x) as.numeric(as.character(x))) #[,2:ncol(trans_omics_count_rm_miss)]
sapply(prot_omics_count_rm_miss, class)

prot_omics_count1_rm_miss[] <- lapply(prot_omics_count1_rm_miss, function(x) as.numeric(as.character(x))) #[,2:ncol(trans_omics_count_rm_miss)]
sapply(prot_omics_count1_rm_miss, class)

lipid_omics_count_rm_miss_norm<-normalize.matrix(lipid_omics_count_rm_miss)
trans_omics_count_rm_miss_norm<-normalize.matrix(trans_omics_count_rm_miss)
prot_omics_count_rm_miss_norm<-normalize.matrix(prot_omics_count_rm_miss)
prot_omics_count1_rm_miss_norm<-normalize.matrix(prot_omics_count1_rm_miss)


write.table(lipid_omics_count_rm_miss_norm, paste(prepped_datas,"lipid_omics_count_rm_miss_norm_Overmyer.tsv", sep = "/"),sep="\t",quote=F,col.names=T, row.names = T)
write.table(trans_omics_count_rm_miss_norm,file=paste(prepped_datas,"trans_omics_count_2_norm_transposed_Su.tsv",sep = "/"),sep = '\t', row.names = T)
write.table(prot_omics_count_rm_miss_norm,file=paste(prepped_datas,"prot_omics_count_2_norm_transposed_mapped_Su.tsv",sep = "/"),sep = '\t', row.names = T)
write.table(prot_omics_count1_rm_miss_norm,file=paste(prepped_datas,"prot_omics_count_norm_transposed_mapped_Overmyer.tsv",sep = "/"),sep = '\t', row.names = T)
