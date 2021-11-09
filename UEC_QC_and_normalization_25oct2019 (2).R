
# Written Date: June 17, 2019
# Author: Paul Tran
# Edit Date: October 25, 2019
# Editor: Lynn Tran
# Title: Nanostring Sample QC and Normalization

rm(list=ls(all=TRUE))
library(pheatmap)
library(tidyverse)
library(limma)
library(gplots)

#geometric mean function
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

#read in and reshape dataframe
setwd("C:\\Users\\lytran\\Box Sync\\Nanostring\\UEC")
ncounter_data<-read.csv("25oct2019_UEC_IGS_19batch_raw.csv")
colnames(ncounter_data)

ncounter_data1<-as.matrix(ncounter_data[,8:dim(ncounter_data)[2]]) #colnames: "Probe name", "Annotation", "accession","class name", "avg count", "%cv
rownames(ncounter_data1)<-ncounter_data$Probe.Name
samp_ids<-colnames(ncounter_data1)
length(samp_ids)

samp_ids_split<-strsplit(samp_ids,"_")
x3 = unlist(lapply(samp_ids_split, function(l) l[[3]]))
x1 = unlist(lapply(samp_ids_split, function(l) l[[1]]))

colnames(ncounter_data1)<-samp_ids

### Step 0: Visualize raw data to detect batch effects and outliers  #############
#make heatmap
my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 204)
mydf <- data.frame(row.names = paste(samp_ids), category = paste(x1))

  pheatmap(log2(ncounter_data1), color=rev(my_palette),cluster_cols = F, cluster_rows = F, 
         annotation_col = mydf,main = "Nanostring raw data",
         fontsize_row = 5,fontsize_col = 5)


"
Batch effects noted on visual inspection. Samples before August 2019 had higher raw counts, corresponding to protocol changes (limit of 3 cores per sample extraction, without pooling).
"

########      Step 1: Remove samples and genes with low counts  ################
endog_gene_index<-which(ncounter_data$Class.Name=="Endogenous")
myexpdat<-as.matrix(ncounter_data1[endog_gene_index,])
background_val<-max(ncounter_data1[grep("NEG_",ncounter_data$Probe.Name),])

#samples
sample_cutoff_nogenes_belowbackground<-84
no_genes_below_background<-apply(myexpdat,2,function(x){length(which(x<background_val))})
hist(no_genes_below_background)
samps_to_keep<-names(which(no_genes_below_background<sample_cutoff_nogenes_belowbackground))
length(samps_to_keep)
no_genes_below_background[which(no_genes_below_background>=sample_cutoff_nogenes_belowbackground)]

"
Loose cutoff of >= 50% probe binding above background. Eliminates higher proportion of sample after August as compared to proportion before August.
"

ncounter_data_lcrm<-ncounter_data1[,samps_to_keep]

rowstokeep<-which(rownames(mydf)%in%samps_to_keep)
mydf1<-data.frame(mydf[rowstokeep,])
rownames(mydf1)<-rownames(mydf)[rowstokeep]
colnames(mydf1)<-"category"

pheatmap(log2(ncounter_data_lcrm), color=rev(my_palette),cluster_cols = F, cluster_rows = F, 
         annotation_col = mydf1,main = "Nanostring raw data low count removed",
         fontsize_row = 5,fontsize_col = 5)
dim(ncounter_data_lcrm)


########      Step 2: background thresholding  ################
ncounter_data_lcrm_bksub<-ncounter_data_lcrm-background_val
ncounter_data_lcrm_bksub[which(ncounter_data_lcrm_bksub<=0)]<-1

pheatmap(log2(ncounter_data_lcrm_bksub), color=rev(my_palette),cluster_cols = T, cluster_rows = F, 
         annotation_col = mydf1,main = "Nanostring raw data low count removed",
         fontsize_row = 3,fontsize_col = 5)

"
Batch effect highlighted by background threshold.
"

########  Step 3: Normalize samples using Positive Control and Housekeeping Genes  ##########
##pos control norm
pos_names<-as.vector(ncounter_data$Probe.Name[ncounter_data$Class.Name=="Positive"])
myposdat<-as.matrix(ncounter_data_lcrm_bksub[pos_names,])

#visual inspection of pos control counts
data.frame(t(myposdat)) %>% 
  gather(Probe, Count, c("POS_A","POS_B","POS_C","POS_D","POS_E","POS_F")) -> 
  posdat_resh
ggplot(posdat_resh,aes(x=Probe,y=log2(Count)))+geom_jitter()+coord_cartesian(ylim=c(0,20))

# some qc work showed POS_F wasn't detect for all samples, so will rm from normalization
myposdat<-as.matrix(ncounter_data_lcrm_bksub[pos_names[-length(pos_names)],])
geomeans_persamp_prenorm<-apply(myposdat,2,gm_mean)
NF_PC <- geomeans_persamp_prenorm/gm_mean(geomeans_persamp_prenorm)
posdat_norm <- sweep(myposdat, 2, NF_PC, "/")

data.frame(t(posdat_norm)) %>% 
  gather(Probe, Count, c("POS_A","POS_B","POS_C","POS_D","POS_E")) -> 
  posdat_norm_resh
ggplot(posdat_norm_resh,aes(x=Probe,y=log2(Count)))+geom_jitter()+coord_cartesian(ylim=c(0,20))

#Normalize data
ncounter_lcrm_posnorm <- sweep(ncounter_data_lcrm_bksub, 2, NF_PC, "/")

## HK norm
HK_gene_name<-ncounter_data$Probe.Name[which(ncounter_data$Class.Name=="Housekeeping")]
HK_genes<-rownames(ncounter_lcrm_posnorm)[rownames(ncounter_lcrm_posnorm)%in%as.vector(HK_gene_name)]
myHKdat<-as.matrix(ncounter_lcrm_posnorm[as.vector(HK_gene_name),])

data.frame(t(myHKdat)) %>% gather(Probe, Count, as.vector(HK_gene_name)) -> HKdat_resh
ggplot(HKdat_resh,aes(x=Probe,y=log2(Count)))+
  geom_jitter()+
  coord_cartesian(ylim=c(0,20))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



#Remove bad HK genes
badHK<-c("MRPL30")
badHK_index<-which(rownames(myHKdat)%in%badHK)
myHKdat1<-myHKdat[-badHK_index,]

#calculate normalization vectors
myHKdat1_geomeans_persamp_prenorm<-apply(myHKdat1,2,gm_mean) #myHKdat1 or myHKdat
NF_HK <- myHKdat1_geomeans_persamp_prenorm/gm_mean(myHKdat1_geomeans_persamp_prenorm)
myHKdat1_norm <- sweep(myHKdat1, 2, NF_HK, "/") #myHKdat1 or myHKdat


data.frame(t(myHKdat1_norm)) %>% 
  gather(Probe, Count, as.vector(HK_gene_name)[-badHK_index]) -> 
  myHKdat1_norm_resh
ggplot(myHKdat1_norm_resh,aes(x=Probe,y=log2(Count)))+
  geom_jitter()+
  coord_cartesian(ylim=c(0,20))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5))

#HK norm data
ncounter_data_pos_HK_norm1 <- sweep(ncounter_lcrm_posnorm, 2, NF_HK, "/")

#check distribution of normalization factors
summary(NF_PC)
summary(NF_HK)

#Data after low count rm genes and samples, background thresholding, pos control norm, and HK gene norm
pheatmap(log2(ncounter_data_pos_HK_norm1), color=rev(my_palette),
         cluster_cols = T, cluster_rows = F, 
         annotation_col = mydf1,main = "Nanostring norm data low count removed",
         fontsize_row = 3,fontsize_col = 3)

"
Most of the steps listed above are standard and recommended by nanostring
"

#################  Step 4: Visualize ####################
pos_ind<-grep("POS_",rownames(ncounter_data_pos_HK_norm1))
neg_ind<-grep("NEG_",rownames(ncounter_data_pos_HK_norm1))
HK_ind<-which(rownames(ncounter_data_pos_HK_norm1)%in%HK_genes)
endo_genes<-rownames(ncounter_data_pos_HK_norm1)[-c(pos_ind,neg_ind,HK_ind)]

mydata<-log2(ncounter_data_pos_HK_norm1[endo_genes,]+1)
mydata_withHK<-log2(ncounter_data_pos_HK_norm1[c(endo_genes,HK_genes),]+1)

#PCA
#before
pca1 = prcomp(t(mydata), scale. = TRUE)
scores1 = as.data.frame(pca1$x)
ggplot(data = scores1, aes(x = PC1, y = PC2, label = rownames(scores1))) +
  geom_point(aes(color = mydf1$category)) +
  ggtitle("PCA plot of Nanostring Processed Before Batch effect removal")


#write dataframes for downstream analyses
date<-gsub("-","",Sys.Date())
# write.csv(mydata,paste0("UEC_IGS_Nanostring_Data_PostQC_",date,".csv"))
# write.csv(mydata_withHK,paste0("UEC_IGS_Nanostring_Data_PostQC_",date,"_withHK.csv"))


#################  Step 5: Batch Effect Removal ####################

library(sva)
batch = unlist(lapply(strsplit(colnames(mydata_withHK),"_"), function(l) l[[1]]))

mydata_withHK_batch <- ComBat(dat = mydata_withHK, batch = batch)
rownames(mydata_withHK_batch)

mydata_batch <- mydata_withHK_batch[rownames(mydata) %in% HK_genes == F, ]

#PCA
#after
pca2 = prcomp(t(mydata_batch), scale. = TRUE)
scores2 = as.data.frame(pca2$x)
ggplot(data = scores2, aes(x = PC1, y = PC2, label = rownames(scores2))) +
  geom_point(aes(color = as.factor(batch))) +
  ggtitle("PCA plot of Nanostring Processed After Batch effect removal")

#Data after low count rm genes and samples, background thresholding, pos control norm, and HK gene norm, along with batch effect removal
pheatmap(mydata_batch, color=rev(my_palette),
         cluster_cols = T, clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
         cluster_rows = F, 
         annotation_col = mydf1,main = "Nanostring norm data batch effect removed",
         fontsize_row = 3,fontsize_col = 3)

#scaled by gene
pheatmap(mydata_batch, color=rev(my_palette),
         cluster_cols = T, clustering_distance_cols = "euclidean", clustering_method = "ward.D2",
         cluster_rows = F, scale = "row",
         annotation_col = mydf1,main = "Nanostring norm data batch effect removed",
         fontsize_row = 3,fontsize_col = 3)

#write dataframes for downstream analyses
date<-gsub("-","",Sys.Date())
write.csv(mydata_batch,paste0("UEC_IGS_Nanostring_Data_PostQC_postbatchrem_",date,".csv"))
write.csv(mydata_withHK_batch,paste0("UEC_IGS_Nanostring_Data_PostQC_postbatchrem_",date,"_withHK.csv"))

preprocess_no_genes<-dim(ncounter_data1)[1]
preprocess_no_samples<-dim(ncounter_data1)[2]
postprocess_no_genes<-dim(mydata_batch)[1]
postprocess_no_samples<-dim(mydata_batch)[2]

preprocess_no_genes
preprocess_no_samples
postprocess_no_genes
postprocess_no_samples

sessionInfo()

