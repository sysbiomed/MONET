# Title: 2021 classification 
# Author: Roberta Coletti
#
#
#library section
library("TCGAbiolinks")
library("RTCGAToolbox")



#1st step: download the datasets from TCGA and create the PanGlioma dataset

## 1st OPTION: direct download trough getFirehoseData function 
## **uncomment lines 14 -- 37 to use this **

##LGG
#rccData <- getFirehoseData(dataset="LGG", runDate="20160128", forceDownload=TRUE, clinical=TRUE) 
##removing the column "days_to_last_known_alive", which is not present in GBM
#LGG_clinic <- rccData@clinical[,-9]
##including a column with the patient ID
#LGG_clinic$`Composite Element REF`=toupper(rownames(LGG_clinic))
#colnames(LGG_clinic)[1]=c("Patient_ID")
#LGG_clinic$Patient_ID <- chartr('.', '-', LGG_clinic$Patient_ID)
##sort the dataset following the patient names
#LGG_clinic=LGG_clinic[order(LGG_clinic$Patient_ID),] 
##remove duplicates
#LGG_clinic= LGG_clinic[!duplicated(LGG_clinic$Patient_ID),]
#
##GBM
#rccData <- getFirehoseData(dataset="GBM", runDate="20160128", forceDownload=TRUE, clinical=TRUE) 
#gbm_clinic <- rccData@clinical
#gbm_clinic$`Composite Element REF`=toupper(rownames(gbm_clinic))
#colnames(gbm_clinic)[1]=c("Patient_ID")
#gbm_clinic$Patient_ID <- chartr('.', '-', gbm_clinic$Patient_ID)
#gbm_clinic=gbm_clinic[order(gbm_clinic$Patient_ID),] 
##remove duplicates
#gbm_clinic= gbm_clinic[!duplicated(gbm_clinic$Patient_ID),]

## 2nd OPTION: upload our reference dataset
## ** comment the lines 41-42 to use the 1st option ** 
load("~/TCGA-datasets/TCGA-GBM.RData")
load("~/TCGA-datasets/TCGA-LGG.RData")

GBM_clinic=GBM_clinic[order(GBM_clinic$Patient_ID),] 
LGG_clinic=LGG_clinic[order(LGG_clinic$Patient_ID),] 

#creation of PanGlioma clinical dataset
clinic <- rbind(LGG_clinic,GBM_clinic) #merge the datasets
clinic=clinic[order(clinic$Patient_ID),] #sort the dataset following the patient names

#------------------------------------------------------------------------------------

# 2nd step: downloading the molecular profiles of both GBM and LGG

LGG_molecular_profiles=TCGAquery_subtype(tumor="lgg")
LGG_molecular_profiles$patient=as.character(LGG_molecular_profiles$patient)
#sort the dataset following the patient names
LGG_molecular_profiles=LGG_molecular_profiles[order(LGG_molecular_profiles$patient),] 

GBM_molecular_profiles=TCGAquery_subtype(tumor="gbm")
GBM_molecular_profiles$patient=as.character(GBM_molecular_profiles$patient)
GBM_molecular_profiles=GBM_molecular_profiles[order(GBM_molecular_profiles$patient),] #sort the dataset following the patient names

# match the molecular_profiles with the reference clinical dataset
LGG_molecular_profiles=LGG_molecular_profiles[which(LGG_molecular_profiles$patient %in% LGG_clinic$Patient_ID, dim(LGG_molecular_profiles$patient), arr.ind = T),]
GBM_molecular_profiles=GBM_molecular_profiles[which(GBM_molecular_profiles$patient %in% GBM_clinic$Patient_ID, dim(GBM_molecular_profiles$patient), arr.ind = T),]

# cheching if the patients are the same (and in the same order)
all(as.character(LGG_molecular_profiles$patient)==LGG_clinic$Patient_ID) #must be TRUE
all(as.character(GBM_molecular_profiles$patient)==GBM_clinic$Patient_ID) #must be TRUE

#put all together (PanGlioma molecular profiles)
molecular_profiles=rbind(LGG_molecular_profiles,GBM_molecular_profiles) #merge the datasets
molecular_profiles=molecular_profiles[order(molecular_profiles$patient),] #sort the dataset following the patient names
all(molecular_profiles$patient==clinic$Patient_ID) #must be TRUE


# -----------------------------------------------------------------------
#include the molecular information in our reference dataset
clinic$IDH.codel.subtype = molecular_profiles$IDH.codel.subtype
clinic$TERT.promoter.status=molecular_profiles$TERT.promoter.status
clinic$Chr.7.gain.Chr.10.loss=molecular_profiles$Chr.7.gain.Chr.10.loss

#save the samples without IDHinfo to create the case set for manual search on TCGA 
case_set=clinic$Patient_ID[is.na(clinic$IDH.codel.subtype)]

#write.csv(case_set, file="case_set_manual_search_TCGA.csv") #save it to perform your manual search

# -----------------------------------------------------------------------
# 2021's classification

#Astrocytoma and Oligodendroglioma groups

astro_clinic=clinic[clinic$IDH.codel.subtype=="IDHmut-non-codel",]
#remove NA rows
astro_clinic=astro_clinic[-which(is.na(astro_clinic$Patient_ID)),]

oligo_clinic=clinic[clinic$IDH.codel.subtype=="IDHmut-codel",]
#remove NA rows
oligo_clinic=oligo_clinic[-which(is.na(oligo_clinic$Patient_ID)),]

#glioblastoma group
clinic_IDHwt=clinic[clinic$IDH.codel.subtype=="IDHwt",]
clinic_IDHwt=clinic_IDHwt[-which(is.na(clinic_IDHwt$Patient_ID)),]

#based on genetic parameters
gbm_clinic2021=clinic_IDHwt[c(which(clinic_IDHwt$histological_type=="glioblastoma", dim(clinic_IDHwt$histological_type),arr.ind = T),which(clinic_IDHwt$TERT.promoter.status=="Mutant", dim(clinic_IDHwt$TERT.promoter.status),arr.ind = T),which(clinic_IDHwt$Chr.7.gain.Chr.10.loss=="Gain chr 7 & loss chr 10",  dim(clinic_IDHwt$Chr.7.gain.Chr.10.loss),arr.ind = T)),]
gbm_clinic2021= gbm_clinic2021[!duplicated(gbm_clinic2021$Patient_ID),]

#based on histological features
gbm_clinic2021=rbind(gbm_clinic2021,clinic_IDHwt[which(clinic_IDHwt$Patient_ID %in% GBM_clinic$Patient_ID),])
#remove duplicated
gbm_clinic2021= gbm_clinic2021[!duplicated(gbm_clinic2021$Patient_ID),]

#update the results of the manual search on TCGA (all samples are from TCGA-GBM project + IDHwt)
case_study_results<- read.csv("~/IDHstatus-TCGA-case_study.csv") #all of them are in TCGA-GBA project (GBM histology)
other_gbm_cases=case_study_results$Case.ID[case_study_results$Case.ID=="WT"] 
gbm_clinic2021=rbind(gbm_clinic2021,clinic[which(clinic$Patient_ID %in% other_gbm_cases),])
gbm_clinic2021=gbm_clinic2021[order(gbm_clinic2021$Patient_ID),]

#all the remaining samples are unclassified 
unclassified_samples=clinic[-which(clinic$Patient_ID %in% union(astro_clinic$Patient_ID,union(oligo_clinic$Patient_ID,gbm_clinic2021$Patient_ID))),]
unclassified_samples=unclassified_samples[order(unclassified_samples$Patient_ID),] #sort the dataset following the patient names

# -----------------------------------------------------------------------
#Outputs

#save(astro_clinic, oligo_clinic, gbm_clinic2021, unclassified_samples, file="Integrated_data_2021.RData")



