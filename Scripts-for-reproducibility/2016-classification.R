# 2021 classification 
#

#library

library("TCGAbiolinks")

library(knitr) 
library(dplyr)
library(tibble)
library(SummarizedExperiment)
library("RTCGAToolbox")
library("stringr")


#1st step: download the datasets from TCGA

## 1st OPTION: direct download trough getFirehoseData function 
## **uncomment lines 21 -- 42 to use this **

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
## ** comment the lines 46-47 to use the 1st option ** 
load("~/TCGA-datasets/TCGA-GBM.RData")
load("~/TCGA-datasets/TCGA-LGG.RData")
GBM_clinic=GBM_clinic[order(GBM_clinic$Patient_ID),] 
LGG_clinic=LGG_clinic[order(LGG_clinic$Patient_ID),] 
gbm_clinic=GBM_clinic
rm(GBM_clinic)

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
GBM_molecular_profiles=GBM_molecular_profiles[which(GBM_molecular_profiles$patient %in% gbm_clinic$Patient_ID, dim(GBM_molecular_profiles$patient), arr.ind = T),]

# cheching if the patients are the same (and in the same order)
all(as.character(LGG_molecular_profiles$patient)==LGG_clinic$Patient_ID) #must be TRUE
all(as.character(GBM_molecular_profiles$patient)==gbm_clinic$Patient_ID) #must be TRUE


# -----------------------------------------------------------------------
#include the molecular information in our reference dataset
LGG_clinic$IDH.codel.subtype = LGG_molecular_profiles$IDH.codel.subtype
gbm_clinic$IDH.codel.subtype = GBM_molecular_profiles$IDH.codel.subtype

# -----------------------------------------------------------------------
# 2016's classification

#Astrocytoma, Oligodendroglioma and mixed glioma groups

astro_clinic=LGG_clinic[LGG_clinic$histological_type=="astrocytoma",]
oligo_clinic=LGG_clinic[LGG_clinic$histological_type=="oligodendroglioma",]
mix_clinic=LGG_clinic[LGG_clinic$histological_type=="oligoastrocytoma",]

#Classification of mixed glioma samples
astro_clinic=rbind(astro_clinic,mix_clinic[mix_clinic$IDH.codel.subtype=="IDHmut-non-codel",])
astro_clinic=astro_clinic[order(astro_clinic$Patient_ID),] #sort the dataset following the patient names

oligo_clinic=rbind(oligo_clinic,mix_clinic[mix_clinic$IDH.codel.subtype=="IDHmut-codel",])
oligo_clinic=oligo_clinic[order(oligo_clinic$Patient_ID),] #sort the dataset following the patient names

unclassified_samples=mix_clinic[mix_clinic$IDH.codel.subtype=="IDHwt",]
# -----------------------------------------------------------------------
#Outputs

save(astro_clinic, oligo_clinic, gbm_clinic, unclassified_samples, file="Integrated_data_2016.RData")



