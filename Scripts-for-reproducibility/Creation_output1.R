#Cration of the Output 1:Matrix_WHO2016
#This matrix will contains all the glioma samples and 3 columns: the current TCGA histological subtype classes,
# the nomenclature of the glioma types according to the 2016 WHO CNS guidelines + our simplified labels

#Uploading the updated clinical matrices
load("~/Integrated_data_2016.RData")

#2016 NOMENCLATURE:
#diffuse_astrocytoma_IDHmut (if astrocytoma or mixed glioma with IDH mutation)
#oligodendroglioma_IDHmut_codel  (if oligodendroglioma or mixed-glioma with IDH mutation and 1p/19q codeletion)
#diffuse_astrocytoma_IDHwt (if astrocytoma with IDHwt)
#oligodendroglioma_NOS (if oligodendroglioma with either IDHwt or IDHmut but 1p/19q non-codel, or  missing information on IDH mutation and 1p/19q codeletion status)
#diffuse_astrocytoma_NOS (if astrocytoma with missing information on IDH mutation status)
#glioblastoma_IDHwt (if glioblastoma with IDHwt)
#glioblastoma_IDHmut (if glioblastoma with IDH mutation)
#glioblastoma_NOS (if glioblastoma with missing information on IDH mutation status)
#unclassified (if mixed glioma with IDHwt)

#astrocytoma
sample_astro_2016_IDHmut=astro_clinic$Patient_ID[substr(astro_clinic$IDH.codel.subtype,1,6)=="IDHmut"]
sample_astro_2016_IDHwt=astro_clinic$Patient_ID[substr(astro_clinic$IDH.codel.subtype,1,5)=="IDHwt"]
sample_astro_2016_NOS=astro_clinic$Patient_ID[is.na(astro_clinic$IDH.codel.subtype)]
sample_astro_2016_IDHmut=sample_astro_2016_IDHmut[!is.na(sample_astro_2016_IDHmut)]
sample_astro_2016_IDHwt=sample_astro_2016_IDHwt[!is.na(sample_astro_2016_IDHwt)]


sample_oligo_2016_IDHmut_codel=oligo_clinic$Patient_ID[substr(oligo_clinic$IDH.codel.subtype,1,12)=="IDHmut-codel"]
sample_oligo_2016_NOS=oligo_clinic$Patient_ID[-which(oligo_clinic$Patient_ID %in% sample_oligo_2016_IDHmut_codel)]
sample_oligo_2016_IDHmut_codel=sample_oligo_2016_IDHmut_codel[!is.na(sample_oligo_2016_IDHmut_codel)]

sample_gbm_2016_NOS=gbm_clinic$Patient_ID[is.na(gbm_clinic$IDH.codel.subtype)]
sample_gbm_2016_IDHmut=gbm_clinic$Patient_ID[substr(gbm_clinic$IDH.codel.subtype,1,6)=="IDHmut"]
sample_gbm_2016_IDHmut=sample_gbm_2016_IDHmut[!is.na(sample_gbm_2016_IDHmut)]
sample_gbm_2016_IDHwt=gbm_clinic$Patient_ID[substr(gbm_clinic$IDH.codel.subtype,1,6)=="IDHwt"]
sample_gbm_2016_IDHwt=sample_gbm_2016_IDHwt[!is.na(sample_gbm_2016_IDHwt)]

sample_unclassified_2016=unclassified_samples$Patient_ID

# Creation of the Matrix_WHO2016

Matrix_WHO2016=as.data.frame(rbind(cbind(astro_clinic$Patient_ID,astro_clinic$histological_type),cbind(oligo_clinic$Patient_ID,oligo_clinic$histological_type),cbind(gbm_clinic$Patient_ID,gbm_clinic$histological_type),cbind(unclassified_samples$Patient_ID,unclassified_samples$histological_type)))
colnames(Matrix_WHO2016)=c("Patient_ID","TCGA-histological.type")
Matrix_WHO2016=Matrix_WHO2016[order(Matrix_WHO2016$Patient_ID),] #sort the dataset following the patient names

#2016 WHOCNS 
samples=c(sample_astro_2016_IDHmut, sample_astro_2016_IDHwt,sample_astro_2016_NOS, 
          sample_oligo_2016_IDHmut_codel,sample_oligo_2016_NOS,
          sample_gbm_2016_NOS,sample_gbm_2016_IDHmut,sample_gbm_2016_IDHwt, sample_unclassified_2016)
WHO2016=as.data.frame(samples)
WHO2016$classification.2016=rep(c("diffuse_astrocytoma_IDHmut","diffuse_astrocytoma_IDHwt","diffuse_astrocytoma_NOS","oligodendroglioma_IDHmut_1p19qcodel","oligodendroglioma_NOS","glioblastoma_NOS","glioblastoma_IDHmut","glioblastoma_IDHwt","unclassified"),
                                times=c(length(sample_astro_2016_IDHmut),length(sample_astro_2016_IDHwt),length(sample_astro_2016_NOS),length(sample_oligo_2016_IDHmut_codel),length(sample_oligo_2016_NOS),length(sample_gbm_2016_NOS),length(sample_gbm_2016_IDHmut),length(sample_gbm_2016_IDHwt),length(sample_unclassified_2016)))
WHO2016=WHO2016[order(WHO2016$samples),] #sort the dataset following the patient names

all(as.character(Matrix_WHO2016$Patient_ID)==as.character(WHO2016$samples) )#must be TRUE

Matrix_WHO2016$classification.2016_complete.labels=WHO2016$classification.2016
Matrix_WHO2016$classification.2016_simplified.label=sub("\\_.*", "", Matrix_WHO2016$classification.2016_complete.labels)
Matrix_WHO2016$classification.2016_simplified.label[Matrix_WHO2016$classification.2016_simplified.label=="diffuse"]="astrocytoma"

#write.csv(Matrix_WHO2016,"Matrix_WHO2016.csv", row.names = FALSE)
