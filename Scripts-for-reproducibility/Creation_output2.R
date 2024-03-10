#Creation of the Output 2:Matrix_WHO2021
#This matrix will contains all the glioma samples and 3 columns: the current TCGA histological subtype classes,
# the nomenclature of the glioma types according to the 2021 WHO CNS guidelines + our simplified labels

#Uploading the updated clinical matrices
load("~/Integrated_data_2021.RData")

#Creation of Matrix2
Matrix_WHO2021=as.data.frame(rbind(cbind(astro_clinic$Patient_ID,astro_clinic$histological_type),cbind(oligo_clinic$Patient_ID,oligo_clinic$histological_type),cbind(gbm_clinic2021$Patient_ID,gbm_clinic2021$histological_type),cbind(unclassified_samples$Patient_ID,unclassified_samples$histological_type)))
colnames(Matrix_WHO2021)=c("Patient_ID","TCGA-histological.type")

Matrix_WHO2021$classification.2021_complete.labels=rep(c("astrocytoma_IDHmut","oligodendroglioma_IDHmut_1p19qcodel","glioblastoma_IDHwt","unclassified"),times=c(length(astro_clinic$Patient_ID),length(oligo_clinic$Patient_ID),length(gbm_clinic2021$Patient_ID),length(unclassified_samples$Patient_ID)))
Matrix_WHO2021$classification.2021_simplified.labels=rep(c("astrocytoma","oligodendroglioma","glioblastoma","unclassified"),times=c(length(astro_clinic$Patient_ID),length(oligo_clinic$Patient_ID),length(gbm_clinic2021$Patient_ID),length(unclassified_samples$Patient_ID)))

Matrix_WHO2021=Matrix_WHO2021[order(Matrix_WHO2021$Patient_ID),] #sort the dataset following the patient names

#write.csv(Matrix_WHO2021,"Matrix_WHO2021.csv", row.names = FALSE)
