#Creation of the Output 3: Comparison between the current TCGA histological subtype, and the new 2016 and 2021 WHO classifications of CNS

Matrix_WHO2016 <- read.csv("~/Final-outputs/Matrix_WHO2016.csv")
Matrix_WHO2021 <- read.csv("~/Final-outputs/Matrix_WHO2021.csv")

#check the order
all(Matrix_WHO2016$Patient_ID==Matrix_WHO2021$Patient_ID) #must be TRUE

#Creation of SIMPLIFIED_CLASSIFICATION_TCGA_2016_2021 matrix 

Matrix=Matrix_WHO2016[,c(1,2)]
Matrix$classification.2016=Matrix_WHO2016$classification.2016_simplified.labels
Matrix$classification.2021=Matrix_WHO2021$classification.2021_simplified.labels

#write.csv(Matrix,"SIMPLIFIED_CLASSIFICATION_TCGA_2016_2021.csv", row.names = FALSE)

