# MONET
**Multi-omic networks in gliomas**, ref. PTDC/CCI-BIO/4180/2020

Here are provided the R scripts that reproduce the methodology for updating the tumour classification of TCGA glioma samples in accordance with WHO-2016 and WHO-2021 guidelines, presented in the preprint available at https://doi.org/10.1101/2023.02.19.529134.

- In folder **Scripts-for-reproducibility** are the two main R scripts:
  - "2016-classification.R", and
  - "2021-classification.R",
  
that, respectively, implement the Method-2016 and Method-2021 described.

- In folder **Final-outputs**, are 3 R scripts' files, named:
  - "Creation_output1.R"
  - "Creation_output2.R"
  - "Creation_output3.R"

that create, respectively, the final 3 output files, named:
   - "Matrix_WHO2016.csv",
   - "Matrix_WHO2021.csv",
   - "SIMPLIFIED_CLASSIFICATION_TCGA_2016_2021.csv",

with characteristics as described in the preprint.
