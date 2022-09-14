# somatic_mutation_ml
## Classification and prediction of cancer type using somatic mutation profiles and machine learning approaches

This repository contains a dataset including information about somatic mutations in cancer patients as well as the type of cancer for each patient. 
In particular, data nare contained in three text files:

**snvs.txt** -> mutation file for each patient. Each row is assigned to a specific patient. The first element of the row is the patient ID whereas the following elements are the names of mutated genes in that patient. Each element is separated by a TAB character.

**samples_labels.txt** -> each row is a patient. The first element of the row is the patient ID whereas the second element is the cancer type.

**Compendium_Cancer_Genes.txt** -> it is the list of genes that are considered as relevant in cancer development. This list is useful in case it would not be possible to consider all the genes in the analyses.

Please see the following website for more details:
https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga

The project aimed at answering to three different questions, i.e.:
**1) Is it possible to predict cancer type based on genes with somatic mutation in a patient?**
**2) Is there a «small» set of genes having a good predictive power, or at least as good as the entire set of genes?**
**3) Does the patient grouping based on similarity of mutated genes reflect the grouping based on cancer type?**

Using R software I answered to these three questions. The file **preprocessing.R** contains the R code I wrote to do some data preprocessing whereas the file **analyses.R ** contains the R code I wrote to answer the research questions using machine learning methods.

Finally, the file **project_work_ENG.pdf** contains a powerpoint presentation highlighting the main results of my analyses.

This project was done and defended as a requirement for the 2nd level Master in "Machine learning and big data for precision medicine and biomedical research" of the Università degli Studi di Padova, Italy.



