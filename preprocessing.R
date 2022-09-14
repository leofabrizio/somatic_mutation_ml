#project_work
rm(list = ls())

library(tidyverse)
library(NMF)
library(Hmisc)

#project folder
setwd('C:/Users/Fabrizio/Documents/corso_ml_unipd/project_work/dataset1/')

#import file with sample and class names
samples <- read.table('samples_labels.txt', col.names = c('patient','cancer_type'))

#import mutation file
mutations <- read_tsv('snvs.tsv', col_names = FALSE)

#import important genes file
driver_mut <- read.table('Compendium_Cancer_Genes.txt')

#I identify the missing patient in the mutation file
setdiff(samples$patient,mutations$X1)

#anmd I remove it
samples_clean <- filter(samples, patient != "TCGA-BP-4765")

#ricontrollo
#setdiff(samples_clean$patient,mutations$X1)

#join the two df
mutations <- rename(mutations, patient = X1)
data <- merge(samples_clean,mutations, by = "patient")

#filter out the genes not considered as relevant for each patient
#I create empty df with columns=names of relevant genes
dummy_df <- setNames(data.frame(matrix(ncol = length(driver_mut$V1), nrow = 0)), driver_mut$V1)

nomi_col <- names(dummy_df)

#create binary mutation matrix
for(i in 1:nrow(data)) {
    vettore <- c()
    geni_paz <- as.list(data[i,3:288]) #extract in list the names of mutated genes in patient i
    geni_imp <- intersect(geni_paz,driver_mut$V1) #I check which of thesde genes are likely important
    #print(geni_imp)
    for (n in nomi_col) {
      if (n %in% geni_imp) {
        vettore <- append(vettore,1) #for each col of binary matrix assign 1 if that gene exists
      }
      else {
        vettore <- append(vettore,0) #or 0 if that gene is not present
      }
    }
    dummy_df[i,] <- vettore #append 1 and 0 vector to mutation matrix 
}

#remove cols with only zeros (42 relevant mutated genes are not present in any patient)
dummy_df %>% select_if(colSums(.) != 0) -> dummy_df

#plot histogram
dfmat <- as.matrix(dummy_df) #hist function needs a matrix
hist(dfmat)

#I check the number of cancer types in the dataset
unique(data$cancer_type) # 11

#add column patient ID and cancer type to binary mutation df
dummy_df$patient <- data$patient
dummy_df$cancer_type <- data$cancer_type
dummy_df <- dummy_df[,c(527,528,1:526)] #reorder cols
#table(dummy_df$cancer_type) #n° patients for cancer type. 

#verifico se ci sono pazienti che non hanno geni mutati dall'elenco 'Compendium Cancer Genes'
# r0 <- which(rowSums(dummy_df[,c(3:528)]) == 0)
# length(r0) #ce ne sono 96 
# #vediamo distribuzione di pazienti di cui sopra x tipo tumore
# pazienti <- dummy_df$cancer_type[c(r0)]
# table(pazienti)

#dummy_df <- dummy_df[-c(r0),] #li rimuovo

#I save clean df
write.csv(dummy_df, file="cleaned_data_def.csv", row.names = FALSE) #salvo in csv

#plot n° patients for cancer type
barplot(table(dummy_df$cancer_type), ylim = c(0,800),cex.names=0.7, main = 'Patients by Cancer Type', ylab = 'N° patients')

#compute gene mutation frequency in the whole dataset
df1 <- dummy_df[,-c(1:2)] %>%
  summarise(across(everything(), list(mean)))

ord_df1 <- as.vector(df1)
ord_df1 <- unlist(ord_df1)
ord_df1 <- sort(ord_df1, decreasing = TRUE)
ord_df1 <- ord_df1[ord_df1 > 0.05] #select mutated genes with freq higher than > .05
print(ord_df1)

plot(ord_df1,xaxt="n", ylim = c(0,0.5), main = 'Genes with most frequent mutation in the entire sample', xlab = 'Gene', ylab = 'Freq')
axis(1,at=1:length(ord_df1),labels=names(ord_df1))

#compute mutation frequency for each cancer type
df2 <- dummy_df[,-1] %>%
  group_by(cancer_type) %>%  
  summarise(across(everything(), list(mean)))

#for each group, I order in decreasing fashion the genes for mutation freq
ord_df2 <- as.vector(df2[11,2:ncol(df2)])
ord_df2 <- unlist(ord_df2)
ord_df2 <- sort(ord_df2, decreasing = TRUE)
ord_df2 <- ord_df2[1:10] #select first 10 genes 
print(ord_df2)

plot(ord_df2,xaxt="n", main = 'Genes with most frequent mutation in the UCEC Cancer Type', xlab = 'Gene', ylab = 'Freq')
axis(1,at=1:length(ord_df2),cex.axis=0.8,labels=names(ord_df2))
