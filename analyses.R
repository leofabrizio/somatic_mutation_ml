#analyses

rm(list = ls())

library(tidyverse)
library(NMF)
library(Hmisc)
library(randomForest)
library(caret)
library(Boruta)
library(e1071)


#project folder
setwd('C:/Users/fabri/Documents/corso_ml_unipd/project_work/dataset1/')

#import file with sample and class names
data <- read.csv('cleaned_data_def.csv')


#########################################################################

#QUESTION1: Is it possible to predict cancer type based on genes with somatic
#mutation in a patient?

df <- data
df$patient <- NULL

#create numeric classes
dict_names <- unique(df$cancer_type)
int_names <- match(df$cancer_type,dict_names)
df$class <- int_names
df <- df[,c(528,1:527)]
df$cancer_type <- NULL
df$class <- as.factor(df$class)

#define train and test splits
n_cv <- 5

set.seed(1234)
index <- sample(1:nrow(df),round(0.75*nrow(df)))
train_ <- df[index,]
test_ <- df[-index,]

#create function to evaluate test set accuracy in various models

model_res <- function(model) {
  train_acc <- max(model$results$Accuracy)
  
  predictions <- predict(model, test_ %>% select(-class))
  matrice_confusione <- confusionMatrix(predictions, test_$class)
  test_acc <- matrice_confusione$overall[1]
  
  tibble(
    method = model$method,
    train = train_acc,
    test = test_acc,
    tune_name = names(model$bestTune),
    tune_value = model$bestTune[[1]]
  )
}

# CV randomforest

set.seed(1234)
tune_grid <- expand.grid(
  mtry = seq_len(sqrt(ncol(train_))+10)
)
rf_seed <- c(
  map(seq_len(n_cv),                        
      ~ nrow(tune_grid) * (.x - 1) +  
        seq_len(nrow(tune_grid))      
  ),
  nrow(tune_grid) * n_cv + 1      
)
rf_model <- train(
  x = train_ %>% select(-class),
  y = train_[["class"]],
  method = "rf",
  metric = "Accuracy",
  trControl = trainControl(
    method = "cv",
    number = n_cv,
    seeds = rf_seed
  ),
  tuneGrid = tune_grid
)

rf_model
plot(rf_model)

#info
rf_res <- model_res(rf_model) 
rf_res
varImpPlot(rf_model$finalModel)

# #prova con strata (downsampling)
# table(train_$class)
# nmin <- sum(train_$class == "11")
# 
# set.seed(1234)
# tune_grid <- expand.grid(
#   mtry = seq_len(sqrt(ncol(train_))+10)
# )
# rf_seed <- c(
#   map(seq_len(n_cv),                        
#       ~ nrow(tune_grid) * (.x - 1) +  
#         seq_len(nrow(tune_grid))      
#   ),
#   nrow(tune_grid) * n_cv + 1      
# )
# 
# rf_model_ds <- train(
#   x = train_ %>% select(-class),
#   y = train_[["class"]],
#   method = "rf",
#   metric = "Accuracy",
#   trControl = trainControl(
#     method = "cv",
#     number = n_cv,
#     seeds = rf_seed
#   ),
#   strata = train_$class,
#   sampsize = rep(nmin, 11),
#   tuneGrid = tune_grid
# )
# 
# rf_model_ds
# plot(rf_model_ds)
# 
# #info
# rf_ds_res <- model_res(rf_model_ds) 
# rf_ds_res
# varImpPlot(rf_model_ds$finalModel)


# CV cond tree
set.seed(1234)
tune_grid <- expand.grid(
  mincriterion = seq(from = 0.05, to = 0.99, length.out = 100)
)
cidt_seed <- c(
  map(seq_len(n_cv),                        # list of seeds for cv steps
      ~ nrow(tune_grid) * (.x - 1) +  # starting point
        seq_len(nrow(tune_grid))      # consecutive sequence of integers
  ),
  nrow(tune_grid) * n_cv + 1      # single last seed for the final model
)
cidt_model <- train(
  x = train_ %>% select(-class),
  y = train_$class,
  method = "ctree",
  metric = "Accuracy",
  trControl = trainControl(
    method = "cv",
    number = n_cv,
    seeds = cidt_seed
  ),
  tuneGrid = tune_grid
)

plot(cidt_model)

cidt_model$finalModel
plot(cidt_model$finalModel)

cidt_res <- model_res(cidt_model) 
cidt_res

# CV KNN
set.seed(1234)
tune_grid <- expand.grid(
  k = seq(1,31,by=2)
)
knn_seed <- c(
  map(seq_len(n_cv),                        
      ~ nrow(tune_grid) * (.x - 1) +  
        seq_len(nrow(tune_grid))      
  ),
  nrow(tune_grid) * n_cv + 1      
)
knn_model <- train(
  x = train_ %>% select(-class),
  y = train_[["class"]],
  method = "knn",
  metric = "Accuracy",
  trControl = trainControl(
    method = "cv",
    number = n_cv,
    seeds = knn_seed
  ),
  tuneGrid = tune_grid
)

knn_model
plot(knn_model)

#info
knn_res <- model_res(knn_model) 
knn_res

# CV SVM
set.seed(1234)

tc <- tune.control(cross = 5)
svm_model = tune(svm, class~., data = train_, probability = TRUE, kernel = "linear", ranges=list(cost=c(0.001, 0.01, 0.1, 1,5,10,100)), tunecontrol = tc)

#let's look at the results and get the best model
summary(svm_model)
bestmod=svm_model$best.model
summary(bestmod)

svm_res <- predict(bestmod, test_ %>% select(-class), decision.values = TRUE, probability = TRUE)
svm_res

#table(predict=svm_res, truth=test_$class)
matrice_confusione <- confusionMatrix(svm_res, test_$class)
matrice_confusione
test_acc <- matrice_confusione$overall[1]
train_acc <- 1 - 0.410793 #minus error

svm_res <- tibble(
  method = 'svm',
  train = train_acc,
  test = test_acc,
  tune_name = 'cost',
  tune_value = svm_model$best.parameters$cost
)


#PLOT the results of the 4 models
res <- bind_rows(
  cidt_res,
  rf_res,
  knn_res,
  svm_res
)
res %>%
  mutate(method = as.factor(method)) %>% 
  pivot_longer(
    cols = c(train, test),
    names_to = "set",
    values_to = "Accuracy"
  ) %>% 
  ggplot(aes(method, Accuracy, color = set)) +
  ggtitle("Model Comparison")+
  theme(plot.title = element_text(lineheight=.8, face="bold"))+
  geom_point(size = 4)
#varImpPlot(rf_model$finalModel)


#Plot confusion matrix of the best model on test data
predictions <- predict(rf_model, test_ %>% select(-class))
matrice_confusione <- confusionMatrix(predictions, test_$class)
matrice_confusione

plt <- as.data.frame(matrice_confusione$table)
plt$Prediction <- factor(plt$Prediction, levels=rev(levels(plt$Prediction)))

rf <- prop.table(matrice_confusione$table, 2) #calcolo freq relative
#convert to list for ggplot
my_list2 <- list()
for(i in 1:ncol(rf)) { 
  my_list2[[i]] <- rf[ , i]
}
prova <- unlist(my_list2)  
plt$rel_freq <- prova

#plot
ggplot(plt, aes(Reference,Prediction, fill= rel_freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="#009194") +
  labs(x = "Reference",y = "Prediction") +
  scale_x_discrete(labels=dict_names) +
  scale_y_discrete(labels=rev(dict_names))+
  ggtitle("Confusion Matrix - Test Data - Best Model (RF)")+
  theme(plot.title=element_text(size=13,face="bold"),
        axis.title.x=element_text(size=11,face="bold"),
        axis.title.y=element_text(size=11,face="bold"))

#Plot sensitivity for each class
barplot(matrice_confusione$byClass[,1], ylim = c(0,1.0), 
        main = "Classification Sensitivity- Test Data - Best Model (RF)", xlab = "Cancer Type", ylab = "Sensitivity",
        names.arg = dict_names, cex.names = 0.5)

#plot balanced accuracy for each class
barplot(matrice_confusione$byClass[,11], ylim = c(0,1.0), 
        main = "Balanced Accuracy - Test Data - Best Model (RF)", xlab = "Cancer Type", ylab = "Balanced Accuracy",
        names.arg = dict_names, cex.names = 0.5)

#########################################################################

#QUESTION 2: Is there a small set of genes having a good predictive power, 
#or at least as good as the entire set of genes (or all genes in "Compendium_Cancer_Genes.txt")?

#RandomForest APPROACH

feat_imp <- importance(rf_model$finalModel, type = 2) #variable relevance in rf best model

#create function to create a formula including selected genes
form <- function(thresh) {
  idx <- which(feat_imp > thresh)
  sel_var <- feat_imp[idx,]
  sel_var <- names(sel_var)
  nvar <- length(sel_var)  #n° var selected
  f <- as.formula(paste("class ~", paste(sel_var[!sel_var %in% "Class"], collapse = " + ")))
}

#create function to measure model accuracy in test set
mod_res1 <- function(rf_sel_var) {
  p_ <- predict(rf_sel_var, test_)
  res_rf_sel_var <- confusionMatrix(p_, test_$class)
  res_rf_sel_var <- res_rf_sel_var$overall[1]
  res_rf_sel_var
}

#1st genes selection
thresh <- 2 #select in sequence 2,3,4,5,6,10 as Gini score threshold values

f <- form(thresh)

#new rf with selected genes
rf_sel_var <- randomForest(f, data=train_, proximity=TRUE) 

print(rf_sel_var)

#accuracy in test set
res_rf_sel_var <- mod_res1(rf_sel_var)

#2nd genes selection
thresh <- 3 #select in sequence 2,3,4,5,6,10 as Gini score threshold values

f <- form(thresh)

#new rf with selected genes
rf_sel_var1 <- randomForest(f, data=train_, proximity=TRUE) 

print(rf_sel_var1)

#test set accuracy
res_rf_sel_var1 <- mod_res1(rf_sel_var1)

#3rd genes selection
thresh <- 4 #select in sequence 2,3,4,5,6,10 as Gini score threshold values

f <- form(thresh)

#new rf with selected genes
rf_sel_var2 <- randomForest(f, data=train_, proximity=TRUE) 

print(rf_sel_var2)

#test set accuracy
res_rf_sel_var2 <- mod_res1(rf_sel_var2)


#4th genes selection
thresh <- 5 #select in sequence 2,3,4,5,6,10 as Gini score threshold values

f <- form(thresh)

#new rf with selected genes
rf_sel_var3 <- randomForest(f, data=train_, proximity=TRUE) 

print(rf_sel_var3)

#test set accuracy
res_rf_sel_var3 <- mod_res1(rf_sel_var3)

#5th genes selection
thresh <- 6 #select in sequence 2,3,4,5,6,10 as Gini score threshold values

f <- form(thresh)

#new rf with selected genes
rf_sel_var4 <- randomForest(f, data=train_, proximity=TRUE) 

print(rf_sel_var4)

#test set accuracy
res_rf_sel_var4 <- mod_res1(rf_sel_var4)

#6th genes selection
thresh <- 10 #select in sequence 2,3,4,5,6,10 as Gini score threshold values

f <- form(thresh)

#new rf with selected genes
rf_sel_var5 <- randomForest(f, data=train_, proximity=TRUE) 

print(rf_sel_var5)

#test set accuracy
res_rf_sel_var5 <- mod_res1(rf_sel_var5)

#new rf with all the genes as reference
rf_all_var <- randomForest(class~., data=train_,proximity=TRUE) 

print(rf_all_var)

#test set accuracy
res_rf_all_var <- mod_res1(rf_all_var)


#BORUTA APPROACH: 
#5-folds CV. For each fold I extract the relevant genes with Boruta and, 
#in the end, I select relevant genes in common to all the folds

folds <- createFolds(train_$class, k = 5, list = TRUE, returnTrain = FALSE)

#RUN1
boruta1 <- Boruta(class~., data=train_[folds$Fold1,], doTrace = 2, maxRuns = 500)
print(boruta1)

bor1 <- TentativeRoughFix(boruta1) #forcing the choice for remaining genes
print(bor1)

sel_var1 <- getSelectedAttributes(bor1) #selected variables
print(sel_var1)

#RUN2
boruta2 <- Boruta(class~., data=train_[folds$Fold2,], doTrace = 2, maxRuns = 500)
print(boruta2)

bor2 <- TentativeRoughFix(boruta2) #forcing the choice for remaining genes
print(bor2)

sel_var2 <- getSelectedAttributes(bor2) #selected variables
print(sel_var2)

#RUN3
boruta3 <- Boruta(class~., data=train_[folds$Fold3,], doTrace = 2, maxRuns = 500)
print(boruta3)

bor3 <- TentativeRoughFix(boruta3) #forcing the choice for remaining genes
print(bor3)

sel_var3 <- getSelectedAttributes(bor3) #selected variables
print(sel_var3)

#RUN4
boruta4 <- Boruta(class~., data=train_[folds$Fold4,], doTrace = 2, maxRuns = 500)
print(boruta4)

bor4 <- TentativeRoughFix(boruta4) #forcing the choice for remaining genes
print(bor4)

sel_var4 <- getSelectedAttributes(bor4) #selected variables
print(sel_var4)

#RUN5
boruta5 <- Boruta(class~., data=train_[folds$Fold5,], doTrace = 2, maxRuns = 500)
print(boruta5)

bor5 <- TentativeRoughFix(boruta5) #forcing the choice for remaining genes
print(bor5)

sel_var5 <- getSelectedAttributes(bor5) #selected variables
print(sel_var5)

sel_var <- Reduce(intersect, list(sel_var1,sel_var2,sel_var3,sel_var4,sel_var5)) #I selected only the genes common in all the runs
length(sel_var) #10 selected genes 

f <- as.formula(paste("class ~", paste(sel_var[!sel_var %in% "Class"], collapse = " + ")))

rfboruta1 <- randomForest(f, data=train_,proximity=TRUE) #new rf with selected genes
print(rfboruta1)

#test set accuracy
res_bor1 <- mod_res1(rfboruta1)

#I try with more genes. To do so, I select each gene that appear at least once in a run
sel_var <- c(sel_var1,sel_var2,sel_var3,sel_var4,sel_var5) 
sel_var <- unique(sel_var) #53 selected genes

f <- as.formula(paste("class ~", paste(sel_var[!sel_var %in% "Class"], collapse = " + ")))

rfboruta2 <- randomForest(f, data=train_,proximity=TRUE) #new rf with selected genes
print(rfboruta2)

#test set accuracy
res_bor2 <- mod_res1(rfboruta2)

#plot results
res1 <- data.frame(method = c('rf_all526', 'rf_122','rf_61','rf_40','rf_35','rf_30', 'rf_18','boruta+rf_10','boruta+rf_53'),
                   Accuracy = c(res_rf_all_var, res_rf_sel_var, res_rf_sel_var1, res_rf_sel_var2, res_rf_sel_var3, res_rf_sel_var4, res_rf_sel_var5, res_bor1, res_bor2),
                   se = c(0.0339, 0.0341, 0.0341, 0.0344, 0.0343, 0.0343, 0.0352, 0.0357, 0.0342)
                   )

# res1 <- within(res1, 
#                    method <- factor(method, 
#                                       levels=names(sort(table(method), 
#                                                         decreasing=FALSE))))

ggplot(res1, aes(x = reorder(method, Accuracy), y = Accuracy, color = method)) +
geom_errorbar(aes(ymin=Accuracy-se, ymax=Accuracy+se), width=.1)+
ggtitle("Feature Selection Performance Comparison in Test set")+
xlab('Method')+
theme(plot.title = element_text(lineheight=.8, face="bold"))+
geom_point(size = 4)+
theme(legend.position="none")

#########################################################################
#QUESTION3: Does the patient grouping based on similarity of mutated genes (with somatic mutations) 
#reflect the grouping based on cancer type?

#transposing mutation matrix
mut <- data[,3:ncol(data)]
mut <- t(mut)

#extracting cancer types
groups <- data$cancer_type

#Using non-negative matrix factorization (NMF) approach.
#I perform 10 runs for each value of k in range 2:24. I use 10 runs for computational
#reasons, usually better to use at least 30-50.
#Using Lee method since the default one fails in certain runs
estim.r <- nmf(mut, 2:24,method = 'lee', nrun = 10,.options = "t+v")

#print the various scores
plot(estim.r)

#Increasing the factorization rank would lead to decreasing residuals, as
#more variables are available to better fit the data - even on random data!
#In other words, there is potentially an overfitting problem.

#In this context, the approach from (Frigyesi et al. 2008) may be useful
#to prevent or detect overfitting as it takes into account the results for
#unstructured data.
#However it requires to compute the quality measure(s) for the random data. 
#NMF package package provides a function that shuffles the original data,
#by permuting the rows of each column, using each time a different permutation.

# shuffle original data
V.random <- randomize(mut)
# estimate quality measures from the shuffled data (use
# default NMF algorithm)
estim.r.random <- nmf(V.random, 2:24, method = 'lee', nrun = 10, .options = "t+v")
# plot measures on same graph
plot(estim.r, estim.r.random)

#to apply Frigyesi's approach, I compute the deltas between an rss and the 
#one in the following rank for both data bases

for (i in 1:length(estim.r$measures$rss-1)) {
  #print(i)
  result <- estim.r$measures$rss[i]-estim.r$measures$rss[i+1]
  print(result)
}

for (i in 1:length(estim.r.random$measures$rss-1)) {
  #print(i)
  result <- estim.r.random$measures$rss[i]-estim.r.random$measures$rss[i+1]
  print(result)
}

#The delta rss is still higher in the non-random data when k=24, therefore much
#higher than the real number of classes

#I try anyway the factorization with k=11 to assess how the obtained clusters 
#differ from the true classes
k11 <- nmf(mut, 11,method = 'lee', nrun = 30,.options = "t+v")

# get matrix W
w <- basis(k11)
dim(w)

# get matrix H
h <- coef(k11)
dim(h)

#I find each patient cluster relative to the NMF
clus <- c()
for (i in 1:ncol(h)) {
  #print(i)
  clus <- append(clus,which.max(h[,i]))
}

#I chose the BRCA cancer type to verify how patients with that
#cancer type were clustered
ind_BRCA <- which(groups=='BRCA') 
res <- clus[ind_BRCA]
table(res)

#note how BRCA patients were clusterized in 4 main clusters (1,4,8 and 10) 
#and in 7 minor clusters
