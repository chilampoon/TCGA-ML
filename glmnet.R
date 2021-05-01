##Classification of Multiple Myeloma sub-challenge 2
# on microarray and RNA-Seq data with glmnet

##1. Load datasets and pre-process data
base.dataset <- readRDS("sc2scaleddataset1.rds")
# relevel to FALSE to match probability scores
base.dataset$HR_FLAG <- factor(base.dataset$HR_FLAG)
base.dataset$HR_FLAG <- relevel(base.dataset$HR_FLAG, ref="FALSE")
base.dataset <- base.dataset[,c(1:3, 7:9, 28:ncol(base.dataset))]
base.dataset <- na.omit(base.dataset)

# filter microarray or RNA-Seq datasets
# library(dplyr)
# base.dataset <- filter(base.dataset, Study=="MMRF")

# # create new HR_FLAGs based on D_PFS
# HR_FLAG22 <- (base.dataset$D_PFS_FLAG==1 & base.dataset$D_PFS<=671)
# base.dataset <- cbind.data.frame(base.dataset, HR_FLAG22)
# base.dataset$HR_FLAG <- base.dataset$HR_FLAG22
# base.dataset$HR_FLAG22 <- NULL

# split into 90% training and 10% testing
library(caTools)
set.seed(77)
sample <- sample.split(base.dataset$Patient, SplitRatio=.9)
base.dataset.train <- subset(base.dataset, sample==TRUE)
base.dataset.test <- subset(base.dataset, sample==FALSE)

# prepare subsets of different datasets
base.dataset.test.EMTAB <- filter(base.dataset.test, Study=="EMTAB4032")
base.dataset.test.GSE <- filter(base.dataset.test, Study=="GSE24080UAMS")
base.dataset.test.MMRF <- filter(base.dataset.test, Study=="MMRF")

# random undersampling for low risk flags in base.dataset.train
# first split HR and LR subsets
library(dplyr)
base.dataset.train.HR <- filter(base.dataset.train, HR_FLAG=="TRUE")
base.dataset.train.LR <- filter(base.dataset.train, HR_FLAG=="FALSE")
set.seed(77)
sample <- sample.split(base.dataset.train.LR$Patient, SplitRatio=nrow(base.dataset.train.HR)/nrow(base.dataset.train.LR))
base.dataset.train.LR <- subset(base.dataset.train.LR, sample==TRUE)
base.dataset.train <- rbind(base.dataset.train.HR, base.dataset.train.LR)

# prepare labels
hrFLAG <- base.dataset.train$HR_FLAG
# hrFLAG <- relevel(hrFLAG, ref="TRUE") #relevel labels to TRUE as 1
pfs <- base.dataset.train$D_PFS
pfsFLAG <- base.dataset.train$D_PFS_FLAG
labels <- cbind.data.frame(hrFLAG, pfs, pfsFLAG)
# create subset of relevant features for training
dataset <- base.dataset.train[,c(3, 6:ncol(base.dataset.train))]

# Loop to progressively reduce number of variance as feature vectors
library(caret)
library(glmnet)
library(FSelector)
library(sda)
source("RScripts/classificationFunctions.R")
source("RScripts/metrics.R")

# create 5 folds
folds <- createFolds(labels$hrFLAG, k=5)
results.list <- list() # create list to store results
metrics.list <- list()
parameters.list <- list() # list to store best parameters
EMTAB.list <- list()
GSE.list <- list()
MMRF.list <- list()

# 5-fold cross-validation
# outer loop of cross-validation
# inner loop of feature selection and model building
for (i in 1:5) {
  results.list[[i]] <- list()
  metrics.list[[i]] <- list()
  parameters.list[[i]] <- list()
  EMTAB.list[[i]] <- list()
  GSE.list[[i]] <- list()
  MMRF.list[[i]] <- list()
  f <- 1793 #number of features
  
  print("Preparing folds...")
  # prepare data according to folds
  # 1/5 of the dataset as validation, remaining data for training
  cv.train.base <- dataset[-folds[[i]],]
  cv.trainLabels <- labels[-folds[[i]],]
  cv.train.clinical <- dataset[-folds[[i]],1:2]
  cv.validate.base <- dataset[folds[[i]],]
  cv.validateLabels <- labels[folds[[i]],]
  cv.validate.clinical <- dataset[folds[[i]],1:2]
  
  print("FS filtering...")
  # feature selection only on the genes
  weights <- symmetrical.uncertainty(HR_FLAG~., cv.train.base[,3:ncol(cv.train.base)])
  # ranking.lda <- sda.ranking(data.matrix(cv.train.base[,3:(ncol(cv.train.base)-1)]), cv.trainLabels$hrFLAG, diagonal=FALSE)
  
  # initialize variables to store best set of features achieving the best results
  parameters.list[[i]]$f <- f
  parameters.list[[i]]$fiAUC <- 0
  
  # first inner FS loop
  # filter method using symmetrical uncertainty
  # start from 1793 features onwards
  while(f > 2) {
    cat(paste0("Fold ", i, " Features: ", f, "\n"))
    results.list[[i]][[as.character(f)]] <- list()
    metrics.list[[i]][[as.character(f)]] <- list()
    EMTAB.list[[i]][[as.character(f)]] <- list()
    GSE.list[[i]][[as.character(f)]] <- list()
    MMRF.list[[i]][[as.character(f)]] <- list()
    
    print("Preparing dataset...")
    # Filter out genes with top f features
    weights.subset <- cutoff.k(weights,f)
    # weights.subset <- rownames(ranking.lda[1:f,])
    cv.train <- cv.train.base[,weights.subset]
    # cv.train <- data.matrix(cv.train)
    # # merge back with clinical data
    cv.train <- data.matrix(cbind(cv.train.clinical, cv.train))
    cv.validate <- cv.validate.base[,weights.subset]
    # cv.validate <- data.matrix(cv.validate)
    cv.validate <- data.matrix(cbind(cv.validate.clinical, cv.validate))
    features <- colnames(cv.validate)

    print("Building model...")
    ##BUILD MODELS##
    # fit model using regularized logistic regression
    cat("glmnet\n")
    cv.fit.lasso <- cv.glmnet(x=cv.train, y=cv.trainLabels$hrFLAG, family="binomial")
    predictionLabels <- predict(cv.fit.lasso, newx=cv.validate, s='lambda.min', type="class")
    predictionDecisions <- predict(cv.fit.lasso, newx=cv.validate, s='lambda.min', type="response")
    results.list[[i]][[as.character(f)]][['glmnet']] <- perf.eval(cv.validateLabels$hrFLAG, predictionLabels, predictionDecisions, cv.validateLabels$hrFLAG==TRUE, decreasing=TRUE)
    metrics.list[[i]][[as.character(f)]][['glmnet']] <- calculate.metrics(as.vector(predictionDecisions), as.numeric(as.logical(predictionLabels)), cv.validateLabels$pfs, cv.validateLabels$pfsFLAG)
    print(paste0("AUC ", results.list[[i]][[as.character(f)]][['glmnet']]$AUC))
    print(paste0("iAUC: ", metrics.list[[i]][[as.character(f)]][['glmnet']]$iAUC))
    print(paste0("BAC: ", metrics.list[[i]][[as.character(f)]][['glmnet']]$BAC))
    
    # split into different datasets
    # prepare subsets for testing
    testing.subset.EMTAB <- data.matrix(base.dataset.test.EMTAB[,features])
    # test for EMTAB
    predictionLabels.glmnet.test <- predict(cv.fit.lasso, newx=testing.subset.EMTAB, s='lambda.min', type="class")
    predictionScores.glmnet.test <- predict(cv.fit.lasso, newx=testing.subset.EMTAB, s='lambda.min', type="response")
    EMTAB.list[[i]][[as.character(f)]][['glmnet']] <- calculate.metrics(as.vector(predictionScores.glmnet.test), as.numeric(as.logical(predictionLabels.glmnet.test)), base.dataset.test.EMTAB$D_PFS, base.dataset.test.EMTAB$D_PFS_FLAG)
    print(paste0("EMTAB Test iAUC ",  EMTAB.list[[i]][[as.character(f)]][['glmnet']]$iAUC))
    print(paste0("EMTAB Test BAC ",  EMTAB.list[[i]][[as.character(f)]][['glmnet']]$BAC))
    
    testing.subset.GSE <- data.matrix(base.dataset.test.GSE[,features])
    # test for GSE
    predictionLabels.glmnet.test <- predict(cv.fit.lasso, newx=testing.subset.GSE, s='lambda.min', type="class")
    predictionScores.glmnet.test <- predict(cv.fit.lasso, newx=testing.subset.GSE, s='lambda.min', type="response")
    GSE.list[[i]][[as.character(f)]][['glmnet']] <- calculate.metrics(as.vector(predictionScores.glmnet.test), as.numeric(as.logical(predictionLabels.glmnet.test)), base.dataset.test.GSE$D_PFS, base.dataset.test.GSE$D_PFS_FLAG)
    print(paste0("GSE Test iAUC ", GSE.list[[i]][[as.character(f)]][['glmnet']]$iAUC))
    print(paste0("GSE Test BAC ", GSE.list[[i]][[as.character(f)]][['glmnet']]$BAC))
    
    testing.subset.MMRF <- data.matrix(base.dataset.test.MMRF[,features])
    # test for MMRF
    predictionLabels.glmnet.test <- predict(cv.fit.lasso, newx=testing.subset.MMRF, s='lambda.min', type="class")
    predictionScores.glmnet.test <- predict(cv.fit.lasso, newx=testing.subset.MMRF, s='lambda.min', type="response")
    MMRF.list[[i]][[as.character(f)]][['glmnet']] <- calculate.metrics(as.vector(predictionScores.glmnet.test), as.numeric(as.logical(predictionLabels.glmnet.test)), base.dataset.test.MMRF$D_PFS, base.dataset.test.MMRF$D_PFS_FLAG)
    print(paste0("MMRF Test iAUC ", MMRF.list[[i]][[as.character(f)]][['glmnet']]$iAUC))
    print(paste0("MMRF Test BAC ", MMRF.list[[i]][[as.character(f)]][['glmnet']]$BAC))
  
    # update results best results
    # save best model
    if(metrics.list[[i]][[as.character(f)]][['glmnet']]$iAUC > parameters.list[[i]]$fiAUC) {
      parameters.list[[i]]$fiAUC <- metrics.list[[i]][[as.character(f)]][['glmnet']]$iAUC
      parameters.list[[i]]$f <- f
      results.list[[i]][['0']][['glmnet']] <- perf.eval(cv.validateLabels$hrFLAG, predictionLabels, predictionDecisions, cv.validateLabels$hrFLAG==TRUE, decreasing=TRUE)
      metrics.list[[i]][['0']][['glmnet']] <- calculate.metrics(as.vector(predictionDecisions), as.numeric(as.logical(predictionLabels)), cv.validateLabels$pfs, cv.validateLabels$pfsFLAG)
      parameters.list[[i]]$model <- cv.fit.lasso
      parameters.list[[i]]$features <- colnames(cv.train)
    }
    
    # update index for the loop
    print("Updating index...")
    f = round(f/1.5)
  }
  # print best performance for each fold
  print(paste0("FOLD ", i, " symmetrical uncertainty"))
  print(paste0("Best iAUC: ", parameters.list[[i]]$fiAUC, " at ", parameters.list[[i]]$f, " features"))
  
}

# test models on training dataset and testing dataset
# using best models from each fold
test.list <- list()

for (i in 1:5) {
  test.list[[i]] <- list()
  model <- parameters.list[[i]]$model
  
  # prepare datasets for testing
  training.subset <- data.matrix(base.dataset.train[,parameters.list[[i]]$features])
  training.subset.labels <- factor(base.dataset.train$HR_FLAG)
  training.subset.pfs <- base.dataset.train$D_PFS
  training.subset.pfsFLAG <- base.dataset.train$D_PFS_FLAG
  testing.subset <- data.matrix(base.dataset.test[,parameters.list[[i]]$features])
  testing.subset.labels <- factor(base.dataset.test$HR_FLAG)
  testing.subset.pfs <- base.dataset.test$D_PFS
  testing.subset.pfsFLAG <- base.dataset.test$D_PFS_FLAG
  
  # test on training dataset first
  predictionLabels.train <- predict(model, newx=data.matrix(training.subset), type="class")
  predictionDecisions.train <- predict(model, newx=data.matrix(training.subset), type="response")
  test.list[[i]][['train']] <- calculate.metrics(as.vector(predictionDecisions.train), as.numeric(as.logical(predictionLabels.train)), training.subset.pfs, training.subset.pfsFLAG)
  print(paste0("Train iAUC ", test.list[[i]][['train']]$iAUC))
  
  # then test on independent test set
  predictionLabels.test <- predict(model, newx=data.matrix(testing.subset), type="class")
  predictionDecisions.test <- predict(model, newx=data.matrix(testing.subset), type="response")
  test.list[[i]][['test']] <- calculate.metrics(as.vector(predictionDecisions.test), as.numeric(as.logical(predictionLabels.test)), testing.subset.pfs, testing.subset.pfsFLAG)
  print(paste0("Test iAUC ", test.list[[i]][['test']]$iAUC))
}



saveRDS(results.list, file="RObjects/RNAscaledEXPRAGEISSresultslist-glmnet_su.rds")
saveRDS(metrics.list, file="RObjects/RNAscaledEXPRAGEISSmetricslist-glmnet_su.rds")
saveRDS(parameters.list, file="RObjects/RNAscaledEXPRAGEISSparameterslist-glmnet_su.rds")
saveRDS(test.list, file="RObjects/RNAscaledEXPRAGEISStestlist-glmnet_su.rds")

##PRINT RESULTS##
pdf(file="Results/RNAscaledEXPRAGEISS-glmnet_su.pdf")
## best features and optimized parameters at each fold
par(mfrow=c(1,4))
#BAC
boxplot(unlist(lapply(lapply(lapply(metrics.list, '[[', '0'), '[[', "glmnet"), '[[', 'BAC')),
        main="BAC" , medcol=1, ylim=c(0,1))
abline(h=0.5, col="gray20", lty=2)

#MCC
boxplot(unlist(lapply(lapply(lapply(metrics.list, '[[', '0'), '[[', "glmnet"), '[[', 'MCC')),
        main="MCC", medcol=1, ylim=c(- 1,1))
abline(h=0, col="gray20", lty=2)

#AUC
boxplot(unlist(lapply(lapply(lapply(metrics.list, '[[', '0'), '[[', "glmnet"), '[[', 'AUC')),
        main="AUC", medcol=1, ylim=c(0,1))
abline(h=0.5, col="gray20", lty=2)

#iAUC
boxplot(unlist(lapply(lapply(lapply(metrics.list, '[[', '0'), '[[', "glmnet"), '[[', 'iAUC')),
        main="iAUC", medcol=1, ylim=c(0,1))
abline(h=0.5, col="gray20", lty=2)


## results of each fold at different features
par(mfrow=c(1,2))
# Plot individual BAC for each fold
for(i in 1:length(metrics.list)) {
  if(i==1) {
    plot(as.integer(names(metrics.list[[1]])), unlist(lapply(lapply(metrics.list[[i]], '[[', "glmnet"), '[[', 'BAC')), log="x", ty='l', col=1, main="BAC", ylim=c(0,1), xlab='number of features', ylab='accuracy')
  } else {
    lines(as.integer(names(metrics.list[[1]])), unlist(lapply(lapply(metrics.list[[i]], '[[', "glmnet"), '[[', 'BAC')), ty='l', col=1)
  }
}
abline(h=0.5, col="gray20", lty=2)

# Plot individual MCC
for(i in 1:length(metrics.list)) {
  if(i==1) {
    plot(as.integer(names(metrics.list[[1]])), unlist(lapply(lapply(metrics.list[[i]], '[[', "glmnet"), '[[', 'MCC')), log="x", ty='l', col=1, main="MCC", ylim=c(-1,1), xlab='number of features', ylab='MCC')
  } else {
    lines(as.integer(names(metrics.list[[1]])), unlist(lapply(lapply(metrics.list[[i]], '[[', "glmnet"), '[[', 'MCC')), ty='l', col=1)
  }
}
abline(h=0, col="gray20", lty=2)

# Plot individual AUC
for(i in 1:length(metrics.list)) {
  if(i==1) {
    plot(as.integer(names(metrics.list[[1]])), unlist(lapply(lapply(metrics.list[[i]], '[[', "glmnet"), '[[', 'AUC')), log="x", ty='l', col=1, main="AUC", ylim=c(0,1), xlab='number of features', ylab='AUC')
  } else {
    lines(as.integer(names(metrics.list[[1]])), unlist(lapply(lapply(metrics.list[[i]], '[[', "glmnet"), '[[', 'AUC')), ty='l', col=1)
  }
}
abline(h=0.5, col="gray20", lty=2)

# Plot individual iAUC
for(i in 1:length(metrics.list)) {
  if(i==1) {
    plot(as.integer(names(metrics.list[[1]])), unlist(lapply(lapply(metrics.list[[i]], '[[', "glmnet"), '[[', 'iAUC')), log="x", ty='l', col=1, main="AUPRC", ylim=c(0,1), xlab='number of features', ylab='AUPRC')
  } else {
    lines(as.integer(names(metrics.list[[1]])), unlist(lapply(lapply(metrics.list[[i]], '[[', "glmnet"), '[[', 'iAUC')), ty='l', col=1)
  }
}
abline(h=0.5, col="gray20", lty=2)

# plot individual box plots of AUC, MCC and AUC for each classifier at each feature vector
par(mfrow=c(1,2))
#BAC
boxplot(unlist(lapply(lapply(lapply(metrics.list, '[[', '0'), '[[', "glmnet"), '[[', 'BAC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '3'), '[[', "glmnet"), '[[', 'BAC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '4'), '[[', "glmnet"), '[[', 'BAC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '6'), '[[', "glmnet"), '[[', 'BAC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '9'), '[[', "glmnet"), '[[', 'BAC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '14'), '[[', "glmnet"), '[[', 'BAC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '21'), '[[', "glmnet"), '[[', 'BAC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '31'), '[[', "glmnet"), '[[', 'BAC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '47'), '[[', "glmnet"), '[[', 'BAC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '70'), '[[', "glmnet"), '[[', 'BAC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '105'), '[[', "glmnet"), '[[', 'BAC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '157'), '[[', "glmnet"), '[[', 'BAC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '236'), '[[', "glmnet"), '[[', 'BAC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '354'), '[[', "glmnet"), '[[', 'BAC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '531'), '[[', "glmnet"), '[[', 'BAC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '797'), '[[', "glmnet"), '[[', 'BAC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '1195'), '[[', "glmnet"), '[[', 'BAC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '1793'), '[[', "glmnet"), '[[', 'BAC')),
        main="BAC", names=c(0, 3, 4, 6, 9, 14, 21, 31, 47, 70, 105, 157, 236, 354, 531, 797, 1195, 1793), medcol=1, ylim=c(0,1))
abline(h=0.5, col="gray20", lty=2)

#MCC
boxplot(unlist(lapply(lapply(lapply(metrics.list, '[[', '0'), '[[', "glmnet"), '[[', 'MCC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '3'), '[[', "glmnet"), '[[', 'MCC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '4'), '[[', "glmnet"), '[[', 'MCC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '6'), '[[', "glmnet"), '[[', 'MCC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '9'), '[[', "glmnet"), '[[', 'MCC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '14'), '[[', "glmnet"), '[[', 'MCC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '21'), '[[', "glmnet"), '[[', 'MCC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '31'), '[[', "glmnet"), '[[', 'MCC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '47'), '[[', "glmnet"), '[[', 'MCC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '70'), '[[', "glmnet"), '[[', 'MCC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '105'), '[[', "glmnet"), '[[', 'MCC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '157'), '[[', "glmnet"), '[[', 'MCC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '236'), '[[', "glmnet"), '[[', 'MCC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '354'), '[[', "glmnet"), '[[', 'MCC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '531'), '[[', "glmnet"), '[[', 'MCC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '797'), '[[', "glmnet"), '[[', 'MCC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '1195'), '[[', "glmnet"), '[[', 'MCC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '1793'), '[[', "glmnet"), '[[', 'MCC')),
        main="MCC", names=c(0, 3, 4, 6, 9, 14, 21, 31, 47, 70, 105, 157, 236, 354, 531, 797, 1195, 1793), medcol=1, ylim=c(-1,1))
abline(h=0, col="gray20", lty=2)

#AUC
boxplot(unlist(lapply(lapply(lapply(metrics.list, '[[', '0'), '[[', "glmnet"), '[[', 'AUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '3'), '[[', "glmnet"), '[[', 'AUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '4'), '[[', "glmnet"), '[[', 'AUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '6'), '[[', "glmnet"), '[[', 'AUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '9'), '[[', "glmnet"), '[[', 'AUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '14'), '[[', "glmnet"), '[[', 'AUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '21'), '[[', "glmnet"), '[[', 'AUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '31'), '[[', "glmnet"), '[[', 'AUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '47'), '[[', "glmnet"), '[[', 'AUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '70'), '[[', "glmnet"), '[[', 'AUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '105'), '[[', "glmnet"), '[[', 'AUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '157'), '[[', "glmnet"), '[[', 'AUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '236'), '[[', "glmnet"), '[[', 'AUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '354'), '[[', "glmnet"), '[[', 'AUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '531'), '[[', "glmnet"), '[[', 'AUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '797'), '[[', "glmnet"), '[[', 'AUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '1195'), '[[', "glmnet"), '[[', 'AUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '1793'), '[[', "glmnet"), '[[', 'AUC')),
        main="AUC", names=c(0, 3, 4, 6, 9, 14, 21, 31, 47, 70, 105, 157, 236, 354, 531, 797, 1195, 1793), medcol=1, ylim=c(0,1))
abline(h=0.5, col="gray20", lty=2)

#iAUC
boxplot(unlist(lapply(lapply(lapply(metrics.list, '[[', '0'), '[[', "glmnet"), '[[', 'iAUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '3'), '[[', "glmnet"), '[[', 'iAUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '4'), '[[', "glmnet"), '[[', 'iAUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '6'), '[[', "glmnet"), '[[', 'iAUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '9'), '[[', "glmnet"), '[[', 'iAUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '14'), '[[', "glmnet"), '[[', 'iAUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '21'), '[[', "glmnet"), '[[', 'iAUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '31'), '[[', "glmnet"), '[[', 'iAUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '47'), '[[', "glmnet"), '[[', 'iAUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '70'), '[[', "glmnet"), '[[', 'iAUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '105'), '[[', "glmnet"), '[[', 'iAUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '157'), '[[', "glmnet"), '[[', 'iAUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '236'), '[[', "glmnet"), '[[', 'iAUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '354'), '[[', "glmnet"), '[[', 'iAUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '531'), '[[', "glmnet"), '[[', 'iAUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '797'), '[[', "glmnet"), '[[', 'iAUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '1195'), '[[', "glmnet"), '[[', 'iAUC')),
        unlist(lapply(lapply(lapply(metrics.list, '[[', '1793'), '[[', "glmnet"), '[[', 'iAUC')),
        main="iAUC", names=c(0, 3, 4, 6, 9, 14, 21, 31, 47, 70, 105, 157, 236, 354, 531, 797, 1195, 1793), medcol=1, ylim=c(0,1))
abline(h=0.5, col="gray20", lty=2)

par(mfrow=c(2,4))
# ROC curve
for(j in 1:length(results.list[[1]])) {
  # loop k for each fold
  # ROC curve
  for(k in 1:length(results.list)) {
    if(k==1) {
      plot(results.list[[k]][[names(results.list[[k]][j])]][["glmnet"]]$FPR, 
           results.list[[k]][[names(results.list[[k]][j])]][["glmnet"]]$TPR, 
           ty='l', col=1, main=names(results.list[[1]])[j], xlab='FPR',ylab='TPR', ylim=c(0,1))
    }
    else {
      lines(results.list[[k]][[names(results.list[[k]][j])]][["glmnet"]]$FPR, 
            results.list[[k]][[names(results.list[[k]][j])]][["glmnet"]]$TPR, 
            ty='l', col=1, main=names(results.list[[1]])[j], xlab='FPR',ylab='TPR', ylim=c(0,1))
    }
  }
  abline(0,1, col="gray20", lty=2)
}

par(mfrow=c(2,4))
#PR curve
for(j in 1:length(results.list[[1]])) {
  # loop k for each fold
  # PR curve
  for(k in 1:length(results.list)) {
    if(k==1) {
      plot(results.list[[k]][[names(results.list[[k]][j])]][["glmnet"]]$TPR, 
           results.list[[k]][[names(results.list[[k]][j])]][["glmnet"]]$precision, 
           ty='l', col=1, main=names(results.list[[1]])[j], xlab='recall (TPR)',ylab='precision', ylim=c(0,1))
    }
    else {
      lines(results.list[[k]][[names(results.list[[k]][j])]][["glmnet"]]$TPR, 
            results.list[[k]][[names(results.list[[k]][j])]][["glmnet"]]$precision, 
            ty='l', col=1, main=names(results.list[[1]])[j], xlab='recall (TPR)',ylab='precision', ylim=c(0,1))
    }
  }
  abline(h=0.3, col="gray20", lty=2)
}

# plot results with best model in each fold
# tested with training dataset
# and independent testing dataset

par(mfrow=c(1,2))
#ACC
barplot(rbind(unlist(lapply(lapply(test.list, '[[', "train"), '[[', 'BAC')),
              unlist(lapply(lapply(test.list, '[[', "test"), '[[', 'BAC'))),
        main="BAC" , col=c(4:5), ylim=c(0,1), beside=TRUE, axes=TRUE,
        legend.text=c("train", "test"))
abline(h=0.5, col="gray20", lty=2)

#MCC
barplot(rbind(unlist(lapply(lapply(test.list, '[[', "train"), '[[', 'MCC')),
              unlist(lapply(lapply(test.list, '[[', "test"), '[[', 'MCC'))),
        main="MCC" , col=c(4:5), ylim=c(-1,1), beside=TRUE, axes=TRUE,
        legend.text=c("train", "test"))
abline(h=0, col="gray20", lty=2)

#AUC
barplot(rbind(unlist(lapply(lapply(test.list, '[[', "train"), '[[', 'AUC')),
              unlist(lapply(lapply(test.list, '[[', "test"), '[[', 'AUC'))),
        main="AUC" , col=c(4:5), ylim=c(0,1), beside=TRUE, axes=TRUE,
        legend.text=c("train", "test"))
abline(h=0.5, col="gray20", lty=2)

#iAUC
barplot(rbind(unlist(lapply(lapply(test.list, '[[', "train"), '[[', 'iAUC')),
              unlist(lapply(lapply(test.list, '[[', "test"), '[[', 'iAUC'))),
        main="iAUC" , col=c(4:5), ylim=c(0,1), beside=TRUE, axes=TRUE,
        legend.text=c("train", "test"))
abline(h=0.5, col="gray20", lty=2)

dev.off()

