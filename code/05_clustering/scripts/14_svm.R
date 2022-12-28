#-------------------------------------------------------------------------------

library(DropletUtils)
library(e1071)
# Gareth James Daniela Witten Trevor Hastie Robert Tibshirani An Introduction to Statistical Learning
library(caret)
# https://rstudio-pubs-static.s3.amazonaws.com/506713_93765c6d66074ed392e12c3339959cbf.html
#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_11"]])

sample_numbers <- sample(1:ncol(sce), as.integer(ncol(sce)*0.8))
sce_train <- sce[,sample_numbers]
sce_test <- sce[,-sample_numbers]

set.seed(1234)

#-------------------------------------------------------------------------------

# SVM on Hierarchical clustering
data_train <- as.data.frame(reducedDim(sce_train, type = "PCA"))
data_train$cluster_hierarchical <- sce_train$cluster_hierarchical

# tune includes 10-fold cross validation of cost, gamma and chooses the best model
tune_out <- tune(svm, factor(cluster_hierarchical) ~ ., 
                 data = data_train, kernel="radial",
                 ranges = list(cost =c(0.001, 0.01, 0.1, 1,5,10,100),
                               gamma = c(0.5,1,2,3,4)))
bestmod <- tune_out$best.model

data_test <- as.data.frame(reducedDim(sce_test, type = "PCA"))
data_test$cluster_hierarchical <- sce_test$cluster_hierarchical

pred <- predict(bestmod, data_test)
sce_test$prediction_hierarchical <- pred

object_list <- list(tune_out, pred)
saveRDS(object_list, file = snakemake@output[["objects_hierarchical"]])

#-------------------------------------------------------------------------------

# Seurat clustering
data_train <- as.data.frame(reducedDim(sce_train, type = "PCA"))
data_train$cluster_seurat <- sce_train$cluster_seurat

tune_out <- tune(svm, factor(cluster_seurat) ~ ., 
                 data = data_train, kernel="radial",
                 ranges = list(cost =c(0.001, 0.01, 0.1, 1,5,10,100),
                               gamma = c(0.5,1,2,3,4)))

bestmod <- tune_out$best.model

data_test <- as.data.frame(reducedDim(sce_test, type = "PCA"))
data_test$cluster_seurat <- sce_test$cluster_seurat

pred <- predict(bestmod, data_test)
sce_test$prediction_seurat <- pred

object_list <- list(tune_out, pred)
saveRDS(object_list, file = snakemake@output[["objects_seurat"]])

#-------------------------------------------------------------------------------

# Louvain clustering
data_train <- as.data.frame(reducedDim(sce_train, type = "PCA"))
data_train$cluster_louvain <- sce_train$cluster_louvain

tune_out <- tune(svm, factor(cluster_louvain) ~ ., 
                 data = data_train, kernel="radial",
                 ranges = list(cost =c(0.001, 0.01, 0.1, 1,5,10,100),
                               gamma = c(0.5,1,2,3,4)))

bestmod <- tune_out$best.model

data_test <- as.data.frame(reducedDim(sce_test, type = "PCA"))
data_test$cluster_louvain <- sce_test$cluster_louvain

pred <- predict(bestmod, data_test)
sce_test$prediction_louvain <- pred

object_list <- list(tune_out, pred)
saveRDS(object_list, file = snakemake@output[["objects_louvain"]])

saveRDS(sce_test, file = snakemake@output[["sce_test"]])
