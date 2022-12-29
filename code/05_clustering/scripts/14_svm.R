#-------------------------------------------------------------------------------

library(DropletUtils)
library(caret)
# https://rstudio-pubs-static.s3.amazonaws.com/506713_93765c6d66074ed392e12c3339959cbf.html
library(e1071)
# Gareth James Daniela Witten Trevor Hastie Robert Tibshirani An Introduction to Statistical Learning

#-------------------------------------------------------------------------------

sce <- readRDS(file = snakemake@input[["sce_13"]])

# Reproducible random numbers (https://bookdown.org/rdpeng/rprogdatascience/parallel-computation.html#example-bootstrapping-a-statistic)
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
sample_numbers <- lapply(1, function(i){
  sample_numbers <- sample(1:ncol(sce), as.integer(ncol(sce)*0.8))
  return(sample_numbers)
})

sce_train <- sce[,sample_numbers[[1]]]
sce_test <- sce[,-sample_numbers[[1]]]
use_kernel <- snakemake@params[["use_kernel"]]

if(use_kernel == "linear"){
  ranges_list <- list(cost =c(0.001, 0.01, 0.1, 1, 2, 5, 10))
}else if(use_kernel == "radial"){
  ranges_list <- list(cost =c(0.001, 0.01, 0.1, 1, 2, 5, 10),
                      gamma = c(0.5,1,2,3,4))
}

# make weight lists for each cluster in the training set because the classes 
# (clusters) are very imbalanced 
# https://www.analyticsvidhya.com/blog/2020/10/improve-class-imbalance-class-weights/
# method that provides the same values*10: https://medium.com/grabngoinfo/balanced-weights-for-imbalanced-classification-465f0e13c5ad
clusters_hierarchical <- unfactor(unique(sce_train$cluster_hierarchical))
nr_classes <- length(clusters_hierarchical)
weight_h <- vector()
for(i in 1:length(clusters_hierarchical)){
  nr_samples <- length(grep(clusters_hierarchical[i], unfactor(sce_train$cluster_hierarchical)))
  weight_h[clusters_hierarchical[i]] <- ncol(sce_train)/(nr_samples*nr_classes)
}

clusters_seurat <- unfactor(unique(sce_train$cluster_seurat))
nr_classes <- length(clusters_seurat)
weight_s <- vector()
for(i in 1:length(clusters_seurat)){
  nr_samples <- length(grep(clusters_seurat[i], unfactor(sce_train$cluster_seurat)))
  weight_s[clusters_seurat[i]] <- ncol(sce_train)/(nr_samples*nr_classes)
}

clusters_louvain <- unfactor(unique(sce_train$cluster_louvain))
nr_classes <- length(clusters_louvain)
weight_l <- vector()
for(i in 1:length(clusters_louvain)){
  nr_samples <- length(grep(clusters_louvain[i], unfactor(sce_train$cluster_louvain)))
  weight_l[clusters_louvain[i]] <- ncol(sce_train)/(nr_samples*nr_classes)
}

#-------------------------------------------------------------------------------

# SVM on Hierarchical clustering
data_train <- as.data.frame(reducedDim(sce_train, type = "PCA"))
data_train$cluster_hierarchical <- sce_train$cluster_hierarchical

# includes 10-fold cross validation of cost, gamma and chooses the best model
tune_out <- tune(svm, factor(cluster_hierarchical) ~ ., 
                 data = data_train, kernel=use_kernel,
                 ranges = ranges_list,
                 class.weights = weight_h)
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
                 data = data_train, kernel=use_kernel,
                 ranges = ranges_list,
                 class.weights = weight_s)

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
                 data = data_train, kernel=use_kernel,
                 ranges = ranges_list,
                 class.weights = weight_l)

bestmod <- tune_out$best.model

data_test <- as.data.frame(reducedDim(sce_test, type = "PCA"))
data_test$cluster_louvain <- sce_test$cluster_louvain

pred <- predict(bestmod, data_test)
sce_test$prediction_louvain <- pred

object_list <- list(tune_out, pred)
saveRDS(object_list, file = snakemake@output[["objects_louvain"]])

saveRDS(sce_test, file = snakemake@output[["sce_test"]])
