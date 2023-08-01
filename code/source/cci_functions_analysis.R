# functions for analysis of cell type interactomes

#-------------------------------------------------------------------------------


get_pca <- function(cci, cts = NULL, sample = NULL, variable = NULL){
  # cts = vector of specific cell types of interest
  # sample = name of a qualitative variable (e.g. Age, Species)
  # variable = qualitative variables of each sample (e.g. yng/old for Age)
  
  #subset and transpose, keep cts of interest
  if(is.null(cts) == TRUE){
    temp_df <- data.frame(t(cci$Score))
  }else{
    for(i in 1:length(cts))
      temp_df <- data.frame(t(cci$Score[,grep(cts[i], colnames(cci$Score))]))
  }
  
  # if no sample is given calculation is quite straight forward
  if(is.null(sample) == TRUE){
    
    print_pca <- FactoMineR::PCA(temp_df, graph = FALSE)
    
    # if sample is given, information needs to be added to the temp_df
  }else if(is.null(sample) == FALSE){
    
    #add an empty vector to contain samples or other variables
    for(j in 1:length(sample)){
      temp_df[,(ncol(temp_df)+1)] <- vector()
      colnames(temp_df)[ncol(temp_df)] <- sample[j]
      temp_df[,ncol(temp_df)] <- variable[[j]]
    }
    
    # indicate which columns of temp_df are just variables, not numbers
    quali_cols <- which(!is.na(match(colnames(temp_df), sample)))
    print_pca <- FactoMineR::PCA(temp_df, graph = FALSE, quali.sup = quali_cols)
    print(quali_cols)
  }
  
  return(print_pca)
}
