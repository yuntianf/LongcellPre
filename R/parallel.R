#' @title coreDetect
#'
#' @description Check if the running platform has required number cores
#' @details Check if the running platform has required number cores, if not will use the maximum
#'
#' @param cores Number of cores used to do parallization.
#' @importFrom parallel detectCores
coreDetect = function(cores){
  max_cores = detectCores()
  if(cores > max_cores){
    warning(paste(c("Distributed number of cores exceed the maximum, will use all ",
                    max_core, " cores for parallization!"),sep = ""))
    cores = max_cores
  }
  return(cores)
}

#' @title genes_distribute
#'
#' @description Distribute genes to different cores to make uniform mission distribution
#' @details Distribute genes to different cores given their number of reads to let each
#' core receive similar number of reads
#'
#' @param data The input dataframe, each row is a read from a gene.
#' @param cores The number of cores required for parallel
#' @param gene_col The name of the column in the input data which record the gene name
#' @import dplyr
#'
#' @return A list including split dataframe
genes_distribute = function(data,cores, gene_col = "gene"){
  if(cores == 1){
    return(list(data))
  }
  gene_read_num = data %>%
                  group_by(across(all_of(gene_col))) %>%
                  summarise(count = n()) %>% arrange(desc(count))
  colnames(gene_read_num)[colnames(gene_read_num) == gene_col] = "gene"

  unit = c(1:cores,cores:1)
  gene_read_num$group = rep(unit,(nrow(data) %/% length(unit) + 1))[1:nrow(data)]

  gene_list = (gene_read_num %>% group_by(group) %>% summarise(list(gene)))$gene

  data_split = lapply(gene_list,function(x){
    sub = data[data[,"gene_col"] %in% x,]
    return(sub)
  })
  return(data_split)
}

#' @title dataSplit
#'
#' @description Split a dataframe by row or a vector uniformly given the number of cores
#' @details Split a dataframe by row or a vector uniformly given the number of cores
#'
#' @param data An input dataframe or vector
#' @param cores The number of cores required for parallel
#' @return A list including split dataframe or vector
dataSplit = function(data,cores){
  cores = coreDetect(cores)

  num_splits <- cores
  if(is.vector(data)){
    group_size = ceiling(length(data)/num_splits)
    data_split = split(data,rep(1:num_splits, each = group_size, length.out = length(data)))
  }
  else if(is.data.frame(data)){
    group_size = ceiling(nrow(data)/num_splits)
    data_split <- split(data, rep(1:num_splits, each = group_size, length.out = nrow(data)))
  }
  else{
    stop("The input data should be a vector or a dataframe!")
  }
  return(data_split)
}
