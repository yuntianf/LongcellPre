#' @title saveResult
#'
#' @description Save the result into a file
#' @details Write the dataframe result in a txt table file
#'
#' @param data A dataframe
#' @param filename The output filename
saveResult = function(data,filename){
  write.table(data,file  = filename,sep = "\t",quote = FALSE,row.names = FALSE)
}


#' @title paramExtract
#'
#' @description Extract parameters from the ... for a function
#' @details Extract parameters from the ... for a function
#'
#' @param func functions to be match with parameters
#' @param ... Omitted parameters
#'
paramExtract = function(func,...){
  params = list(...)
  if(length(params) == 0){
    return(NULL)
  }
  
  arg.names = formalArgs(func)
  arg.names = arg.names[arg.names %in% names(params)]
  return(params[arg.names])
}

#' @title paramMerge
#'
#' @description Merge two sets of parameters
#' @details Merge two sets of parameters
#'
#' @inheritParams paramExtract
#' @param neceParam A list of user input parameters
paramMerge = function(func,neceParam,...){
  dotParam = paramExtract(func,...)
  if(length(dotParam) > 0 & length(neceParam) > 0){
    dotParam = dotParam[!names(dotParam) %in% names(neceParam)]
  }
  Param = c(neceParam,dotParam)
}

#' @title long2wide
#'
#' @description transform a long table into a wide matrix
#'
#' @param long The input long table
#' @param row_names_from The name of the column which would be transformed to each row index
#' @param col_names_from The name of the column which would be transformed to each col index
#' @param values_from The name of the column which would be set as values in the matrix
#' @param symmetric The flag to indicate if the transformed matrix is symmetric.
long2wide = function (long, row_names_from, col_names_from, values_from, 
          symmetric = TRUE) 
{
  long = long[, c(row_names_from, col_names_from, values_from)]
  if (symmetric) {
    rlong = long[, c(col_names_from, row_names_from, values_from)]
    colnames(rlong) = c(row_names_from, col_names_from, values_from)
    long = rbind(long, rlong)
    long = long[!duplicated(long), ]
  }
  out = pivot_wider(long, names_from = col_names_from, values_from = values_from)
  out = as.data.frame(out)
  rownames(out) = out[, row_names_from]
  out = out[,-which(colnames(out)==row_names_from)]
  if (symmetric) {
    out = out[rownames(out), rownames(out)]
  }
  return(out)
}

#' @title long2square
#'
#' @description transform a long table into a square matrix
#' @inheritParams long2wide
#' @param na.fill The value to fill the NA in the square matrix
#' @param nodes The pre-defined row and col index
long2square = function(long, row_names_from, col_names_from, values_from, 
                       symmetric = TRUE,na.fill = 0,nodes = NULL){
  long = long[, c(row_names_from, col_names_from, values_from)]
  if (symmetric) {
    rlong = long[, c(col_names_from, row_names_from, values_from)]
    colnames(rlong) = c(row_names_from, col_names_from, values_from)
    long = rbind(long, rlong)
    long = long[!duplicated(long), ]
  }
  if(is.null(nodes)){
    nodes = unique(unlist(long[,row_names_from],long[,col_names_from]))
  }
  mat = as.data.frame(tidyr::pivot_wider(long,
                                         names_from = col_names_from,
                                         values_from = values_from))
  rownames(mat) = mat[, row_names_from]
  mat = mat[,-which(colnames(mat)==row_names_from),drop = FALSE]
  #return(mat)
  diff = setdiff(nodes,colnames(mat))
  if(length(diff) > 0){
    mat[,diff] = NA
  }
  
  mat = mat[nodes,nodes,drop = FALSE]
  rownames(mat) = colnames(mat) = nodes
  
  mat[is.na(mat)] = na.fill
  mat = as.matrix(mat)
  
  return(mat)
}

#' @title save10X
#'
#' @description save a long table into cellRanger outpur format
#' @param long The input long table
#' @param path The output path
#' @param i The name of the column which stores the features
#' @param j The name of the column which stores the cells
#' @param value The name of the column which stores the feature count
#' @importFrom Matrix writeMM
save10X = function(long,path,i = "gene",j = "cell",value = "count"){
  long = long[,c(i,j,value)]
  
  x = names(table(long[,i]))
  y = names(table(long[,j]))
  
  long = long %>% mutate_at(c(i,j), ~as.numeric(as.factor(.)))
  long = as.data.frame(long)
  
  sparse_mat <- Matrix::sparseMatrix(
    i = long[,i],
    j = long[,j],
    x = long[,value]
  )
  
  write.table(x, file = file.path(path,"features.tsv"), sep = "\t", 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(y, file = file.path(path,"barcodes.tsv"), sep = "\t", 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  writeMM(sparse_mat, file = file.path(path,"matrix.mtx"))
}

#' @title saveIsoMat
#'
#' @description save the output from "cells_genes_isos_count" into 
#' cellRanger outpur format, for both gene quantification matrix and isoform quantification
#' matrix.
#' @param iso The output from "cells_genes_isos_count", which is a dataframe with four columns
#' @param path The output path
#' @param gene_col The name of the column which stores the gene name
#' @param cell_col The name of the column which stores the cells
#' @param iso_col The name of the column which stores the isoform name
#' @param count_col The name of the column which stores the isoform count
saveIsoMat = function(iso,path, cell_col = "cell",gene_col = "gene",
                      iso_col = "isoform",count_col= "count"){
  dir.create(file.path(path,"gene"),showWarnings = FALSE)
  dir.create(file.path(path,"isoform"),showWarnings = FALSE)
  
  iso = iso[,c(cell_col,gene_col,iso_col,count_col)]
  colnames(iso) = c("cell","gene","isoform","count")
  
  gene_long = iso %>% group_by(cell,gene) %>% 
              summarise(count = sum(count),.groups = "drop")
  iso_long = iso %>% filter(isoform != "unknown",count > 0) %>% dplyr::select(-gene)
  
  save10X(gene_long,file.path(path,"gene"))
  save10X(iso_long,file.path(path,"isoform"),i = "isoform")
}