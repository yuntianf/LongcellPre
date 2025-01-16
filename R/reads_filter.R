#' @title frNN_dis
#'
#' @description Build frNN object from distance long table.
#' @details Build frNN object from distance long table for cluster by dbscan package.
#' @param dis A distance long table with 3 columns, each means node1, node2 and distance between them
#' @param n The total count of the nodes
#' @param eps The maximum threshold of distance to be preserved
#' @param if_direct If there is a direction between each nodes
#' @importFrom dbscan frNN
#' @return A frNN object
frNN_dis <- function(dis,n,eps = 10,if_direct = FALSE){
    dis <- as.matrix(dis)
    if(!if_direct){
        dis <- rbind(dis,dis[,c(2,1,3)])
    }
    dis <- as.data.frame(dis)
    colnames(dis) <- c("node1","node2","dist")
    dis <- dis[order(dis$node1),]

    out_frNN <- dbscan::frNN(as.dist(1), eps = eps)

    out_frNN$dist <- lapply(1:n,function(i){
      return(dis[dis$node1 == i,"dist"])
    })
    out_frNN$id <- lapply(1:n,function(i){
      return(dis[dis$node1 == i,"node2"])
    })
    out_frNN$sort <- F

    out_frNN <- dbscan::frNN(out_frNN,eps = eps)

    return(out_frNN)
}

#' @title isoform_dis_cluster
#'
#' @description Cluster isoforms using dbscan.
#' @details Cluster isoforms using dbscan.
#' @param isoforms A string vector storing isoforms
#' @param thresh The maxmimum threshold of difference for two isoforms to be clustered
#' @param eps The maximum threshold of distance to be preserved
#' @param split The character to split each exon in the isoform representation
#' @param sep The character to split the start and end for each exon in the isoform representation
#' @importFrom dbscan dbscan
#' @return A vector to record the cluster information of the isoforms
isoform_dis_cluster <- function(isoforms,thresh = 20,eps = 20,
                            split = "|",sep = ","){
    isoforms_count = table(isoforms)
    if(length(isoforms_count) == 1){
        return(rep(1,length(isoforms)))
    }

    dis = isos_dis(names(isoforms_count),
                   thresh = thresh,split = split,sep = sep)
    if(nrow(dis) == 0){
        cluster = 1:length(isoforms_count)
    }
    else{
        #dis[,3] = 1-dis[,3]
        dis_frNN = frNN_dis(dis,length(isoforms_count),eps = eps)
        cluster = dbscan(dis_frNN,minPts =1,weights = isoforms_count)$cluster
    }
    names(cluster) <- names(isoforms_count)
    cluster_expand = cluster[isoforms]
    return(cluster_expand)
}

#' @title size_filter_error
#'
#' @description Filter out singletons and small clusters which are not likely to be
#' correctly clustered, given the distribution of cluster size.
#' @param size A vector to store the size for each cluster.
#' @param ratio The estimated frequency for UMIs which are enriched of sequencing errors and
#' couldn't be clustered correctly and become singletons
#' @return A vector to indicating whether the cluster is filtered out or not
size_filter_error <- function(size,ratio = 0.1){
    size = as.numeric(size)
    if(sum(is.na(size))){
        stop("size must be integer!")
    }
    if(length(size) == 1){
        return(1)
    }
    if(max(size) == 1){
      return(rep(1,length(size)))
    }
    size = as.data.frame(cbind(1:length(size),size))
    colnames(size) = c("id","size")
    size = size[order(size$size),]

    size$weight = size_filter_cpp(size$size,ratio)
    size = size[order(size$id),]
    return(size$weight)
}


#' #' @title isoform_size_filter
#' #'
#' #' @description Filter out singletons and small clusters which are not likely to be
#' #' correctly clustered, given the distribution of cluster size.
#' #' @param size A vector to store the size for each cluster.
#' #' @param ratio The estimated frequency for UMIs which are enriched of sequencing errors and
#' #' couldn't be clustered correctly and become singletons
#' #' @return A numeric vector to indicating whether the cluster is filtered out or not
#' isoform_size_filter <- function(isoforms,size,ratio = 0.1){
#'   cluster = isoform_dis_cluster(isoforms,thresh = 10,eps = 10)
#'
#'   isoform_size = as.data.frame(cbind(cluster,size))
#'   colnames(isoform_size) = c("cluster","size")
#'
#'   isoform_size = isoform_size %>% group_by(cluster) %>%
#'     mutate(weight = size_filter_error(size,ratio))
#'   return(isoform_size$weight)
#' }


#' @title isoform_size_filter
#'
#' @description Filter out singletons and small clusters which are not likely to be
#' correctly clustered, given the distribution of cluster size within a series of isoform clusters.
#' @param isoforms A string vector, each string represent a corrected isoform by collapse of reads
#' from the same UMI
#' @inheritParams isoform_dis_cluster
#' @return A numeric vector to indicating whether the cluster is filtered out or not
isoform_size_filter <- function(isoforms,size,ratio = 0.1,...){
  if(length(isoforms) == 1){

  }
  len = isos_len(isoforms,...)
  cluster = isoform_dis_cluster(isoforms,thresh = 10,eps = 10)

  isoform_size = as.data.frame(cbind(cluster,size,len))
  colnames(isoform_size) = c("cluster","size","len")

  cluster_weight = isoform_size %>% group_by(cluster) %>%
                 mutate(weight = size_filter_error(size,ratio)) %>%
                 group_by(cluster) %>% summarise(weight = sum(weight),.groups = "drop")

  result <- isoform_size %>%
    mutate(id = 1:nrow(isoform_size)) %>%
    inner_join(cluster_weight, by = "cluster") %>%
    group_by(cluster) %>%
    arrange(-size,-len) %>%
    mutate(rank = row_number())%>%
    mutate(count = ifelse(rank <= weight,1,0)) %>%
    ungroup() %>% arrange(id)

    # arrange(desc(len), .by_group = TRUE) %>%
    # mutate(rank = row_number()) %>%
    # filter(rank <= weight) %>% # Keep only top n_longest per cluster
    # dplyr::select(-rank,-n_longest) %>%             # Drop the helper column
    # ungroup()

  return(result$count)
}

#' @title cells_isoforms_size_filter
#'
#' @description Apply cluster size filtering to each isoform in each cell.
#' @param cell_isoform_table A dataframe for reads including its cell source,
#' isoform sequence, splicing sites and umi cluster size.
#' @param ratio The estimated frequency for UMIs which are enriched of sequencing errors and
#' couldn't be clustered correctly and become singletons
#' @return A dataframe with isoform count after cluster size filtering for each cell
cells_isoforms_size_filter <- function(cell_isoform_table,ratio = 0.1){
  cell_isoform_table = cell_isoform_table %>%
                       group_by(cell,mid) %>%
                       mutate(weight = isoform_size_filter(isoform,size,ratio),
                              .groups = "drop")
  cell_isoform_table = cell_isoform_table %>%
                       group_by(cell,isoform) %>%
                       summarise(size = sum(size),cluster = n(),count = sum(weight),
                                 polyA = sum(polyA*weight)/sum(weight),
                                 .groups = "drop")
  return(cell_isoform_table)
}

# isoforms_size_filter <- function(isoform_table,relation,
#                                  alpha = 0.05,ratio = 0.2,
#                                  isoform = "isoform",mid = "mid",
#                                  size = "size",
#                                  polyA = "polyA"){
#     reads_filter <- lapply(unique(isoform_table[,mid]),function(i){
#         sub = isoform_table[isoform_table[,mid] == i,]
#         cluster = isoform_dis_cluster(sub[,isoform],thresh = 10,eps = 10)
#         weight = lapply(unique(cluster),function(x){
#             id = which(cluster == x)
#             sub_weight = 1 - size_filter(sub[id,size],relation = relation,
#                                      alpha = alpha, ratio = ratio)
#             return(cbind(id,sub_weight))
#         })
#         weight = as.data.frame(do.call(rbind,weight))
#         colnames(weight) = c("id","weight")
#         weight = weight[order(weight$id),]
#         sub$weight = weight$weight

#         sub = sub %>%
#               group_by(across(all_of(c(isoform)))) %>%
#               summarise(size = sum(across(all_of(size))),
#                         ### just for benchmark ###
#                         cluster = n(),
#                         ### just for benchmark ###
#                         count = sum(weight),
#                         polyA = mean(!! rlang::sym(polyA)),
#                        .groups = "drop")
#         #if(sum(sub$count == 0) > 0){
#         #    sub = sub[-which(sub$count == 0),]
#         #}
#         return(sub)
#     })
#     reads_filter <- as.data.frame(do.call(rbind,reads_filter))
#     return(reads_filter)
# }

# cells_isoforms_size_filter <- function(cell_isoform_table,
#                                        relation,alpha = 0.05,ratio = 0.2,
#                                        cell = "cell",isoform = "isoform",
#                                        mid = "mid",size = "size",
#                                        polyA = "polyA"){
#     reads_filter = lapply(unique(cell_isoform_table[,cell]),function(i){
#         sub = cell_isoform_table[cell_isoform_table[,cell] == i,]
#         filter = isoforms_size_filter(isoform_table = sub,relation = relation,
#                                       alpha = alpha,ratio = ratio,
#                                       mid = mid,isoform = isoform,
#                                       size = size,polyA = polyA)
#         filter = cbind(i,filter)
#         colnames(filter)[1] = cell
#         return(filter)
#     })
#     reads_filter = do.call(rbind,reads_filter)
#     return(reads_filter)
# }
