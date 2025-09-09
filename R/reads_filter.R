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
size_filter_error <- function(size, ratio = 0.1) {
  size <- as.numeric(size)
  if (any(is.na(size))) stop("size must be integer!")
  n <- length(size)
  if (n <= 1L) return(1)
  if (max(size) == 1) return(rep(1, n))

  ord <- order(size)                   # ascending, like your code
  w_sorted <- size_filter_cpp(size[ord], ratio)  # your existing C++ fn
  w <- numeric(n)
  w[ord] <- w_sorted
  w
}


#' @title isoform_size_filter
#'
#' @description Filter out singletons and small clusters which are not likely to be
#' correctly clustered, given the distribution of cluster size within a series of isoform clusters.
#' @param isoforms A string vector, each string represent a corrected isoform by collapse of reads
#' from the same UMI
#' @inheritParams isoform_dis_cluster
#' @importFrom dplyr row_number
#' @return A numeric vector to indicating whether the cluster is filtered out or not
isoform_size_filter <- function(isoforms, size, ratio = 0.1, ...,
                                   thresh = 10, eps = 10) {
  n <- length(isoforms)
  if (n == 0L) return(integer())

  len     <- isos_len_cpp(isoforms)                      # now fast (C++)
  cluster <- isoform_dis_cluster(isoforms, thresh, eps)   # micro-optimized

  dt <- data.table(cluster = cluster,
                   size    = as.numeric(size),
                   len     = as.numeric(len),
                   ord     = seq_len(n))

  cw <- dt[, .(weight = sum(size_filter_error(size, ratio))), by = cluster]

  setorder(dt, cluster, -size, -len)
  dt[, rank := seq_len(.N), by = cluster]
  dt[cw, weight := i.weight, on = "cluster"]
  dt[, count := as.integer(rank <= weight)]
  setorder(dt, ord)
  dt[["count"]]
}

#' @title cells_isoforms_size_filter
#'
#' @description Apply cluster size filtering to each isoform in each cell.
#' @param cell_isoform_table A dataframe for reads including its cell source,
#' isoform sequence, splicing sites and umi cluster size.
#' @param ratio The estimated frequency for UMIs which are enriched of sequencing errors and
#' couldn't be clustered correctly and become singletons
#' @return A dataframe with isoform count after cluster size filtering for each cell
cells_isoforms_size_filter <- function(cell_isoform_table, ratio = 0.1, ...,
                                          thresh = 10, eps = 10) {
  dt <- as.data.table(cell_isoform_table)
  dt[, weight := isoform_size_filter(isoform, size, ratio, ...,
                                        thresh = thresh, eps = eps),
     by = .(cell, mid)]
  dt[, {
    sw <- sum(weight)
    .(size    = sum(size),
      cluster = .N,
      count   = sw,
      polyA   = if (sw > 0) sum(polyA * weight) / sw else NA_real_)
  }, by = .(cell, isoform)]
}


