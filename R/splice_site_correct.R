#' @title splice_site_table
#'
#' @description extract splice sites from isoform sequence.
#' @details extract splice sites from isoform sequence and filter out artifacts.
#'
#' @param isoform The string which represents the isoform
#' @param split The character to seperate different exons in the isoform representation.
#' @param sep The character to seperate the start and end position of an exon in the isoform representation.
#' @param splice_site_thresh The minimum number of occurance for a splice site to be preserved.
#' @return A dataframe, the first and the last column stores the start and end position of each read,
#' the middle columns store if each read has the splicing sites
splice_site_table <- function(isoforms,
                              #polyA,strand,
                              split = "|",sep = ",",
                              splice_site_thresh = 10){
  out = splice_site_table_cpp(isoforms,split,sep,splice_site_thresh)

  out$start = as.numeric(out$start)
  out$end = as.numeric(out$end)

  out = as.data.frame(do.call(cbind,out))
  if(ncol(out) > 3){
    end_thresh= 5
    mid_sites = as.numeric(colnames(out)[3:(ncol(out)-1)])
    out = lapply(1:nrow(out),function(i){
      x = unlist(out[i,])
      id = as.numeric(x["id"])
      start = as.numeric(x["start"])
      end = as.numeric(x["end"])
      x[as.character(mid_sites[(mid_sites > (end-end_thresh)) | (mid_sites < (start+end_thresh))])] = NA

      #if(polyA[id] == 0){
      #  x[as.character(mid_sites[(mid_sites > (end+end_thresh)) | (mid_sites < (start-end_thresh))])] = NA
      #}
      #else{
      #  if(strand == "+"){
      #    x[as.character(mid_sites[mid_sites < (start-end_thresh)])] = NA
      #  }
      #  else if(strand == "-"){
      #    x[as.character(mid_sites[mid_sites > (end + end_thresh)])] = NA
      #  }
      #}
      return(x)
    })
    out = as.data.frame(do.call(rbind,out))
    mid_sites = as.character(sort(mid_sites))
    out = as.data.frame(out[,c("id","start",mid_sites,"end")])
  }
  return(out)
}

#' @title mid_len
#'
#' @description count the number of middle splicing sites covered by a read
#'
#' @param mid A vector of strings, each string denotes the existence of the middle splicing sites
#' @param sep The character to seperate splicing sites in the string
#' @return A dataframe with two columns, the first column denotes the middle splicing sites,
#' the second denotes the number of sites which covered.
mid_len = function(mid,sep = ","){
  mat = strsplit(mid,split = sep,fixed = TRUE)
  mat = do.call(rbind,mat)
  suppressWarnings(storage.mode(mat) <- "numeric")

  out = as.data.frame(cbind(mid,rowSums(!is.na(mat))))
  colnames(out) = c("mid","size")
  out = out %>% mutate(size = as.numeric(size))
  return(out)
}

#' @title mid_group
#'
#' @description group middle splicing sites which can transform to each other by truncations.
#'
#' @inheritParams mid_len
#' @return A dataframe with two columns, the first column denotes the middle splicing site string,
#' the second denotes the splicing site string that can transform to this splicing site string by
#' truncations.
mid_group = function(mid,sep = ","){
  mat = strsplit(mid,split = sep,fixed = TRUE)
  mat = do.call(rbind,mat)
  suppressWarnings(storage.mode(mat) <- "numeric")


  mid = mid[order(rowSums(is.na(mat)),decreasing = TRUE)]
  mat = mat[order(rowSums(is.na(mat)),decreasing = TRUE),,drop = FALSE]


  if(nrow(mat) == 1){
    result = as.data.frame(cbind(mid,mid))
    colnames(result) = c("c","p")
    return(result)
  }

  is_subset = function(child, parent){
    if(length(child) != length(parent)){
      stop("The length of two input vectors should be the same!")
    }
    if(sum(is.na(child)|is.na(parent))!=sum(is.na(child))){
      return(FALSE)
    }
    diff = xor(child,parent)
    diff[is.na(diff)] = 0
    if(sum(diff) == 0){
      return(TRUE)
    }
    else{
      return(FALSE)
    }
  }

  result = lapply(1:(nrow(mat)-1),function(i){
    sub = lapply((i+1):nrow(mat),function(j){
      if(is_subset(mat[i,],mat[j,])){
        return(c(i,j))
      }
    })
    sub = do.call(rbind,sub)
  })
  result = as.data.frame(do.call(rbind,result))
  if(length(result) > 0){
    colnames(result) = c("c","p")
    result = result %>% mutate(c = mid[c],p = mid[p])
    orphan = setdiff(mid,result$c)
    if(length(orphan) > 0){
      orphan = as.data.frame(cbind(orphan,NA))
      colnames(orphan) = c("c","p")
      result = rbind(result,orphan)
    }
  }
  else{
    result = as.data.frame(cbind(mid,NA))
    colnames(result) = c("c","p")
  }

  return(result)
}

#' @title mid_count
#'
#' @description Count the frequency of middle splicing sites for each UMI cluster.
#' @details Count the frequency of middle splicing sites for each UMI cluster, the highest one will be chosen
#' as the representative and other splicing sites will be attributed to the representative.
#'
#' @param mid A string vector of middle splicing sites within a UMI cluster
#' @param total All types of middl splicing sites vector
#' @param parent The middle splicing site group output from "mid_group"
#' @param len The number of covered splicing sites for each read output from mid_len.
#' @return A dataframe, the first column is the original middle splicing sites, the second column is the chosen
#' representative within this UMI cluster, the third column stores the count of the original.
mid_count = function (mid, total, parent, len) {
  count_list = table(mid)
  count_mat = as.data.frame(count_list)

  colnames(count_mat) = c("mid","count")
  count_mat = count_mat %>% mutate(mid = as.character(mid))
  if(nrow(count_mat) < 2){
    mid_coexist = cbind(count_mat$mid,
                        count_mat$mid,
                        count_mat$count)
    return(mid_coexist)
  }
  count_mat = left_join(count_mat,len,by = "mid")
  count_mat = count_mat %>% arrange(-count,-size)
  count_mat = count_mat[1:2,]

  #start = Sys.time()
  count_mat = suppressWarnings(left_join(count_mat,
                                         parent %>% filter(p %in% count_mat$mid),
                                         by = c("mid" = "c")))

  if(sum(is.na(count_mat$p)) == 1){
    mode_mid = count_mat$mid[is.na(count_mat$p)]
    mid_coexist = cbind(mode_mid, mode_mid, sum(count_mat$count))
  }
  else{
    mode_mid = count_mat$mid[1]
    mid_coexist = cbind(count_mat$mid, mode_mid, count_mat$count)
  }
  rownames(mid_coexist) = NULL
  #end = Sys.time()
  #print(end - start)
  return(mid_coexist)
}

#' @title disagree_sites
#'
#' @description Identify the different sites between different isoforms within the same
#' UMI cluster
#'
#' @param from,to A string vector of middle splicing sites. From is the original middle splicing
#' sites, to is the middle splicing sites the from mapped to within the same UMI cluster
#' @return A dataframe, recording the number of splicing sites existing in the correctly mapped
#' and wrongly mapped reads.
disagree_sites = function(from,to){
  if(length(from) != length(to)){
    stop("There should be a one to one correspondence!")
  }
  from_table = do.call(rbind,strsplit(from,split = ","))
  to_table = do.call(rbind,strsplit(to,split = ","))

  from_table <- suppressWarnings(matrix(as.numeric(from_table),ncol = ncol(from_table)))
  to_table <- suppressWarnings(matrix(as.numeric(to_table),ncol = ncol(to_table)))
  #return(list(from_table,to_table))
  disagree = xor(from_table,to_table) & (from_table == 1)

  filter = as.data.frame(cbind(colSums(disagree,na.rm = TRUE),
                               colSums(from_table,na.rm = TRUE),
                               colSums(to_table,na.rm = TRUE)))
  colnames(filter) = c("disagree","wrong","correct")

  return(filter)
}

#' @title site_correct
#'
#' @description correct wrongly mapped splicing sites.
#'
#' @inheritParams disagree_sites
#' @param sep The character to seperate splicing sites in the middle splicing site string.
#' @return A datafrom recording the mapping between wrongly mapped isoform and correct isoform.
site_correct = function(from,to,sep = ","){
  if(length(from) != length(to)){
    stop("The length of from and to should match!")
  }
  out = sapply(1:length(from),function(i){
    if(is.na(to[i])){
      return(to[i])
    }
    if(from[i] == to[i]){
      return(to[i])
    }
    else{
      from_vec = suppressWarnings(as.numeric(unlist(strsplit(from[i],split = sep,fixed = TRUE))))
      to_vec = suppressWarnings(as.numeric(unlist(strsplit(to[i],split = sep,fixed = TRUE))))

      from_vec[!is.na(from_vec) & !is.na(to_vec)] = to_vec[!is.na(from_vec) & !is.na(to_vec)]
      return(paste(from_vec,collapse = sep))
    }
  })
  return(out)

}

#' @title cells_mid_filter
#'
#' @description Filter out middle splice sites pattern caused by wrong mapping.
#' @details Filter out middle splice sites pattern caused by wrong mapping given the information from UMI cluster as
#' reads within the same UMI cluster should come from the same transcript.
#'
#' @param cells The string vector indicating the cell
#' @param cluster The factor vector indicating the UMI cluster
#' @param isoform_mid The string vector indicating middle splicing sites.
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr reframe
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr select
#' @importFrom tidyr pivot_wider
#' @return A string vector indicating preserved patterns of middle splicing sites

mid_correct_input = function(cells,cluster,gene_isoform){
  if(nrow(gene_isoform) != length(cells)){
    stop("The size of isoforms and cells don't match!")
  }
  if(length(cluster) != length(cells)){
    stop("The size of clusters and cells don't match!")
  }
  mid = apply(gene_isoform[,2:(ncol(gene_isoform)-1)],
              1,function(x) paste(x,collapse = ","))
  out = as.data.frame(cbind(cells,cluster,
                            gene_isoform$start,mid,gene_isoform$end))
  colnames(out) = c("cell","cluster","start","mid","end")
  return(out)
}

mid_coexist = function(data){
  total = sort(table(data$mid),decreasing = TRUE)
  total = names(total)

  len = mid_len(total)
  parent = mid_group(total)

  data = data %>% group_by(cell,cluster) %>%
    reframe(coexist = mid_count(mid,total,parent,len),.groups = "drop")
  coexist = as.data.frame(data$coexist)
  colnames(coexist) = c("from","to","count")
  coexist$count = as.numeric(coexist$count)
  concensus = as.data.frame(cbind(data[,c("cell","cluster")],
                                  coexist[,"to",drop = FALSE])) %>%
                                  group_by(cell,cluster) %>%
                                  summarise(concensus = unique(to),.groups = "drop")

  coexist = coexist %>% group_by(from,to) %>%
    summarise(count = sum(count),.groups = "drop")
  return(list(concensus,coexist))
}

isoform_corres = function(coexist_matrix){
  isoform_coexist_filter <- cbind(diag(coexist_matrix),
                                  rowSums(coexist_matrix)-diag(coexist_matrix))
  isoform_coexist_filter <- isoform_coexist_filter[order(isoform_coexist_filter[,1],
                                                         decreasing = T),]
  isoform_coexist_filter <- names(which(isoform_coexist_filter[,1] >
                                          2*isoform_coexist_filter[,2]))

  corres = as.data.frame(coexist_matrix[,isoform_coexist_filter])
  corres = sapply(1:nrow(corres),function(i){
    x = unlist(corres[i,])
    if(sum(x) == 0){
      return(NA)
    }
    return(isoform_coexist_filter[which(x == max(x))[1]])
  })

  corres = as.data.frame(cbind(rownames(coexist_matrix),corres))
  colnames(corres) = c("from","to")
  rownames(corres) = corres$from

  disagree = corres %>% filter(from != to,!is.na(to))
  if(nrow(disagree) == 0){
    return(corres)
  }
  ds = disagree_sites(disagree$from,disagree$to)
  id = which(ds$wrong > ds$correct & ds$disagree*3 > ds$wrong)

  if(length(id) == 0){
    return(corres)
  }

  correct_iso = na.omit(unique(corres$to))
  to_table = do.call(rbind,strsplit(correct_iso,split = ","))
  to_table <- suppressWarnings(matrix(as.numeric(to_table),ncol = ncol(to_table)))
  to_table[is.na(to_table)] = 0
  if(length(id) == 1){
    correct_iso = correct_iso[to_table[,id] == 0]
  }
  else if(nrow(to_table) == 1){
    correct_iso = correct_iso[sum(to_table[,id]) == 0]
  }
  else{
    correct_iso = correct_iso[rowSums(to_table[,id]) == 0]
  }

  corres_new = as.data.frame(coexist_matrix[,correct_iso])
  corres_new = sapply(1:nrow(corres),function(i){
    x = unlist(corres_new[i,])
    if(sum(x) == 0){
      return(NA)
    }
    return(correct_iso[which(x == max(x))[1]])
  })

  corres_new = as.data.frame(cbind(rownames(coexist_matrix),corres_new))
  colnames(corres_new) = c("from","to")
  rownames(corres_new) = corres_new$from
  return(corres_new)
}

#' @title cluster_isoform_correct
#'
#' @description Corret for wrong mapping and truncation within each UMI cluster
#' @details Corret for wrong mapping and truncation within each UMI cluster
#'
#' @param start A numeric vector indicating start position of reads within the UMI cluster
#' @param mid A string vector indicating patterns of middle splicing sites
#' @param end A numeric vector indicating end position of reads within the UMI cluster
#' @param polyA A boolean vector indicating the polyA existence
#' @param preserve_mid A string vector indicating preserved patterns of middle splice sites
#' @return A vector including the corrected start, middle splicing sites, end position, polyA existence
#' and the number of reads within this UMI cluster
cluster_isoform_correct <- function(start,mid,end,concensus,polyA,preserve_mid){
    if(length(start) != length(mid) || length(start) != length(end)){
      stop("The size of isoforms representation don't match!")
    }
    if(length(start) != length(polyA)){
        stop("The size of isoforms and polyA don't match!")
    }
    cluster_size = length(start)
    mode_isoform = NA
    mode_start = NA
    mode_end = NA
    mode_polyA = NA

    if(length(preserve_mid) > 0 && nrow(preserve_mid) > 0){
      mid_corres = na.omit(preserve_mid[concensus,"to"])
      if(length(mid_corres) > 0){
        mode_isoform = mid_corres
        mode_preserve = which(mid == mode_isoform)
        if(length(mode_preserve) == 0){
          mode_start = -1
          mode_end = -1
        }
        else{
          mode_start = min(start[mode_preserve])
          #mode_start = table(start[mode_preserve])
          #mode_start = min(as.numeric(names(mode_start[which(mode_start == max(mode_start))])))

          mode_end = max(end[mode_preserve])
          #mode_end = table(end[start == mode_start])
          #mode_end = max(as.numeric(names(mode_end[which(mode_end == max(mode_end))])))
        }
        polyA = as.logical(polyA)
        mode_polyA = mean(as.numeric(polyA[mid == concensus]))
      }
    }

    return(c(mode_start,mode_isoform,mode_end,cluster_size,mode_polyA))
}

#' @title isoform_correct
#'
#' @description Corret for wrong mapping and truncation for each UMI cluster in each gene in each cell
#' @details Corret for wrong mapping and truncation for each UMI cluster in each gene in each cell
#'
#' @param gene_isoform A data frame containing cell, gene,and isoform informaion including
#' start, middle splicing sites, end position and polyA existence
#' @param preserve_mid A string vector indicating preserved patterns of middle splice sites
#' @importFrom tidyr replace_na
#' @importFrom dplyr mutate_at
#' @importFrom dplyr group_by
#' @return A dataframe including the corrected start, middle splicing sites, end position, polyA existence
#' and the number of reads of each UMI cluster for each gene in each cell
isoform_correct <- function(gene_isoform,preserve_mid){
  gene_isoform_adjust = gene_isoform %>% group_by(cell,cluster) %>%
                        summarise(adjust = list(cluster_isoform_correct(start,mid,end,unique(concensus),
                                                                        polyA,preserve_mid)),.groups = "drop")
  adjust = as.data.frame(do.call(rbind,gene_isoform_adjust$adjust))
  colnames(adjust) = c("start","mid","end","size","polyA")

  adjust = adjust %>% mutate_at(c("start","end","size","polyA"),as.numeric)
  #adjust = adjust %>% group_by(mid) %>%
  #  mutate_at(c("start","end"),~tidyr::replace_na(., median(., na.rm=TRUE)))

  out = as.data.frame(cbind(gene_isoform_adjust[,c("cell")],adjust))
  return(out[!is.na(out$start),])
}

cells_mid_correct <- function(cells,cluster,gene_isoform,polyA){
  data = mid_correct_input(cells,cluster,gene_isoform)

  out = mid_coexist(data)
  concensus = out[[1]]
  coexist = out[[2]]

  if(nrow(coexist) == 1){
    corres = as.data.frame(coexist %>% dplyr::select(-count))
    rownames(corres) = corres$from
  }
  # return(coexist)
  else{
    coexist_matrix = long2square(coexist,"from","to","count",symmetric = FALSE)
    corres = isoform_corres(coexist_matrix)
  }

  data = left_join(data,concensus,by = c("cell","cluster"))
  data$polyA = polyA
  # return(data)
  data = isoform_correct(data,corres)
  return(data)
}

cells_nomid_correct <- function(cells,cluster,gene_isoform,polyA){
  data = as.data.frame(cbind(cells,cluster,gene_isoform,polyA))
  colnames(data) = c("cell","cluster","start","end","polyA")

  data = data %>% group_by(cell,cluster) %>%
    summarise(start = min(start),end = max(end),
              size = n(),
              polyA = mean(polyA),.groups = "drop")
  return(data)
}

#' @title site_recover
#'
#' @description transform the middle splicing sites pattern from boolean string into original sites
#' @details transform the middle splicing sites pattern from boolean string into original sites
#'
#' @param start A numeric value indicating start position of reads within the UMI cluster
#' @param mid A string indicating patterns of middle splicing sites
#' @param end A numeric value indicating end position of reads within the UMI cluster
#' @param sites A string vector indicating middle splicing sites
#' @param sep The character to seperate the start and end for each exon bin in the isoform representation
#' @param split The character to seperate each exon bin in the isoform representation
#' @return A string representing isoform
site_recover <- function(start,end,mid = NULL,sites = NULL,flank = 5,sep = ",",split = "|"){
    if(is.null(sites)){
        splice_sites = paste(start,end,sep = sep)
    }
    else{
      mid = suppressWarnings(as.logical(as.numeric(unlist(strsplit(mid,split = ",")))))
      mid[is.na(mid)] = FALSE
      if(length(mid) != length(sites)){
        stop("The size of splicing sites and binary indicator don't match!")
      }
      mid_sites = sites[mid]
      if(as.numeric(start) == -1){
        start = as.numeric(mid_sites[1])-flank
      }
      if(as.numeric(end) == -1){
        end = as.numeric(mid_sites[length(mid_sites)])+flank
      }
      if(is.na(start)){
        return(NA)
      }
      splice_sites = c(start,sites[mid],end)
      splice_sites = t(matrix(splice_sites,nrow = 2))
      splice_sites = apply(splice_sites,1,function(x) paste(x,collapse = sep))
      splice_sites = paste(splice_sites,collapse = split)
    }
    splice_sites = gsub(" ","",splice_sites,fixed = TRUE)
    return(splice_sites)
}

cells_isoform_correct <- function(cells,cluster,gene_isoform,polyA){
  if(ncol(gene_isoform) > 2){
    splice_sites = colnames(gene_isoform)[2:(ncol(gene_isoform)-1)]
    data = cells_mid_correct(cells,cluster,gene_isoform,polyA)
    if(nrow(data) == 0){
      return(NULL)
    }
    data$isoform <- apply(data,1,function(x){
      site_recover(x["start"],x["end"],x["mid"],splice_sites)
    })
  }
  else{
    splice_sites = NULL
    data = cells_nomid_correct(cells,cluster,gene_isoform,polyA)
    if(nrow(data) == 0){
      return(NULL)
    }
    data$isoform <- apply(data,1,function(x){
      site_recover(x["start"],x["end"])
    })
    data$mid = "null"
  }

  data <- na.omit(data) %>% dplyr::select(-c(start,end))
  return(data)
}
