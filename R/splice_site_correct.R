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
splice_site_table <- function(isoforms,polyA,strand,
                              split = "|",sep = ",",
                              splice_site_thresh = 10){
  out = splice_site_table_cpp(isoforms,split,sep,splice_site_thresh)

  out$start = as.numeric(out$start)
  out$end = as.numeric(out$end)

  out = as.data.frame(do.call(cbind,out))
  if(ncol(out) > 3){
    mid_sites = as.numeric(colnames(out)[3:(ncol(out)-1)])
    out = lapply(1:nrow(out),function(i){
      x = unlist(out[i,])
      id = as.numeric(x["id"])
      start = as.numeric(x["start"])
      end = as.numeric(x["end"])
      if(polyA[id] == 0){
        x[as.character(mid_sites[mid_sites > end | mid_sites < start])] = NA
      }
      else{
        if(strand == "+"){
          x[as.character(mid_sites[mid_sites < start])] = NA
        }
        else if(strand == "-"){
          x[as.character(mid_sites[mid_sites > end])] = NA
        }
      }
      return(x)
    })
    out = as.data.frame(do.call(rbind,out))
    mid_sites = as.character(sort(mid_sites))
    out = as.data.frame(out[,c("id","start",mid_sites,"end")])
  }
  return(out)
}

#' @title mid_count
#'
#' @description Count the frequency of middle splicing sites for each UMI cluster.
#' @details Count the frequency of middle splicing sites for each UMI cluster, the highest one will be chosen
#' as the representative and other splicing sites will be attributed to the representative.
#'
#' @param mid A string vector of middle splicing sites within a UMI cluster
#' @return A dataframe, the first column is the original middle splicing sites, the second column is the chosen
#' representative within this UMI cluster, the third column stores the count of the original.
mid_count = function(mid,total){
  count = table(mid)
  name = names(count)
  #count = count*log(sum(count),2)

  if(length(name) == 1){
    mid_coexist = cbind(name,name,count)
  }
  else{
    mode_mid = name[which(count == max(count))]
    mode_mid = total[total %in% mode_mid][1]
    mid_coexist = cbind(name,mode_mid,count)
  }
  rownames(mid_coexist) = NULL
  return(mid_coexist)
}


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
cells_mid_filter <- function(cells,cluster,isoform_mid){
  if(length(isoform_mid) != length(cells)){
    stop("The size of isoforms and cells don't match!")
  }
  if(length(cluster) != length(cells)){
    stop("The size of clusters and cells don't match!")
  }
  total = sort(table(isoform_mid),decreasing = TRUE)
  total = names(total)

  temp = as.data.frame(cbind(cells,cluster,isoform_mid))
  colnames(temp) = c("cell","cluster","mid")

  temp = temp %>% group_by(cell,cluster) %>% reframe(coexist = mid_count(mid,total))
  coexist = as.data.frame(temp$coexist)
  colnames(coexist) = c("from","to","count")
  coexist$count = as.numeric(coexist$count)
  coexist = coexist %>% group_by(from,to) %>% summarise(count = sum(count),.groups = "drop")
  if(nrow(coexist) == 1){
    isoform_corres = as.data.frame(coexist %>% dplyr::select(-count))
    rownames(isoform_corres) = isoform_corres$from
    return(isoform_corres)
  }
  mid_uniq = unique(c(coexist$from,coexist$to))
  coexist_matrix = as.data.frame(tidyr::pivot_wider(coexist,names_from = "to",values_from = "count"))
  rownames(coexist_matrix) = coexist_matrix$from
  coexist_matrix = coexist_matrix %>% dplyr::select(-from)
  coexist_matrix[,setdiff(mid_uniq,colnames(coexist_matrix))] = NA
  coexist_matrix = coexist_matrix[mid_uniq,mid_uniq]
  coexist_matrix[is.na(coexist_matrix)] = 0
  coexist_matrix = as.matrix(coexist_matrix)
  #return(coexist_matrix)

  isoform_coexist_filter <- cbind(diag(coexist_matrix),
                                  rowSums(coexist_matrix)-diag(coexist_matrix))
  isoform_coexist_filter <- isoform_coexist_filter[order(isoform_coexist_filter[,1],
                                                         decreasing = T),]
  isoform_coexist_filter <- names(which(isoform_coexist_filter[,1] > 2*isoform_coexist_filter[,2]))

  isoform_corres = as.data.frame(coexist_matrix[,isoform_coexist_filter])
  isoform_corres = sapply(1:nrow(isoform_corres),function(i){
    x = unlist(isoform_corres[i,])
    if(sum(x) == 0){
      return(NA)
    }
    return(isoform_coexist_filter[which(x == max(x))[1]])
  })

  isoform_corres = as.data.frame(cbind(rownames(coexist_matrix),isoform_corres))
  colnames(isoform_corres) = c("from","to")
  rownames(isoform_corres) = isoform_corres$from
  #return(isoform_corres)

  disagree = isoform_corres %>% filter(from != to,!is.na(to))
  if(nrow(disagree) == 0){
    return(isoform_corres)
  }
  ds = disagree_sites(disagree$from,disagree$to)
  id = which(ds$wrong > ds$correct & ds$disagree*3 > ds$wrong)

  if(length(id) == 0){
    return(isoform_corres)
  }

  correct_iso = na.omit(unique(isoform_corres$to))
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

  isoform_corres_new = as.data.frame(coexist_matrix[,correct_iso])
  isoform_corres_new = sapply(1:nrow(isoform_corres_new),function(i){
    x = unlist(isoform_corres_new[i,])
    if(sum(x) == 0){
      return(NA)
    }
    return(correct_iso[which(x == max(x))[1]])
  })

  isoform_corres_new = as.data.frame(cbind(rownames(coexist_matrix),isoform_corres_new))
  colnames(isoform_corres_new) = c("from","to")
  #print(nrow(isoform_corres_new))
  #isoform_corres_new$to = site_correct(isoform_corres_new$from,isoform_corres_new$to)

  rownames(isoform_corres_new) = isoform_corres_new$from
  return(isoform_corres_new)
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
cluster_isoform_correct <- function(start,mid,end,polyA,preserve_mid){
    if(length(start) != length(mid) || length(start) != length(end)){
      stop("The size of isoforms representation don't match!")
    }
    if(length(start) != length(polyA)){
        stop("The size of isoforms and polyA don't match!")
    }
    cluster_size = length(start)

    if(length(preserve_mid) > 0 && nrow(preserve_mid) > 0 & !is.na(mid[1])){
        mid_corres = na.omit(preserve_mid[mid,"to"])

        if(length(mid_corres) == 0){
          mode_isoform = NA
          mode_start = NA
          mode_end = NA
          polyA = NA
          return(c(mode_start,mode_isoform,mode_end,cluster_size,polyA))
        }
        else{
          mode_isoform = table(mid_corres)
          mode_isoform = names(mode_isoform[which(mode_isoform == max(mode_isoform))])
          if(length(mode_isoform) > 1){
            loc= which(mid_corres %in% mode_isoform)
            mode_isoform = mode_isoform[which(loc == min(loc))]
          }
          mode_preserve = which(mid == mode_isoform)
        }
    }
    else{
        mode_preserve = 1:length(start)
        mode_isoform = NA
    }


    if(length(mode_preserve) == 0){
      mode_start = -1
      mode_end = -1
      polyA = 0
    }
    else{
      mode_start = table(start[mode_preserve])
      mode_start = min(as.numeric(names(mode_start[which(mode_start == max(mode_start))])))

      mode_end = table(end[start == mode_start])
      mode_end = max(as.numeric(names(mode_end[which(mode_end == max(mode_end))])))

      polyA = as.logical(polyA)
      polyA = mean(as.numeric(polyA[mode_preserve]))
    }

    return(c(mode_start,mode_isoform,mode_end,cluster_size,polyA))
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
                        summarise(adjust = list(cluster_isoform_correct(start,mid,end,polyA,preserve_mid)),.groups = "drop")
  adjust = as.data.frame(do.call(rbind,gene_isoform_adjust$adjust))
  colnames(adjust) = c("start","mid","end","size","polyA")

  adjust = adjust %>% mutate_at(c("start","end","size","polyA"),as.numeric)
  adjust = adjust %>% group_by(mid) %>%
    mutate_at(c("start","end"),~tidyr::replace_na(., median(., na.rm=TRUE)))

  out = as.data.frame(cbind(gene_isoform_adjust[,c("cell","cluster")],adjust))
  return(out[!is.na(out$start),])
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
site_recover <- function(start,mid,end,sites,sep = ",",split = "|"){
    if(is.na(mid)){
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
        start = as.numeric(mid_sites[1])-1
      }
      if(as.numeric(end) == -1){
        end = as.numeric(mid_sites[length(mid_sites)])+1
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
