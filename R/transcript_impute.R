#' @title isos2exonids_per_gene
#' @description map the exon to the annotated exons gene bed for one gene
#' @details transform the exon representation from start end position to exon id in gene bed
#' @param isoform A string vector storing isoforms
#' @param bed gene bed annotation for the specific gene
#' @param start the name of the column which stores the start position for each exon in the gene bed.
#' @param end the name of the column which stores the end position for each exon in the gene bed.
#' @param id the name of the column which stores the exon id in the gene bed.
#' @param mid_bias The tolerance for the bias for mapping of exons in the middle of the isoform
#' @param end_bias The tolerance for the bias for mapping of exons at the end of the isoform
#' @param end_overlap The minimum length of the end exon
#' @param nonsense_label The label for exons which couldn't be mapped to gene bed
#' @param split The character to split the exons in the isoform
#' @param sep The character to split the start and end position for each exon in the isoform.
#' @return A named vector storing isoforms in exon id form.
iso_corres = function(transcripts,gene,gtf,thresh = 3,overlap_thresh = 0.5,
                      end_bias = 200,
                      gtf_gene_col = "gene",gtf_iso_col = "transname",
                      gtf_start_col = "start",gtf_end_col = "end",
                      sep = ",",split = "|"){
  sub_gtf = gtf %>% filter_at(gtf_gene_col,~.==gene) %>%
    arrange_at(c(gtf_iso_col,gtf_start_col,gtf_end_col))
  sub_gtf_iso = sub_gtf %>% group_by_at(gtf_iso_col) %>%
    summarise(iso = paste(paste(!!sym(gtf_start_col),
                                !!sym(gtf_end_col),sep = sep),collapse = split))
  sub_gtf_iso = as.data.frame(sub_gtf_iso)

  transcripts_uniq = unique(transcripts)

  #return(list(transcripts_uniq,sub_gtf_iso))
  transcripts_iso_corres = isoset_mid_diff(transcripts_uniq,sub_gtf_iso$iso,
                                           thresh,overlap_thresh,end_bias,split,sep)
  if(transcripts_iso_corres[1,1] == -1){
    return(NULL)
  }

  transcripts_iso_corres = as.data.frame(transcripts_iso_corres)
  colnames(transcripts_iso_corres) = c("isoform","transname","dis","overlap")
  transcripts_iso_corres = transcripts_iso_corres[,c("isoform","transname","overlap","dis")]
  transcripts_iso_corres = transcripts_iso_corres %>% group_by(isoform) %>%
                           summarise_all(~.[dis == min(dis)]) %>%
                           dplyr::select(isoform,transname,overlap) %>%
                           arrange(isoform,-overlap)
  transcripts_iso_corres = as.data.frame(transcripts_iso_corres)
  #return(transcripts_iso_corres)
  transcripts_iso_corres$isoform = transcripts_uniq[transcripts_iso_corres$isoform+1]
  transcripts_iso_corres$transname = sub_gtf_iso[transcripts_iso_corres$transname+1,gtf_iso_col]

  return(transcripts_iso_corres)
}

cell_iso_count_impute = function(data){
  colnames(data) = c("isoform","overlap","level","count")

  level_uniq = sort(unique(data$level))

  isoform_uniq = unique(unlist(data$isoform))
  isoform_count = rep(0,length(isoform_uniq))
  names(isoform_count) = isoform_uniq

  for(i in level_uniq){
    sub_data = data %>% filter(level == i)
    level_freq = lapply(1:nrow(sub_data),function(j){
      iso = unlist(sub_data$isoform[j])
      overlap = unlist(sub_data$overlap[j])

      freq = rep(0,length(isoform_uniq))
      names(freq) = isoform_uniq
      if(sum(isoform_count[iso]) == 0){
        iso_filter = iso[overlap == max(overlap)]
        freq[iso_filter] = 1/length(iso_filter)
      }
      else{
        freq[iso] = isoform_count[iso]/sum(isoform_count[iso])
      }
      return(freq)
    })
    level_freq = do.call(rbind,level_freq)
    level_freq = colSums(level_freq*sub_data$count)
    isoform_count = isoform_count+level_freq
  }
  return(isoform_count)
}

iso_count_impute = function(data,cell_col = "cell",
                            iso_col = "isoform",overlap_col = "overlap",
                            count_col = "count"){
  data$level = sapply(data[,iso_col],length)
  data[data$level == 1,"overlap"] = list(1)

  data = data %>% group_by_at(c(cell_col,iso_col,overlap_col,"level")) %>%
    summarise_at(count_col,~sum(.))
  colnames(data) = c("cell","isoform","overlap","level","count")
  data = ungroup(data)

  cell_uniq = unique(data$cell)
  count_impute = lapply(cell_uniq,function(x){
    sub_data = data %>% filter(cell == x)
    sub_count_impute = cell_iso_count_impute(sub_data %>% dplyr::select(-cell))
    out = as.data.frame(cbind(x,names(sub_count_impute),sub_count_impute))
    colnames(out) = c("cell","isoform","count")
    return(out)
  })
  count_impute = do.call(rbind,count_impute)
  count_impute$count = as.numeric(count_impute$count)
  return(count_impute)
}

cells_genes_isos_count = function(data,gtf,thresh = 3,overlap_thresh = 0.5,
                                  cell_col = "cell",gene_col = "gene",
                            transcript_col = "isoform",count_col = "count",
                            gtf_gene_col = "gene",gtf_iso_col = "transname",
                            gtf_start_col = "start",gtf_end_col = "end",
                            split = "|",sep = ","){
  data = as.data.frame(data)
  gene_uniq = unique(data[,gene_col])

  out = lapply(gene_uniq,function(i){
    #print(i)
    sub_data = data %>% filter_at(gene_col,~.==i)
    iso_index = iso_corres(sub_data[,transcript_col],
                           gene = i,gtf = gtf,thresh = thresh,
                           overlap_thresh = overlap_thresh,
                           gtf_gene_col = gtf_gene_col,
                           gtf_iso_col = gtf_iso_col,
                           gtf_start_col = gtf_start_col,
                           gtf_end_col = gtf_end_col,
                           sep = sep,split = split)
    if(is.null(iso_index)){
      sub_data$transname = "unknown"
    }
    else{
      iso_index = iso_index %>% group_by_at(transcript_col) %>%
        summarise(transname = list(transname),overlap = list(overlap))
      sub_data = left_join(sub_data,iso_index,by = transcript_col)
      sub_data$transname[sapply(sub_data$transname,function(i) is.null(i))] = "unknown"
    }

    sub_out = iso_count_impute(sub_data,cell_col = cell_col,
                               iso_col = "transname",overlap_col = "overlap",
                               count_col = count_col)
    sub_out$gene = i
    return(sub_out)
  })
  out = do.call(rbind,out)
  out = as.data.frame(tidyr::pivot_wider(out,names_from = "cell",values_from = "count"))

  out[is.na(out)] = 0
  return(out)
}
