#' @title intron_only
#' @description judge if a read only contains intron sequence
#' @details This function compares the coverage between the read and the annotated isoform, if the read
#' has no overlap with the annotated isoform, this read is classified as intron.
#' @param reads A string vector contains strings for reads
#' @param gtf The gtf annotation, each row is an exon for an isoform
#' @param gtf_start_col the name of the column which stores the start position of exon in the gtf.
#' @param gtf_end_col the name of the column which stores the end position of exon in the gtf.
#' @param sep The character to split the start and end position for each exon in the isoform.
#' @param split The character to split the exons in the isoform
#' @importFrom valr bed_subtract
#' @return A boolean vector.
intron_only = function(reads,gtf,gtf_start_col = "start",gtf_end_col = "end",
                       sep = ",",split = "|"){
  exon_bin = as.data.frame(gtf[,c(gtf_start_col,gtf_end_col)])
  colnames(exon_bin) = c("start","end")
  exon_bin$chrom = "chr1"

  intron_flag = sapply(reads,function(x){
    bins = read2bins(x,sep,split)
    bins$chrom = "chr1"
    diff = valr::bed_subtract(bins,exon_bin)

    if(binsum(diff) == binsum(bins)){
      return(TRUE)
    }
    else{
      return(FALSE)
    }
  })

  return(intron_flag)
}


#' @title iso_full_dis
#' @description calculate the difference in the alignment between two reads in base pair
#' @param read1,read2 A string or matrix to record the exons in a read
#' @param sep The character to split the start and end position for each exon in the isoform.
#' @param split The character to split the exons in the isoform
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomicRanges width
#' @importFrom GenomicRanges pintersect
#' @importFrom GenomicRanges queryHits
#' @importFrom GenomicRanges subjectHits
#' @return A numerical number to denote the difference.
iso_full_dis = function(read1,read2,split = "|",sep = ","){
  read1 = read2bins(read1,split = "|",sep = ",")
  read2 = read2bins(read2,split = "|",sep = ",")

  bins1 <- GRanges(seqnames = "chr1", ranges = IRanges(start = read1$start, end = read1$end))
  bins2 <- GRanges(seqnames = "chr1", ranges = IRanges(start = read2$start, end = read2$end))

  overlaps <- findOverlaps(bins1, bins2)
  overlap_lengths <- width(pintersect(bins1[queryHits(overlaps)], bins2[subjectHits(overlaps)]))
  total_difference <- sum(width(bins1)) + sum(width(bins2)) - 2 * sum(overlap_lengths)

  return(total_difference)
}

#' @title iso_full_dis_v
#' @description The vectorized version of the iso_full_dis function
#' @inheritParams iso_full_dis
iso_full_dis_v <- Vectorize(iso_full_dis, c("read1", "read2"))


#' @title iso_corres
#' @description map the reads to the annotated isoforms for one gene
#' @details map the reads to the annotated isoforms for one gene
#' @inheritParams intron_only
#' @param transcripts A string vector storing reads
#' @param gene The gene name, should exist in the gtf annotation
#' @param thresh The maximum threshold for the total offset of middle splicing sites.
#' @param overlap_thresh The minimum threshold for the coverage of the annotated isoform.
#' @param end_bias The maximum threshold for the offset of start and end position.
#' @param gtf_gene_col the name of the column which stores gene name in the gtf.
#' @param gtf_iso_col the name of the column which stores isoform name in the gtf.
#' @importFrom dplyr reframe
#' @importFrom dplyr arrange_at
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise_all
#' @importFrom rlang sym
#' @return A dataframe with two columns, the first column records the reads, the second column records the isoform
#' the read could align to.
#' @export
iso_corres = function(transcripts,gene,gtf,thresh = 3,overlap_thresh = 0.25,
                             end_bias = 200,
                             gtf_gene_col = "gene",gtf_iso_col = "transname",
                             gtf_start_col = "start",gtf_end_col = "end",
                             sep = ",",split = "|"){
  sub_gtf = gtf %>% filter_at(gtf_gene_col,~.==gene) %>%
    arrange_at(c(gtf_iso_col,gtf_start_col,gtf_end_col))

  sub_gtf_iso = sub_gtf %>% group_by_at(gtf_iso_col) %>%
    reframe(iso = paste(paste(!!sym(gtf_start_col),
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
  # transcripts_iso_corres = transcripts_iso_corres[,c("isoform","transname","overlap","dis")]
  # return(transcripts_iso_corres)
  suppressWarnings({
    transcripts_iso_corres = transcripts_iso_corres %>% group_by(isoform) %>%
      filter(dis == min(dis)) %>%
      dplyr::select(isoform,transname,overlap) %>%
      arrange(isoform,-overlap)
  })
  transcripts_iso_corres = as.data.frame(transcripts_iso_corres)
  transcripts_iso_corres = transcripts_iso_corres %>%
    mutate(diff = iso_full_dis_v(transcripts_uniq[isoform + 1],sub_gtf_iso[transname+1,"iso"]),
           isoform = transcripts_uniq[isoform+1],
           transname = sub_gtf_iso[transname+1,gtf_iso_col]) %>%
    arrange(isoform,diff)

  transcripts_iso_corres = transcripts_iso_corres[!duplicated(transcripts_iso_corres$isoform),]
  #return(transcripts_iso_corres)
  #transcripts_iso_corres$isoform = transcripts_uniq[transcripts_iso_corres$isoform+1]
  #transcripts_iso_corres$transname = sub_gtf_iso[transcripts_iso_corres$transname+1,gtf_iso_col]


  return(transcripts_iso_corres)
}

#' @title cell_iso_count_impute
#' @description imputate the isoform count for one gene in a cell.
#' @details Some truncated reads can have ambigous alignment to multiple annotated isoforms, this function will
#' iterative imputate the isoform count based on non-ambigous alignment.
#' @param data The input dataframe, including the reads, reads count and its aligned annotated isoform.
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

#' @title iso_count_impute
#' @description imputate the isoform count for one gene in multiple cells.
#' @details Some truncated reads can have ambigous alignment to multiple annotated isoforms, this function will
#' iterative imputate the isoform count based on non-ambigous alignment.
#' @inheritParams cell_iso_count_impute
#' @param cell_col The name of the column in the input data which records the cell barcode.
#' @param iso_col The name of the column in the input data which records the isoform.
#' @param overlap_col The name of the column in the input data which records the ratio of coverage.
#' @param count_col The name of the column in the input data which records the reads count.
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom dplyr group_by_at
#' @importFrom dplyr summarise_at
#' @importFrom dplyr ungroup
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

#' @title cells_genes_isos_count
#' @description imputate the isoform count for multiple genes in multiple cells.
#' @details Some truncated reads can have ambigous alignment to multiple annotated isoforms, this function will
#' iterative imputate the isoform count based on non-ambigous alignment.
#' @inheritParams iso_corres
#' @inheritParams iso_count_impute
#' @param gene_col The name of the column in the input data which records the gene name.
#' @param transcript_col The name of the column in the input data which records the isoform.
#' @param filter_only_intron A boolean to indicate if read only cover intron part should be preserved.
#' @importFrom magrittr %>%
#' @importFrom dplyr left_join
#' @importFrom dplyr group_by_at
#' @importFrom dplyr filter_at
#' @importFrom dplyr summarise
#' @importFrom tidyr pivot_wider
#' @export
cells_genes_isos_count = function(data,gtf,thresh = 3,overlap_thresh = 0.25,
                                  filter_only_intron = TRUE,
                                  cell_col = "cell",gene_col = "gene",
                            transcript_col = "isoform",count_col = "count",
                            gtf_gene_col = "gene",gtf_iso_col = "transname",
                            gtf_start_col = "start",gtf_end_col = "end",
                            split = "|",sep = ","){
  data = as.data.frame(data)
  gene_uniq = unique(data[,gene_col])

  out = lapply(gene_uniq,function(i){
    print(i)
    sub_data = data %>% filter_at(gene_col,~.==i)
    sub_gtf = gtf %>% filter_at(gtf_gene_col,~.==i) %>%
      arrange_at(c(gtf_iso_col,gtf_start_col,gtf_end_col))

    if(filter_only_intron){
      transcripts_uniq = unique(sub_data[,transcript_col])
      intron_flag = intron_only(transcripts_uniq,sub_gtf,
                                gtf_start_col,gtf_end_col,
                                sep,split)
      transcripts_uniq = transcripts_uniq[!intron_flag]
      if(length(transcripts_uniq) == 0){
        return(NULL)
      }
      sub_data = sub_data %>% filter_at(transcript_col,~.%in% transcripts_uniq)
    }

    iso_index = iso_corres(sub_data[,transcript_col],
                           gene = i,gtf = sub_gtf,thresh = thresh,
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
      #iso_index = iso_index %>% group_by_at(transcript_col) %>%
      #  summarise(transname = list(transname),overlap = list(overlap))
      sub_data = left_join(sub_data,iso_index,by = transcript_col)
      sub_data$transname[sapply(sub_data$transname,function(i) is.null(i))] = "unknown"
    }

    sub_out = sub_data %>% group_by(cell,transname) %>% summarise(count = sum(count),.groups = "drop")
    #sub_out = iso_count_impute(sub_data,cell_col = cell_col,
    #                           iso_col = "transname",overlap_col = "overlap",
    #                           count_col = count_col)
    sub_out$gene = i
    return(sub_out)
  })
  out = do.call(rbind,out)
  #out = as.data.frame(tidyr::pivot_wider(out,names_from = "cell",values_from = "count"))

  #out[is.na(out)] = 0
  return(out)
}
