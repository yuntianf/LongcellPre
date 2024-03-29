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
isos2exonids_per_gene = function(isoform,bed,
                        start = "start",end = "end",id = "id",
                        mid_bias = 0,end_bias = 10,
                        end_overlap =10,nonsense_label = "N",
                        split = "|",sep = ","){
  iso_uniq = unique(isoform)
  bed = as.data.frame(bed)
  bed[,id] = as.character(bed[,id])
  iso_uniq_index = isos2exonids_index(iso_uniq,bed[,start],bed[,end],bed[,id],
                                      mid_bias,end_bias,
                                      end_overlap,nonsense_label,
                                      split,sep)
  return(iso_uniq_index)
}

#' @title isos2exonids
#' @description map the exon to the annotated exons gene bed for all genes
#' @param data A dataframe storing reads information, the gene and isoform column are required
#' @param gene_col The name of the column which stores the gene name in the gene bed
#' @inheritParams isos2exonids_per_gene
#' @importFrom dplyr left_join
#' @return A dataframe as input data, but replace the isoforms by exon id sequences
isos2exonids = function(data,gene_bed,
                        gene_col = "gene",start = "start",end = "end",id = "id",
                        mid_bias = 0,end_bias = 10,
                        end_overlap =10,nonsense_label = "N",
                        split = "|",sep = ","){
  gene_isoform = unique(data %>% dplyr::select(gene,isoform))

  gene_uniq = unique(gene_isoform$gene)
  out = lapply(gene_uniq,function(i){
    sub_data = gene_isoform %>% filter(gene == i)
    sub_bed = gene_bed %>% filter(!!as.symbol(gene_col) == i)

    iso_uniq_index = isos2exonids_per_gene(sub_data$isoform,sub_bed,
                                           start = start,end = end,id = id,
                                           mid_bias = mid_bias,end_bias = end_bias,
                                           end_overlap =end_overlap,nonsense_label = nonsense_label,
                                           split = split,sep = sep)
    temp = as.data.frame(cbind(i,names(iso_uniq_index),unlist(iso_uniq_index)))
    rownames(temp) = NULL
    return(temp)
  })

  out = as.data.frame(do.call(rbind,out))
  colnames(out) = c("gene","isoform","exon_id")
  data = left_join(data,out,by = c("gene","isoform"))
  data = data %>% dplyr::select(-isoform)

  return(data)
}

#' @title bins2iso
#' @description Paste the exon bins (start and end) together as an isoform sequence
#' @param start A numerical vector storing the start positions for all exons in the isoform
#' @param end A numerical vector storing the end positions for all exons in the isoform
#' @param split The character to split the exons in the isoform
#' @param sep The character to split the start and end position for each exon in the isoform.
#' @return A string representing isoform
bins2iso= function(start,end,sep = ",",split = "|"){
  if(length(start) != length(end)){
    stop("The start and end should be paired, but the input size are different!")
  }
  bins = paste(start,end,sep = sep)
  iso = paste(bins,collapse = split)
  return(iso)
}

#' @title gtf_bed_corres
#' @description mapping canonical isoforms from gtf annotation to exon parts in gene bed
#' @param gtf A dataframe of gtf annotation, each row is an exon in an isoform
#' @param gene_bed A dataframe of gene_bed annotation, each row is an exon part
#' @param gtf_gene_col The name of the column storing the gene in gtf
#' @param gtf_start_col The name of the column storing the exon start position in gtf
#' @param gtf_end_col The name of the column storing the exon end position in gtf
#' @param transname The name of the column storing the transcript id in gtf
#' @param bed_gene_col The name of the column storing the gene in gene bed
#' @param bed_start_col The name of the column storing the exon start position in gene bed
#' @param bed_end_col The name of the column storing the exon end position in gene bed
#' @param id The name of the column storing the exon part id in gene bed
#' @return A data frame including 2 columns, one is the transcript name from gtf, the other is its exon part form.
gtf_bed_corres = function(gtf,gene_bed,genes = NULL,
                          gtf_gene_col = "gene",gtf_start_col = "start",gtf_end_col = "end",transname = "transname",
                          bed_gene_col = "gene",bed_start_col = "start",bed_end_col = "end",id = "id",
                          ...){
  gtf = as.data.frame(gtf)
  gene_bed = as.data.frame(gene_bed)
  all_genes = intersect(unique(gtf[,gtf_gene_col]),unique(gene_bed[,bed_gene_col]))

  if(is.null(genes)){
    genes = all_genes
  }
  else{
    if(length(setdiff(genes,all_genes)) > 0){
      warning("Some input genes don't exsit in the gene bed or gtf, those genes will be ignored!")
      genes = genes[genes %in% all_genes]
    }
  }

  gtf = gtf %>% filter(!!as.symbol(gtf_gene_col) %in% genes)
  gene_bed = gene_bed %>% filter(!!as.symbol(bed_gene_col) %in% genes)

  gtf_iso = gtf %>%
            group_by(!!as.symbol(gtf_gene_col),!!as.symbol(transname)) %>%
            summarise(isoform = bins2iso(!!as.symbol(gtf_start_col),
                                         !!as.symbol(gtf_end_col),
                                         sep = ",",split = "|"),.groups = "drop")
  colnames(gtf_iso)[colnames(gtf_iso) == gtf_gene_col] = "gene"

  gtf_iso = isos2exonids(gtf_iso,gene_bed,
                         gene_col = bed_gene_col,start = bed_start_col,end = bed_end_col,id = id,...)

  return(gtf_iso)
}


exon_ids2transname_per_gene = function(exon_id,gtf_exon_id,gene_i,
                                       gtf_gene_col = "gene",
                                       transname = "transname"){
  exon_id_uniq = unique(exon_id)
  gtf_exon_id = gtf_exon_id %>% filter(!!as.symbol(gtf_gene_col) == gene_i)
  corres = lapply(exon_id_uniq,function(x){
    out = gtf_exon_id$transname[grepl(x,gtf_exon_id$exon_id,fixed = TRUE)]
    if(length(out) == 0){
      out = NA
    }
    return(out)
  })
  names(corres) = exon_id_uniq
  out = corres[exon_id]
  names(out) = NULL
  return(out)
}

exon_ids2transname = function(data,gtf,gene_bed,
                              gtf_gene_col = "gene",gtf_start_col = "start",gtf_end_col = "end",transname = "transname",
                              bed_gene_col = "gene",bed_start_col = "start",bed_end_col = "end",id = "id",
                              ...){
  genes = unique(data$gene)
  gtf_exon_id = gtf_bed_corres(gtf,gene_bed,genes,
                               gtf_gene_col = gtf_gene_col,gtf_start_col = gtf_start_col,
                               gtf_end_col = gtf_end_col,transname = transname,
                               bed_gene_col = bed_gene_col,bed_start_col = bed_start_col,
                               bed_end_col = bed_end_col,id = id,
                               ...)
  data = data %>% group_by(gene) %>% mutate(transname = exon_ids2transname_per_gene(exon_id,gtf_exon_id,unique(gene)))

  return(data)
}


freq_local = function(freq,p){
  if(is.na(p[1])){
    return(1)
  }

  local = freq[p]
  local[is.na(local)] = 0

  total = sum(local)

  out = c(local/total)
  out[is.na(out)] = 0
  return(out)
}


transcript_count_impute = function(transcript,group,count){
  group_count = table(group)
  singleton = names(group_count)[group_count == 1]

  freq = table(transcript[group %in% singleton])
  freq = freq/sum(freq)

  data = as.data.frame(cbind(transcript,group))
  data = data %>% group_by(group) %>% mutate(prob = freq_local(freq,transcript))

  count = count * data$prob
  return(count)
}

transcript_impute = function(data,gtf,gene_bed,...){
  data = exon_ids2transname(data,gtf,gene_bed,...)
  data$id = 1:nrow(data)
  data = data %>% tidyr::unnest(transname)

  data = data %>% group_by(gene) %>% mutate(count_impute = transcript_count_impute(transname,id,count))
  data$transname[is.na(data$transname)] = "Unknown"

  out = data %>% dplyr::select(cell,gene,transname,count_impute)
  colnames(out)[3:4] = c("transcript","count")
  out = out %>% group_by(cell,gene,transcript) %>% summarise(count = sum(count))
  return(out)
}

save_cell_iso_count = function(data,mode = c("exon_id","transcript","both"),
                               gtf = NULL,gene_bed = NULL,path = "./",...){
  mode = match.arg(mode)

  cell_uniq = names(table(data$cell))
  gene_uniq = names(table(data$gene))

  write.table(cell_uniq,file = paste(path,"barcodes.tsv",sep = "/"),
              row.names = FALSE,col.names = FALSE,quote = FALSE)
  write.table(gene_uniq,file = paste(path,"features.tsv",sep = "/"),
              row.names = FALSE,col.names = FALSE,quote = FALSE)


  if(mode %in% c("exon_id","both")){
    exon_id_uniq = names(table(data$exon_id))

    out = data %>% dplyr::select(cell,gene,exon_id,size,cluster,count,polyA)

    out$cell = as.numeric(as.factor(out$cell))
    out$gene = as.numeric(as.factor(out$gene))
    out$exon_id = as.numeric(as.factor(out$exon_id))

    write.table(exon_id_uniq,file = paste(path,"exon_id.tsv",sep = "/"),
                row.names = FALSE,col.names = FALSE,quote = FALSE)
    write.table(out,file = paste(path,"exon_id_count.tsv",sep = "/"),
                row.names = FALSE,col.names = TRUE,quote = FALSE)
  }
  if(mode %in% c("transcript","both")){
    if(is.null(gene_bed) | is.null(gtf)){
      stop("To save the data as transcript count, the gtf and gene bed annotation should be provided!")
    }
    data = transcript_impute(data,gtf = gtf,gene_bed = gene_bed,...)
    transcript_uniq = names(table(data$transcript))

    out = data %>% dplyr::select(cell,gene,transcript,count)
    out$cell = as.numeric(as.factor(out$cell))
    out$gene = as.numeric(as.factor(out$gene))
    out$transcript = as.numeric(as.factor(out$transcript))

    write.table(transcript_uniq,file = paste(path,"transcript.tsv",sep = "/"),
                row.names = FALSE,col.names = FALSE,quote = FALSE)
    write.table(out,file = paste(path,"transcript_count.tsv",sep = "/"),
                row.names = FALSE,col.names = TRUE,quote = FALSE)
  }

  return(NULL)
}
