#' @title isoform_correct_filter
#'
#' @description correct wrong mapped and truncated reads and filter out UMI singletons.
#' @details correct wrong mapped and truncated reads according to other complete reads in the same UMI cluster,
#' and filter out UMI singletons (UMIs enriched with sequencing errors) according to UMI cluster size.
#' @param gene_cells_cluster The UMI cluster in each cell for a gene
#' @param filter_ratio The ratio of singletons to be removed.
#' @param split The character to seperate different exons in the isoform representation.
#' @param sep The character to seperate the start and end position of an exon in the isoform representation.
#' @param splice_site_thresh The minimum threshold of occurance  for an splice site to be preserved

isoform_correct_filter <- function(gene_cells_cluster,filter_ratio,strand,
                            split = "|",sep = ","){

  gene_isoform = splice_site_table(gene_cells_cluster$isoform,
                                   # gene_cells_cluster$polyA,strand,
                                   split,sep,splice_site_thresh=0)
  #return(gene_isoform)
  if(length(gene_isoform) == 0 || nrow(gene_isoform) == 0){
    return(NULL)
  }

  gene_isoform = gene_isoform %>% dplyr::select(-id)

  if(ncol(gene_isoform) > 2){
    filter = as.data.frame(gene_isoform %>% dplyr::select(-c(start,end)))
    filter = colSums(filter,na.rm = TRUE) > 0
    gene_isoform = gene_isoform[,c(TRUE,filter,TRUE)]
  }

  # return(list(gene_cells_cluster,gene_isoform))
  gene_isoform = cells_isoform_correct(gene_cells_cluster$cell,
                                       gene_cells_cluster$cluster,
                                       gene_isoform,
                                       gene_cells_cluster$polyA)
  # return(gene_isoform)
  if(nrow(gene_isoform) == 0){
    return(gene_isoform)
  }

  gene_isoform = cells_isoforms_size_filter(cell_isoform_table = gene_isoform,
                                            ratio = filter_ratio)
  return(gene_isoform)
}

#' @title gene_umi_count
#'
#' @description get the umi count for each cell in a gene.
#' @details get the umi count for each cell in a gene, the process includes clustering for UMI and correction for the isoform.
#' @param cell_exon the input dataframe which includes cell barcode, UMI sequence, isoform and polyA status.
#' @param gene_cells_cluster The UMI cluster in each cell for a gene
#' @param qual The adapter distance from barcode matching to represent the data quality
#' @param bar The name of the column to represent cell barcode.
#' @param isoform The name of the column to represent the isoform.
#' @param polyA The name of the column to represent the polyA status.
#' @param sim_thresh the minimum threshold of needleman score to connect to different UMIs.
#' @param split The character to seperate different exons in the isoform representation.
#' @param sep The character to seperate the start and end position of an exon in the isoform representation.
#' @param splice_site_thresh The minimum threshold of occurance  for an splice site to be preserved
#' @param relation The function to describe the relationship between the mean and variance of UMI cluster size
#' @param alpha the alpha of confidence interval for the variance of UMI cluster size.
#' @param verbose if print information during UMI clustering for each cell.
#' @export
gene_umi_count <- function(cell_exon,qual,strand,bar = "barcode",
                      isoform = "isoform",polyA = "polyA",
                      sim_thresh = NULL,split = "|",sep = ",",
                      splice_site_thresh = 10,verbose = FALSE){
    colnames(cell_exon)[which(colnames(cell_exon) == bar)] = "cell"
    colnames(cell_exon)[which(colnames(cell_exon) == isoform)] = "isoform"
    colnames(cell_exon)[which(colnames(cell_exon) == polyA)] = "polyA"
    cell_exon$polyA = as.numeric(cell_exon$polyA)

    gene_isoform = splice_site_table(cell_exon$isoform,
                                     # cell_exon$polyA,strand,
                                     split,sep,splice_site_thresh)
    #return(gene_isoform)
    if(length(gene_isoform) == 0 || nrow(gene_isoform) == 0){
      return(NULL)
    }
    cell_exon = cell_exon[gene_isoform$id,]

    cells = unique(cell_exon[,"cell"])
    if(is.null(sim_thresh)){
      sim_thresh = nchar(cell_exon$umi)[1]/2+1
    }

    gene_cells_cluster <- lapply(cells,function(i){
        cell_i = cell_exon[cell_exon[,"cell"] == i,]
        if(verbose){
          cat(nrow(cell_i), " reads in cell ",i,"\n")
        }

        cell_i$cluster = 0

        # temporarily filter out large cells, will optimize later
        if(nrow(cell_i) > 400000 || length(unique(cell_i$umi)) > 20000){
          warning("Too many reads in cell ", i, "which exceeds the max limit of memory")
          return(NULL)
        }
        if(nrow(cell_i) != 1){
          cell_i$cluster = umi_cluster(cell_i$umi,
                           iso = cell_i$isoform,thresh = sim_thresh)
        }
        else{
          cell_i$cluster = 1
        }
        return(cell_i)
    })

    gene_cells_cluster = as.data.frame(do.call(rbind,gene_cells_cluster))
    # return(gene_cells_cluster)
    #cat("Clustering for UMI finished!\n")
    filter_ratio = mean(c(sum(qual[qual$needle < sim_thresh,"count"]),
                          sum(qual[qual$needle < sim_thresh+1,"count"])))/sum(qual$count)
    gene_isoform = isoform_correct_filter(gene_cells_cluster,filter_ratio,strand,
                                   split = split,sep = sep)


    return(gene_isoform)
}

#' @title umi_count
#'
#' @description get the umi count for each cell for each gene.
#' @details get the umi count for each cell for each gene, the process includes clustering
#' for UMI and correction for the isoform.
#' @inheritParams gene_umi_count
#' @param qual The adapter distance from barcode matching to represent the data quality
#' @param gene The name of column which stores the gene names
#' @export
umi_count <- function(cell_exon,qual, gene_strand,
                      bar = "barcode",gene = "gene",
                      isoform = "isoform",polyA = "polyA",
                      sim_thresh = NULL, split = "|",sep = ",",
                      splice_site_thresh = 10,verbose = FALSE){
    genes <- unique(cell_exon[,gene])
    genes_umi_count <- lapply(genes,function(i){
      if(verbose){
        cat(i,"\n")
      }
        sub_cell_exon = cell_exon[cell_exon[,gene] == i,]
        strand = unique(gene_strand[gene_strand$gene == i,"strand"])
        if(nrow(sub_cell_exon) < splice_site_thresh){
          if(verbose){
            cat("too few reads, will be filtered out\n")
          }
            return(NULL)
        }
        sub_umi_count = gene_umi_count(sub_cell_exon,qual = qual,strand = strand,
                                       bar = bar,isoform = isoform,
                                       polyA = polyA,sim_thresh = sim_thresh,
                                       split = split,sep = sep,
                                       splice_site_thresh = splice_site_thresh,
                                       verbose = verbose)
        if(is.null(sub_umi_count) || nrow(sub_umi_count) == 0){
            return(NULL)
        }
        sub_umi_count$gene = i
        return(sub_umi_count)
    })

    genes_umi_count <- do.call(rbind,genes_umi_count)
    return(genes_umi_count)
}
