#' @title BarcodeMatchUnit
#'
#' @description Find the cell barcode within the softclips given the cell barcode whitelist
#' @details Find the cell barcode within the softclips given the cell barcode whitelist and the extract the UMI
#' sequence after cell barcodes for further UMI deduplication
#'
#' @param seq A vector of strings recording the tag regions for each sequence.
#' @param barcodes A character vector of cell barcode
#' @param mu The mean end location of cell barcode in the softclip, using the middle sites of the softclip as default
#' @param sigma The start variance for the end position distribution
#' @param sigma_start The variance to buffer at first to avoid sudden convergence of variance
#' @param k The length of sub string match between softclips and barcodes,
#' only softclips and barcodes with high enough kmer match will be further validated
#' @param batch The number of softclips to search in one run
#' @param top Barcodes with top n cos similarity will be used for further validation
#' @param cos_thresh The minimum threshold of the cos similariy to preserve the cell barcode for further validation
#' @param alpha The probability thresh to filter out barcode match, matches
#' with too high edit distance or too far position have low probability to be correct
#' @param edit_thresh The highest tolerance of edit distance for the barcode match
#' @param UMI_len The length of the UMI sequence
#' @param UMI_flank The length of the flank to add to tolerate the insertions and deletions around UMI
#' @import Longcellsrc
#' @return A dataframe contains matched cell barcode, the position of the barcode, UMI and the edit distance of the match.
#'
BarcodeMatchUnit = function(seq, barcodes,
                            mu = 20, sigma = 10, sigma_start = 10,
                            k = 6, batch = 100,top = 5, cos_thresh = 0.25, alpha = 0.05,
                            edit_thresh = 3,mean_edit_thresh = 1.5,
                            UMI_len = 10,UMI_flank = 1){
  barcode = barcodeMatch(seq, barcode,
                         mu, sigma, sigma_start,
                         k, batch,
                         top, cos_thresh, alpha, edit_thresh,
                         UMI_len, UMI_flank)

  return(data)
}

#' @title BarcodeFilter
#'
#' @description Filter out barcode match with low quality
#' @details Filter out barcode match with low quality according to the mean edit distance of the cell barcode
#'
#' @param data The output from the barcode match step. Should be a dataframe contains the read name, cell barcode and its matching edit distance
#' @param mean_edit_thresh The maximum tolerance of mean edit distance for a cell barcode to be preserved
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
BarcodeFilter = function(data,mean_edit_thresh = 1.5){
  #data = data %>% arrange(edit,-nchar(isoform))
  #data <- data[!duplicated(data$qname),]

  bc_edit = data %>% group_by(barcode) %>% summarise(edit = mean(edit))
  preserve = bc_edit[bc_edit$edit < mean_edit_thresh,]$barcode
  data = data[data$barcode %in% preserve,]
  #data = data %>% select(barcode,gene,isoform,umi,polyA,edit)
  #data = data %>% dplyr::select(qname,barcode,gene,isoform,end,softclip,umi,polyA,edit)

  return(data)
}


#' @title BarcodeMatch
#'
#' @description BarcodeMatch with parallelization
#' @details BarcodeMatch for large dataset which is paralleled by future.apply.
#' Large data will be split according to cores.
#'
#' @inheritParams BarcodeMatchUnit
#' @inheritParams BarcodeFilter
#' @param cores The number of cores to use for parallelization
#' @import Longcellsrc
#' @importFrom future.apply future_lapply
#' @importFrom parallel detectCores
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @export
BarcodeMatch = function(seq, barcodes,
                        mu = 20, sigma = 10, sigma_start = 10,
                        k = 6, batch = 100,top = 5, cos_thresh = 0.25, alpha = 0.05,
                        edit_thresh = 3,mean_edit_thresh = 1.5,
                        UMI_len = 10,UMI_flank = 1,cores = 1){
  n = length(seq)

  seq = as.data.frame(cbind(1:length(seq),seq))
  colnames(seq) = c("id","seq")
  seq = seq %>% mutate(is = as.numeric(id))

  seq_split <- dataSplit(seq,cores)

  bm = future_lapply(seq_split,function(x){
    out = barcodeMatch(x$seq, barcodes,
                       mu, sigma, sigma_start, k, batch,
                       top, cos_thresh, alpha, edit_thresh,
                       UMI_len, UMI_flank)
    out$id = x$id[out$id + 1]
    return(out)
  },future.packages = c("Longcellsrc"),future.seed=TRUE)

  bm = as.data.frame(do.call(rbind,bm))
  cat(nrow(bm)," out of ",n, "reads are identified with a vaild cell barcode.\n")
  #bm = bm[,-which(colnames(bm) == "softclip")]

  # low quality barcode alignment filtering
  bm = BarcodeFilter(data = bm,mean_edit_thresh = mean_edit_thresh)
  cat("After filtering, ",nrow(bm),
      " reads with barcodes are preserved, with mean edit distance as ",mean(bm$edit),".\n")
  return(bm)
}
