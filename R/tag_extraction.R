#' @title extractTagBc
#'
#' @description Extract the tag region and identify cell barcode in the ta region
#' @details Extract the tag region for each reads and trim the tag region out to generate a
#' polished read. Then the tag region is used to identify cell barcode an UMI.
#'
#' @param fastq_path The path of the input fastq
#' @param out_name The path for the polished fastq
#' @param barcode_path The path for the barcode whitelist
#' @param toolkit The kit to build the library, should be 5 or 3
#' @param adapter The sequence of the adapter which is aside the cell barcode
#' @param window The window size to search the substring of adapter in the read
#' @param step The step size to search the substring of adapter in the read
#' @param len The length of tag region to be extratced
#' @param polyA_bin The window size to search for polyA.
#' @param polyA_base_count The minimum threshold of the number of A within the polyA search window
#' @param polyA_len The maximum length of polyA to be preserved in the read
#' @param mu The initial expectation of the start positions for cell barcodes in the tag region
#' @param sigma The initial variance of the start positions for cell barcodes in the tag region
#' @param k The length of barcode substring search in the tag region
#' @param batch The number of reads to be processed to update the distribution of
#' start position of cell barcodes
#' @param top The number of candidate barcodes which are further compared with edit distance
#' @param cos_thresh The minimum threshold for the cosine similarity between the barcode and
#' the tag region, only barcode with cosine similarity larger than the thresh can be preserved as
#' candidate barcodes.
#' @param alpha The size of the confidence interval for the start position of cell barcode alignment
#' @param edit_thresh The maximum threshold for the edit distance of the barcode alignment, alighment with
#' edit distance over the threshold would be filtered out.
#' @param UMI_len The length of the UMI
#' @param flank The length of flank when extract UMI, which is used to be tolerant of indels and dels.
#' @param cores The number of cores to use for parallization
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom stringi stri_reverse
#' @importFrom spgs reverseComplement
#' @importFrom future sequential
#' @importFrom future multisession
#' @importFrom future multicore
#' @importFrom future.apply future_lapply
#' @export
#'
extractTagBc = function(fastq_path,barcode_path,out_name,
                        # parameters to extract the tag region
                        toolkit,adapter = NULL,
                        window = 10,step = 2,len = 55,
                        polyA_bin = 20,polyA_base_count = 15,polyA_len = 10,
                        # parameters for barcode match
                        mu = 20, sigma = 10, k = 6, batch = 100,
                        top = 5, cos_thresh = 0.25, alpha = 0.05,
                        edit_thresh = 3,mean_edit_thresh = 1.5,
                        UMI_len = 10, UMI_flank = 1,
                        # parameter for parallel
                        cores = 1){
  if(is.null(adapter)){
    if(toolkit == 5){
      adapter = "TTTCTTATATGGG"
    }
    else if(toolkit == 3){
      adapter = "GCGTCGTGTAG"
    }
  }

  reads = extractTagFastq(fastq_path,out_name,
                        adapter,toolkit,window,step,len,
                        polyA_bin,polyA_base_count,polyA_len)

  barcode = read.table(barcode_path)[,1]
  barcode = substr(barcode,1,16)

  if(toolkit == 5){
    reads$tag = stringi::stri_reverse(reads$tag)
    barcode = stringi::stri_reverse(barcode)
  }
  else if(toolkit == 3){
    barcode = spgs::reverseComplement(barcode,case = "upper")
  }

  bc = BarcodeMatch(reads$tag, barcode,
                    mu = mu, sigma = sigma, sigma_start = 10,
                    k = k, batch = batch,top = top, cos_thresh = cos_thresh, alpha = alpha,
                    edit_thresh = edit_thresh,mean_edit_thresh = mean_edit_thresh,
                    UMI_len = UMI_len,UMI_flank = UMI_flank,cores = cores)

  if(toolkit == 5){
    bc$barcode = stringi::stri_reverse(bc$barcode)
  }
  else if(toolkit == 3){
    bc$barcode = spgs::reverseComplement(bc$barcode,case = "upper")
  }

  bc = bc %>% filter(nchar(umi) == UMI_len+2*UMI_flank)
  orig_col = setdiff(colnames(bc),"id")
  bc = cbind(bc,reads[bc$id,c("name","tag","polyA")]) %>% dplyr::select(-id)

  bc = bc[,c("name","tag",orig_col,"polyA")]
  return(bc)
}


#' @title adapter_dis
#'
#' @description Calculate the needleman score between original adapter and adapter in reads.
#' @details Calculate the edit distance between original adapter and adapter in reads to represent the data quality.
#' Only adapters upstream confidently identified barcodes will be used.
#'
#' @param data The output from the barcode match step. Should be a dataframe contains edit distance and
#' the adapter sequence
#' @param flank The length of flank when extract UMI, which is used to be tolerant of indels and dels.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr filter
#' @importFrom dplyr summarise
#' @importFrom NameNeedle needleScores
#' @return A dataframe with two columns, the first is the edit distance and the secons is
#' how many adapters have such an edit distance to the original adapter
adapter_dis = function(data,UMI_len = 10,flank = 1){
  adapter = data %>% filter(edit == 0) %>%
    filter(nchar(adapter) == UMI_len+2*flank)

  adapter_count = table(adapter$adapter)
  adapter_seq = names(adapter_count)[adapter_count == max(adapter_count)]

  adapter_dis = cbind(needleScores(adapter_seq,names(adapter_count)),adapter_count)
  adapter_dis = as.data.frame(adapter_dis)
  colnames(adapter_dis) = c("needle","count")
  adapter_dis = adapter_dis %>% group_by(needle) %>% summarise(count = sum(count))

  return(adapter_dis)
}
