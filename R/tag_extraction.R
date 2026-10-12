#' @title seq_config
#' @description Preset parameters for different sequencing protocol
#' @param protocol The sequence protocol
#' @param toolkit The kit to build the library, should be 5 or 3
#' @param adapter The sequence of the adapter which is aside the cell barcode
#' @param left_flank,right_flank The length of the left/right part aside the adapter to be extracted
#' @param drop_adapter A flag to indicate if the adapter region should be removed from the tag region.
#' @param barcode_len,UMI_len The length of the cell barcode/UMI
seq_config = function(protocol,toolkit,
                      adapter,
                      left_flank,right_flank,drop_adapter,
                      barcode_len,UMI_len){
  if(protocol == "10X"){
    if(toolkit == 5){
      adapter = "TTTCTTATATGGG"
      UMI_len = 10
    }
    else{
      adapter = "GCGTCGTGTAG"
      UMI_len = 12
    }
    left_flank = 55
    right_flank = -10
    drop_adapter =FALSE
    barcode_len = 16

  }
  else if(protocol == "VISIUM"){
    if(toolkit == 5){
      adapter = "TTTCTTATATGGG"
    }
    else{
      adapter = "GCGTCGTGTAG"
    }
    left_flank = 55
    right_flank = -10
    drop_adapter =FALSE
    barcode_len = 16
    UMI_len = 12
  }
  else if(protocol == "Curio"){
    adapter = "TCTCGGGAACGCTGAAGA"
    left_flank = 25
    right_flank = 25
    drop_adapter =TRUE
    barcode_len = 14
    UMI_len = 9
  }
  else if(protocol == "other"){
    if(is.null(adapter)){
      stop("Please sepcify the adapter, left_flank, right_flank for your protocol to search and extract the adapter sequence!")
    }
  }
  else{
    stop("The protocol you specified is not incorporated yet, please specify your config!")
  }

  out = list(adapter,left_flank,right_flank,drop_adapter,barcode_len,UMI_len)
  names(out) = c("adapter","left_flank","right_flank","drop_adapter","barcode_len","UMI_len")
  return(out)
}


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
#' @param left_flank,right_flank The length of the left/right part aside the adapter to be extracted
#' @param drop_adapter A flag to indicate if the adapter region should be removed from the tag region.
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
#' @param barcode_len,UMI_len The length of the cell barcode/UMI
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
                        toolkit,protocol = "10X", adapter = NULL,
                        window = 10,step = 2,
                        left_flank = 0, right_flank = 0, drop_adapter = FALSE,
                        polyA_bin = 20,polyA_base_count = 15,polyA_len = 10,
                        # parameters for barcode match
                        barcode_len = 16,mu = 20, sigma = 10, k = 6, batch = 100,
                        top = 5, cos_thresh = 0.25, alpha = 0.05,
                        edit_thresh = 3,mean_edit_thresh = 1.5,
                        UMI_len = 10, UMI_flank = 1,
                        # parameter for parallel
                        cores = 1){
  fastq_format(fastq_path)

  config = seq_config(protocol,toolkit,
                      adapter,
                      left_flank,right_flank,drop_adapter,
                      barcode_len,UMI_len)
  barcode_len = config$barcode_len
  UMI_len = config$UMI_len
  barcode = read_barcodes(barcode_path,barcode_len)

  reads = extractTagFastq(fastq_path,out_name,
                          config$adapter,toolkit,window,step,
                          config$left_flank,config$right_flank,config$drop_adapter,
                          polyA_bin,polyA_base_count,polyA_len)

  if(length(reads) == 0 || nrow(reads) == 0){
    stop("No usable adapter/tag regions were found. Use untrimmed reads and verify toolkit (3 or 5), protocol and adapter; trimmed-away barcodes cannot be recovered.")
  }

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
  if(!nrow(bc)) stop("No barcodes with complete UMIs passed filtering. Check whitelist, chemistry and untrimmed input.")
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
adapter_dis = function(data){
  len = unique(nchar(data$umi))
  if(length(len) != 1L || !is.finite(len) || len == 0L){
    stop("Cannot estimate adapter quality: no consistent, non-empty UMIs remain. Check untrimmed input, chemistry, and whitelist.")
  }
  adapter = data %>%
            filter(edit == 0) %>%
            filter(nchar(adapter) == len)

  adapter_count = table(adapter$adapter)
  if(!length(adapter_count)){
    stop("Cannot estimate adapter quality: no exact barcode matches with a complete adjacent adapter. Use untrimmed reads and check toolkit/protocol and sample-specific barcodes.")
  }
  adapter_seq = names(adapter_count)[adapter_count == max(adapter_count)][1]

  adapter_dis = data.frame(needle = as.numeric(unlist(needleScores(adapter_seq,names(adapter_count)))),
                           count = as.numeric(adapter_count))
  adapter_dis = adapter_dis %>% group_by(needle) %>% summarise(count = sum(count))

  return(adapter_dis)
}
