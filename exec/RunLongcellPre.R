#!/usr/bin/env Rscript
library(LongcellPre)
library(argparse)


argparse = function(){

  args <- commandArgs(trailingOnly = TRUE)  # Fetch only arguments passed to the script
  help = "-h" %in% args || "--help" %in% args
  full_help <- help & "--full" %in% args

  p = ArgumentParser(description = "Pipeline for single cell Nanopore RNA-seq preprocessing")

  ### Necessary parameters
  ##### sequencing related args
  p$add_argument("-f", "--fastq", help="The path for the input fastq file",required =TRUE)
  p$add_argument("-b", "--barcode", help="The path for the cell barcode whitelist",required =TRUE)
  p$add_argument("-t", "--toolkit", help="The toolkit used in sequencing, should be 5 or 3",choices = c("5","3"),required =TRUE)
  p$add_argument("-q", "--protocol", help="The sequencing protocol, ex. '10X', 'VISIUM', 'Curio'",
                 choices = c("10X","VISIUM","Curio","other"),default = "10X")

  ##### reference related args
  p$add_argument("-g", "--genome_path", help="The path of the genome reference", required =TRUE)
  p$add_argument("-n", "--genome_name", help="the genome name used for mapping, ex. hg38", required =TRUE)

  p$add_argument("--gtf", help="The path of the gtf annotation",required = FALSE)
  p$add_argument("--gene_bed_path", help="The path of the gene bed annotation",required = FALSE)
  p$add_argument("--minimap_bed_path", help="The path of your bed annotation for minimap2, can be generated from GTF/GFF3 with ‘paftools.js gff2bed anno.gtf’",
                 required = FALSE)

  ##### output related args
  p$add_argument("-o","--work_dir", help="The output directory", default = "./")
  p$add_argument("--to_isoform", help="A flag to indicate if the cell by isoform matrix should be generated", default = TRUE)

  ##### parallel related args
  p$add_argument("-c","--cores", help="The number of cores used for parallization", default = 1)
  p$add_argument("-m","--mode", help="The mode for parallization. The parallization is implemented with future.apply,
                 the feasible modes can be 'sequential','multicore','cluster'",
                 choices = c("sequential","multisession","multicore","cluster"), default = "sequential")

  ##### tool related args
  p$add_argument("--minimap2", help="The path of the minimap2", default = "minimap2")
  p$add_argument("--samtools", help="The path of the samtools", default = "samtools")
  p$add_argument("--bedtools", help="The path of the bedtools", default = "bedtools")


  if(help & !full_help){
      message("Here only shows the necessary parameters to run this pipeline. To get a full list of parameters, please run Rscript RunLongcellPre.R -h --full")
  }

  if (full_help | !help){
    ### Unnecessary parameters
    ##### tag extraction
    p$add_argument("--window", help="The kmer size used to search for adapter sequence", default = 10)
    p$add_argument("--step", help="The step size when search the adapter in kmer way", default = 2)
    p$add_argument("--left_flank", help="After the adapter is found, the length of the left flank sequence to be preserved", default = 0)
    p$add_argument("--right_flank", help="After the adapter is found, the length of the right flank sequence to be preserved", default = 0)
    p$add_argument("--drop_adapter", help="After the adapter is found, the molecular tag aside the adapter would be extracted and returned. In this step, 'drop_adapter' indicates
                 if the adapter sequence should be droped.", default = FALSE)
    p$add_argument("--polyA_bin", help="The window size to search for the polyA", default = 20)
    p$add_argument("--polyA_base_count", help="The minimum number of base A's required within the search window to confirm the presence of a polyA sequence", default = 15)

    # parameters for barcode match
    p$add_argument("--barcode_len", help="The length of the cell barcode", default = 16)
    p$add_argument("--batch", help="The number of reads to be processed for cell barcode alignment for one time.
                   After one round of process, the parameters for the distribution of the start position of cell barcode would be updated", default = 100)
    p$add_argument("--cos_thresh", help="The lower bound of the cos similarity for the cell barcodes to be preserved as candidates for another round of
                   edit distance comparison", default = 0.25)
    p$add_argument("--top", help="If there are multiple cell barcodes passing the threshold of cos similarity, the number of top barcodes to be preserved for the
                   edit distance comparison", default = 5)
    p$add_argument("--edit_thresh", help="The higher bound of the edit distance for the cell barcode match", default = 3)
    p$add_argument("--UMI_len", help="The length of the unique molecule identifier", default = 10)
    p$add_argument("--UMI_flank", help="To be tolerant of the insertions and deletions, the length of the flank sequence when extracting the UMI", default = 1)
    p$add_argument("--mean_edit_thresh", help="After cell barcode matching, any cell barcode with an average edit distance across the mapped reads exceeding the mean_edit_thresh will be filtered out.", default = 1.5)

    # parameters for reads extraction
    p$add_argument("--map_qual", help="The lower bound for the mapping quality of a read", default = 30)
    p$add_argument("--end_flank", help="The maximum allowable length for a read to exceed the boundaries of a gene during mapping.", default = 200)

    # parameters for rerun
    p$add_argument("--force_barcode_match", help="The flag to indicate if the cell barcode match should be rerun if there already exist the output files", default = FALSE)
    p$add_argument("--force_map", help="The flag to indicate if the read mapping should be rerun if there already exist the output files", default = FALSE)
    p$add_argument("--force_isoform_extract", help="The flag to indicate if the isoform extraction should be rerun if there already exist the output files", default = FALSE)
    p$add_argument("--force_UMI_dedup", help="The flag to indicate if the UMI deduplication should be rerun if there already exist the output files", default = FALSE)
    p$add_argument("--force_fastq_out", help="The flag to indicate if the fastq output should be rerun if there already exist the output files", default = FALSE)

    # parameters for reads filtering
    p$add_argument("--splice_site_thresh", help="If the number of the appearance for a splice site is lower than this threshold, reads with this splice site would be filtered out", default = 10)
    p$add_argument("--filter_only_intron", help="A bool flag to indicate if reads only introns should be preserved", default = TRUE)

    # parameters for mapping reads to isoform
    p$add_argument("--mid_offset_thresh", help="When mapping a read to an annotated isoform, the maximum allowable offset for an internal splice site.", default = 3)
    p$add_argument("--overlap_thresh", help="When mapping a read to an annotated isoform, the read must have a minimum overlap ratio with the annotated isoform", default = 0)

    # parameters for verbose
    p$add_argument("--verbose",help = "The flag to indicate the information print of umi count process", default = FALSE)
  }

  return(p$parse_args())
}

args = argparse()

cache = do.call(RunLongcellPre, args)
