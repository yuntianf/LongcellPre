#' @title init
#'
#' @description Initiation of the LongcellPre pipeline
#' @details Build the directories for the output of LongcellPre
#'
#' @param work_dir The specified directory for the output of LongcellPre
#' @export
init = function(work_dir){
  dir.create(work_dir,showWarnings = FALSE)
  dir.create(file.path(work_dir,"annotation"),showWarnings = FALSE)
  dir.create(file.path(work_dir,"BarcodeMatch"),showWarnings = FALSE)
  dir.create(file.path(work_dir,"bam"),showWarnings = FALSE)
  dir.create(file.path(work_dir,"out"),showWarnings = FALSE)

  out_path = normalizePath(work_dir)

  message = paste("The output path would be set in",out_path,sep = " ")
  return(message)
}



#' @title createAnnotation
#'
#' @description Build annotation based on the gtf or gene bed input
#' @details Build annotation based on the gtf or gene bed input
#'
#' @param gtf_path The path for the gtf annotation.
#' @param gene_bed_path The path for the gene bed annotation.
#' @param word_dir The path to save the output
#' @param bed_chr_col The name of the column indicating the chromosome in the gene bed
#' @param bed_start_col The name of the column indicating the start position of the gene in the gene bed
#' @param bed_end_col The name of the column indicating the end position of the gene in the gene bed
#' @param bed_strand_col The name of the column indicating the gene strand in the gene bed
#' @param bed_gene_col The name of the column indicating the gene name in the gene bed
#' @return A list with two dataframes, the first one is the gene bed and the second one is the gtf exon annotation
createAnnotation = function(gtf_path = NULL,gene_bed_path = NULL,work_dir = "./",
                            bed_chr_col = "chr",bed_start_col = "start",bed_end_col = "end",
                            bed_strand_col = "strand",bed_gene_col = "gene"){
  if(!is.null(gtf_path)){
    print("Start to build exon annotation based on the gtf file.")
    cache = gtf2bed(gtf_path,file.path(work_dir,"annotation/"),
                    if_store_binary = TRUE)
    gene_bed = cache[[1]]
    gtf = cache[[2]]
  }
  else if(!is.null(gene_bed_path)){
    gene_bed = read.table(gene_bed_path,header = TRUE)
    gene_bed = gene_bed[,c(bed_chr_col,bed_start_col,bed_end_col,
                           bed_strand_col,bed_gene_col)]
    colnames(gene_bed) = c("chr","start","end","strand","gene")
    write.table(gene_bed,file = file.path(work_dir,"annotation/gene_bed.txt"),sep = "\t",
                col.names = TRUE,row.names = FALSE)
    saveRDS(gene_bed,file = file.path(work_dir,"annotation/gene_bed.rds"))
    gtf = NULL
  }
  else{
    stop("The gene bed annotation is required, it should at least contains the chromsome positions for target genes and gene strand.\n Or a gtf annotation could be provided and LongcellPre would generate gene bed annotation based on gtf.")
  }
  return(list(gene_bed,gtf))
}

#' @title annotation
#'
#' @description Build annotations given if there already exists an annotation
#' @details Build annotations given if there already exists an annotation
#'
#' @inheritParams createAnnotation
#' @param overwrite The flag to indicate if the annotation file needs to be rewritten.
#' @export
annotation = function(gtf_path = NULL,gene_bed_path = NULL,work_dir = "./",overwrite = FALSE,
                      bed_chr_col = "chr",bed_start_col = "start",bed_end_col = "end",
                      bed_strand_col = "strand",bed_gene_col = "gene"){
  if(file.exists(file.path(work_dir,"annotation/gene_bed.rds"))){
    if(!overwrite){
      warning("The annotation already exists. If you want to overwrite it please set overwrite as TRUE!")
      gene_bed = readRDS(file.path(work_dir,"annotation/gene_bed.rds"))
      if(file.exists(file.path(work_dir,"annotation/exon_gtf.rds"))){
        gtf = readRDS(file.path(work_dir,"annotation/exon_gtf.rds"))
      }
      else{
        gtf = NULL
      }
      return(list(gene_bed,gtf))
    }
  }
  out = createAnnotation(gtf_path = gtf_path,gene_bed_path = gene_bed_path,work_dir = work_dir,
                         bed_chr_col = bed_chr_col,bed_start_col = bed_start_col,bed_end_col = bed_end_col,
                         bed_strand_col = bed_strand_col,bed_gene_col = bed_gene_col)
  return(out)
}


#' @title fastqMap
#'
#' @description Map the fastq file to genome via minimap2
#' @details Map the fastq file to genome via minimap2
#'
#' @param fastq The filename of the input fastq
#' @param out_name The filename of the output polished fastq
#' @param genome_path The filename of the reference genome
#' @param bed_path The filename of the reference bed for annotated splicing sites if provided,
#' BED12 file can be converted from GTF/GFF3 with ‘paftools.js gff2bed anno.gtf’.
#' @param minimap2 The path of the minimap2
#' @param samtools The path of the samtools
#' @param minimap2_thread The number of threads for minimap2 parallization
#' @param samtools_thread The number of threads for samtools parallization
#' @export
fastqMap = function(fastq,out_name,genome_path,bed_path = NULL,
                    minimap2 = "minimap2",samtools = "samtools",
                    minimap2_thread = 4,samtools_thread = 4){
  minimap_command = paste(c(minimap2,"-ax splice -uf --sam-hit-only -t",
                            minimap2_thread,genome_path,fastq),collapse = " ")
  if(!is.null(bed_path)){
    minimap_command = paste(c(minimap_command,"--junc-bed",bed_path),collapse = " ")
  }
  samtools_command1 = paste(c(samtools,"view -bS -@", samtools_thread,"-"),collapse = " ")
  samtools_command2 = paste(c(samtools,"sort - -@", samtools_thread,"-o",out_name,
                              "&&",samtools,"index",out_name),collapse = " ")

  command = paste(c(minimap_command,samtools_command1,samtools_command2),collapse = " | ")
  print("The command for mapping is:")
  print(command)
  system(command)
  return(NULL)
}

#' @title bamGeneCoverage
#'
#' @description detected genes with alignment in the bam file.
#' @details detected genes with alignment in the bam file.
#'
#' @param bam The filename of the input bam.
#' @param gene_range_bed The filename of the gene bed file, each row records the start and end position for the gene
#' @param outdir The folder to store the output file.
#' @param bedtools The path for the bedtools.
#' @return A bed file including all genes without read alignment in the bam file.
bamGeneCoverage = function(bam,gene_range_bed,outdir,bedtools = "bedtools"){
  command1 = paste(c(bedtools, "bamtobed -i", bam, "|",
                     bedtools, "merge -i - > ", file.path(outdir,"cover.bed")),collapse = " ")
  command2 = paste(c(bedtools, "subtract -a", gene_range_bed, "-b", file.path(outdir,"cover.bed"),
                     "-A >", file.path(outdir,"noncover.bed")),collapse = " ")

  system(command1)
  system(command2)

  filename = file.path(outdir,"noncover.bed")
  if(file.size(filename) == 0L){
    return(NULL)
  }
  noncover = read.table(filename)
  colnames(noncover) = c("chr","start","end","strand","gene")
  return(noncover)
}

#' @title reads_extract_bc
#'
#' @description A wrapper function to include extractTagBc, fastqMap and reads_extraction
#' @details A wrapper function to include extractTagBc, fastqMap and reads_extraction. And this function
#' would integrate the output from above functions by read name.
#'
#' @param work_dir The folder to store intermediate files and output
#' @inheritParams extractTagBc
#' @inheritParams fastqMap
#' @inheritParams reads_extraction
#' @inheritParams BarcodeFilter
#' @param force_barcode_match flag to force to redo the barcode match step
#' @param force_map flag to force to redo the mapping step
#' @param force_isoform_extract flag to force to redo the isoform extraction step
#' @import Longcellsrc
#' @importFrom BSgenome getBSgenome
#' @importFrom peakRAM peakRAM
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr inner_join
#' @return A list with two dataframes: the first one records the barcode, UMI and isoform information
#' for each read, the second one stores the needleman score distirbution for adapters as evaluation of
#' data quality.
#' @export
reads_extract_bc = function(fastq_path,barcode_path,
                            gene_bed,adapter,
                            genome_path,genome_name,
                            toolkit,protocol = "10X",
                            minimap_bed_path = NULL,
                            work_dir = "./",
                            # parameters for tag extraction
                            window = 10,step = 2,
                            left_flank = 0, right_flank = 0, drop_adapter = FALSE,
                            polyA_bin = 20,polyA_base_count = 15,polyA_len = 10,
                            # parameters for barcode match
                            barcode_len = 16,mu = 15, sigma = 10, k = 6, batch = 100,
                            top = 5, cos_thresh = 0.25, alpha = 0.05, edit_thresh = 3,
                            UMI_len = 10, UMI_flank = 1,
                            # parameters for fastq mapping
                            minimap2 = "minimap2",samtools = "samtools",
                            # parameters for reads extraction
                            bedtools = "bedtools",
                            map_qual = 30,end_flank = 200,
                            splice_site_bin = 2,
                            #parameters for barcode filtering
                            mean_edit_thresh = 1.5,
                            # parameters for cache
                            force_barcode_match = FALSE,
                            force_map = FALSE,
                            force_isoform_extract = FALSE,
                            # parameters for parallel
                            cores = 1){
  cores = coreDetect(cores)

  # barcode match
  cat("Start to do barcode match:")
  do_bc_flag = TRUE
  if(file.exists(file.path(work_dir,"BarcodeMatch/BarcodeMatch.txt")) &
     file.exists(file.path(work_dir,"polish.fq.gz"))){
    if(!force_barcode_match){
      warning("The barcode match output already exist,
              if you want to redo it please set force_barcode_match to be TRUE")
      bc = read.table(file.path(work_dir,"BarcodeMatch/BarcodeMatch.txt"),header = TRUE,sep = "\t")
      do_bc_flag = FALSE
    }
  }
  if(do_bc_flag){
    start_time <- Sys.time()
    mem = peakRAM::peakRAM({
      bc = extractTagBc(fastq_path = fastq_path,barcode_path = barcode_path,
                        out_name = file.path(work_dir,"polish.fq.gz"),
                        # parameters to extract the tag region
                        toolkit = toolkit,protocol = protocol,adapter = adapter,
                        window = window,step = step,
                        left_flank = left_flank, right_flank = right_flank,
                        drop_adapter = drop_adapter,
                        polyA_bin = polyA_bin,
                        polyA_base_count = polyA_base_count,polyA_len = polyA_len,
                        # parameters for barcode match
                        barcode_len = barcode_len,mu = mu, sigma = sigma,
                        k = k, batch = batch,top = top, cos_thresh = cos_thresh,
                        alpha = alpha, edit_thresh = edit_thresh,
                        UMI_len = UMI_len, UMI_flank = UMI_flank,
                        # parameter for parallel
                        cores = cores)
    })
    saveResult(bc,file.path(work_dir,"BarcodeMatch/BarcodeMatch.txt"))
    end_time <- Sys.time()
    duration = end_time-start_time
    log = sprintf('Barcode match took %.2f %s\n', duration, units(duration))
    cat(log,"\n")
    print(mem[,2:4])
  }


  # fastq mapping
  cat("Start to map polished fastq to genome:")
  do_map_flag = TRUE
  if(file.exists(file.path(work_dir,"bam/polish.bam"))){
    if(!force_map){
      warning("The mapping result already exists,
              if you want to redo it please set force_map to be TRUE")
      do_map_flag = FALSE
    }
  }
  if(do_map_flag){
    start_time <- Sys.time()
    cache = fastqMap(fastq = file.path(work_dir,"polish.fq.gz"),
                     out_name = file.path(work_dir,"bam/polish.bam"),
                     genome_path = genome_path,bed_path = minimap_bed_path,
                     minimap2 = minimap2,samtools = samtools,
                     minimap2_thread = cores,samtools_thread = cores)
    end_time <- Sys.time()
    duration = end_time-start_time
    sprintf('Genome mapping took %.2f %s\n', duration, units(duration))
  }

  # isoform extraction
  cat("Start to extract isoforms:")
  do_extract_flag = TRUE
  if(file.exists(file.path(work_dir,"BarcodeMatch/BarcodeMatchIso.txt")) &
     file.exists(file.path(work_dir,"BarcodeMatch/adapterNeedle.txt"))){
    if(!force_isoform_extract){
      warning("The barcode match and isoform extraction result already exist,
              if you want to redo it please set force_isoform_extract to be TRUE")
      reads_bc = read.table(file.path(work_dir,"BarcodeMatch/BarcodeMatchIso.txt"),header = TRUE,sep = "\t")
      qual = read.table(file.path(work_dir,"BarcodeMatch/adapterNeedle.txt"),header = TRUE,sep = "\t")
      do_extract_flag = FALSE
    }
  }

  if(do_extract_flag){
    gene_range = gene_bed %>% group_by(gene) %>% summarise(chr = unique(chr),start = min(start),
                                                           end = max(end),strand = unique(strand))
    gene_range = gene_range[,c("chr","start","end","strand","gene")]
    write.table(gene_range,file.path(work_dir,"annotation/gene_range.txt"),sep = "\t",quote = FALSE,
                row.names = FALSE,col.names = FALSE)

    noncover = bamGeneCoverage(bam = file.path(work_dir,"bam/polish.bam"),
                                 gene_range_bed = file.path(work_dir,"annotation/gene_range.txt"),
                                 outdir = file.path(work_dir,"annotation"),
                                 bedtools = bedtools)
    if(!is.null(noncover)){
      gene_bed = gene_bed %>% filter(!gene %in% noncover$gene)
    }

    start_time <- Sys.time()
    mem = peakRAM::peakRAM({
      genome<-BSgenome::getBSgenome(genome_name)
      reads = reads_extraction(bam_path = file.path(work_dir,"bam/polish.bam"),
                               gene_bed = gene_bed,genome = genome,
                               toolkit = toolkit,
                               map_qual = map_qual,end_flank = end_flank,
                               splice_site_bin = splice_site_bin)

      reads_bc = inner_join(bc,reads,by = c("name" = "qname"))
      reads_bc = reads_bc %>%
        mutate(polyA.x = as.numeric(polyA.x),polyA.y = as.numeric(polyA.y)) %>%
        mutate(polyA = polyA.x & polyA.y) %>% dplyr::select(-polyA.x,-polyA.y)
    })
    end_time <- Sys.time()
    duration = end_time-start_time
    log = sprintf('Isoform extraction took %.2f %s\n', duration, units(duration))
    cat(log,"\n")
    print(mem[,2:4])

    if(nrow(reads_bc) > 0){
      # evaluate data quality
      qual = adapter_dis(data = reads_bc,UMI_len = UMI_len,flank = UMI_flank)

      #reads_bc = reads_bc %>% dplyr::select(qname,barcode,gene,isoform,umi,polyA)
      saveResult(reads_bc,file.path(work_dir,"BarcodeMatch/BarcodeMatchIso.txt"))
      saveResult(qual,file.path(work_dir,"BarcodeMatch/adapterNeedle.txt"))
    }
    else{
      stop("No read is found with valid barcode, please check if your barcode and fastq file match!")
    }
  }

  return(list(reads_bc,qual))
}


#' @title umi_count_corres
#'
#' @description A wrapper function to include all steps in the UMI deduplication.
#' @details A wrapper function to include all steps in the UMI deduplication, including UMI cluster,
#' isoform correction, isoform imputation and quantification.
#'
#' @param data The first output dataframe from reads_extract_bc.
#' @inheritParams umi_count
#' @inheritParams cells_genes_isos_count
#' @inheritParams genes_distribute
#' @importFrom dplyr bind_rows
#' @import Longcellsrc
#' @importFrom future.apply future_lapply
#' @importFrom peakRAM peakRAM
#' @export
#'

umi_count_corres = function(data,qual,dir,gene_bed,gtf = NULL,
                            # parameter for umi count
                            bar = "barcode",gene = "gene",
                            isoform = "isoform",polyA = "polyA",
                            sim_thresh = NULL, split = "|",sep = ",",
                            splice_site_thresh = 10,verbose = FALSE,
                            bed_gene_col = "gene",bed_strand_col = "strand",
                            #parameter for mapping to transcript
                            to_isoform = TRUE,filter_only_intron = TRUE,
                            mid_offset_thresh = 3,overlap_thresh = 0,
                            gtf_gene_col = "gene",gtf_start_col = "start",
                            gtf_end_col = "end",gtf_iso_col = "transname",
                            #parameter for parallel
                            cores = 1){
  cores = coreDetect(cores)
  data_split = genes_distribute(data,cores,gene)

  gene_strand = unique(gene_bed[,c(bed_gene_col,bed_strand_col)])

  cat("Start to do UMI deduplication:")
  start_time <- Sys.time()
  mem = peakRAM::peakRAM({
    count = future_lapply(data_split,function(x){
      sub_count = umi_count(x,qual,gene_strand,
                            bar = bar,gene = gene,
                            isoform = isoform,polyA = polyA,
                            sim_thresh = sim_thresh,
                            split = split,sep = sep,
                            splice_site_thresh = splice_site_thresh,
                            verbose = verbose)

      if(length(sub_count) == 0 || nrow(sub_count) == 0){
        return(NULL)
      }
      return(sub_count)
    },future.packages = c("Longcellsrc"),future.seed=TRUE)
  })
  end_time <- Sys.time()
  duration = end_time-start_time
  log = sprintf('UMI deduplication took %.2f %s\n', duration, units(duration))
  cat(log,"\n")
  print(mem[,2:4])

  if(to_isoform){
    if(!is.null(gtf)){
      cat("Start to do isoform alignment:")
      start_time <- Sys.time()
      mem = peakRAM::peakRAM({
      count_mat = future_lapply(count,function(x){
        if(is.null(x)){
          return(NULL)
        }
        sub_count_mat = cells_genes_isos_count(x,gtf,
                                               thresh = mid_offset_thresh,
                                               overlap_thresh = overlap_thresh,
                                               filter_only_intron = filter_only_intron,
                                               gtf_gene_col = gtf_gene_col,
                                               gtf_iso_col = gtf_iso_col,
                                               gtf_start_col = gtf_start_col,
                                               gtf_end_col = gtf_end_col,
                                               split = split,sep = sep)
        if(length(sub_count_mat) == 0 || nrow(sub_count_mat) == 0){
          return(NULL)
        }
        return(sub_count_mat)
      },future.packages = c("Longcellsrc"),future.seed=TRUE)
        count_mat = as.data.frame(do.call(dplyr::bind_rows,count_mat))
        #count_mat[is.na(count_mat)] = 0
        #saveResult(count_mat,file.path(dir,"iso_count_mat.txt"))
        saveIsoMat(count_mat,dir)
      })
      end_time <- Sys.time()
      duration = end_time-start_time
      log = sprintf('Isoform alignment took %.2f %s\n', duration, units(duration))
      cat(log,"\n")
      print(mem[,2:4])
    }
    else{
      warning("The gtf annotation is not provided for the isoform imputation, will skip this step!")
    }
  }

  count = as.data.frame(do.call(rbind,count))
  count = count %>% dplyr::select(cell,gene,isoform,count,polyA)
  saveResult(count,file.path(dir,"iso_count.txt"))
  return(0)
}


#' @title RunLongcellPre
#'
#' @description A wrapper function to include the reads_extract_bc and umi_count_corres.
#' @details A wrapper function to include the reads_extract_bc and umi_count_corres.
#'
#' @inheritParams reads_extract_bc
#' @inheritParams umi_count_corres
#' @param force_barcode_match decide whether to redo the barcode match if there already exist the output
#' @param ... All other tuning parameters
#'
#' @importFrom future plan
#' @importFrom future sequential
#' @importFrom future multisession
#' @importFrom future multicore
#' @export
#'
RunLongcellPre = function(fastq_path,barcode_path,
                          toolkit,
                          genome_path,genome_name,
                          protocol = "10X",
                          adapter = NULL,
                          gtf_path = NULL,
                          gene_bed_path = NULL,
                          minimap_bed_path = NULL,
                          to_isoform = TRUE,
                          work_dir = "./",cores = 1,
                          mode = c("sequential","multisession","multicore","cluster"),
                          minimap2 = "minimap2",
                          samtools = "samtools",
                          bedtools = "bedtools",
                          ...){
  # initialization
  cache = init(work_dir)

  # annotation
  neceParam = list(gtf_path = gtf_path,gene_bed_path = gene_bed_path,work_dir = work_dir)
  Param = paramMerge(func = annotation,neceParam = neceParam,...)
  cache = do.call(annotation,Param)
  gene_bed = cache[[1]]
  gtf = cache[[2]]

  # parallel
  cores = coreDetect(cores)
  mode = match.arg(mode,c("sequential","multisession","multicore","cluster"))

  if(cores > 1 & mode == "sequential"){
    warning("Set more than one core to use, will use multisession mode instead!")
    mode = "multisession"
  }
  if(mode == "sequential"){
    plan(strategy = mode)
  }
  else{
    plan(strategy = mode,workers = cores)
  }
  cat("LongcellPre would be applied with ",cores," threads in ",mode," mode.\n")

  # barcode match and reads extraction
  neceParam = list(fastq_path = fastq_path,barcode_path = barcode_path,
                   gene_bed = gene_bed,adapter = adapter,
                   genome_path = genome_path,genome_name = genome_name,
                   toolkit = as.numeric(toolkit),protocol = protocol,
                   minimap_bed_path = minimap_bed_path,
                   work_dir = work_dir,
                   minimap2 = minimap2,samtools = samtools,bedtools = bedtools,
                   cores = cores)
  Param = paramMerge(reads_extract_bc,neceParam,...)
  bc = do.call(reads_extract_bc,Param)
  bc_out = bc[[1]]
  qual = bc[[2]]

  # UMI deduplication
  neceParam = list(data = bc_out,qual = qual,
                   dir = file.path(work_dir,"out"),
                   gene_bed = gene_bed,gtf = gtf,
                   to_isoform = to_isoform,
                   cores = cores)
  Param = paramMerge(umi_count_corres,neceParam,...)
  uc = do.call(umi_count_corres,Param)


  print("LongcellPre pipeline finished!")
}


