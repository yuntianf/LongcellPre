#' @title extractReadSeq
#' @description extract the sequence from the parent sequence given the coordinates.
#' @param seq The string of the parent sequence
#' @param start,end The start and the end position of the parent sequence in the genome
#' @param exon_bins A matrix with two columns, each row is a bin and each element represents the start and end position of the bin.
#' @return A string of the extracted sequence.
extractReadSeq = function(seq,start,end,exon_bins){
  coord = exon_bins-start+1

  if(ncol(coord) != 2L || nrow(coord) == 0L || anyNA(coord) ||
     any(!is.finite(as.matrix(coord))) || any(coord < 1) ||
     any(coord > length(seq)) || any(coord[,1] > coord[,2])){
    stop("Exon coordinates must lie within the parent sequence and have start <= end.")
  }
  read = lapply(seq_len(nrow(coord)),function(i){
    start = coord[i,1]
    end = coord[i,2]

    exon = seq[start:end]
    return(exon)
  })
  read = do.call(Biostrings::xscat,read)
  return(read)
}


#' @title extractSeq
#' @description extract the sequence from the genome given the splicing sites.
#' @param genome The name of the genome, will be used to retrieve the genome object from BSgenome.
#' @param chr The chromosome index of the sequence to be extracted.
#' @param isoforms A string vector to indicate the isoforms.
#' @param strand The strand of the sequence to be extracted.
#' @param split The character to seperate the exons in the isoform string.
#' @param sep The character to seperate the start and end position of an exon in the isoform string.
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings DNAStringSet
#' @return A string of the extracted sequence.
extractSeq = function(genome,chr,isoforms,strand,split = "|",sep = ","){
  exons_list = lapply(isoforms,function(x){
    bins = read2bins(x,split = split,sep = sep)
  })

  sides = do.call(rbind,exons_list)
  start = min(sides[,1])
  end = max(sides[,2])

  whole_seq = reference_seq(genome,chr,start = start,end = end,strand = "+")

  reads = lapply(exons_list,function(x){
    read = extractReadSeq(whole_seq,start,end,x)
    return(read)
  })

  reads = Biostrings::DNAStringSet(reads)
  if(strand == "-"){
    reads = Biostrings::reverseComplement(reads)
  }
  return(reads)
}

#' @title isoformCount2ReadsForGene
#' @description Transform the isoform quantification into corresponding read sequences for a gene.
#' @param mat The dataframe which records the isoform quantification, should have at least two columns, including
#' the isoform and the count of the isoform, there can also be other columns about the meta data, which would be saved
#' in the read name.
#' @inheritParams extractSeq
#' @param isoform_col The name of the column which records the isoforms
#' @param count_col The name of the column which records the read count
#' @param name_col The name of the columns which record the meta information of each read. The information from those columns would be saved in the read name.
#' @param quality Currently the quality for the UMI-collapsed reads are all set to the highest as '~'.
#' @param read_offset Number of reads already emitted for this gene; used for stable identifiers across batches.
#' @importFrom dplyr filter_at across everything cur_column
#' @importFrom ShortRead ShortReadQ
#' @importFrom Biostrings BStringSet
#' @return A ShortReadQ object storing the fastq file.
isoformCount2ReadsForGene = function(mat,genome,chr,strand,
                                     isoform_col = "isoform",count_col = "count",name_col = c("cell","gene","polyA"),
                                     quality = "~",read_offset = 0,...){
  if(length(quality) != 1L || is.na(quality) ||
     nchar(quality, type = "bytes") != 1L ||
     utf8ToInt(quality) < 33L || utf8ToInt(quality) > 126L){
    stop("quality must be one printable FASTQ quality character (ASCII 33-126).")
  }
  mat = as.data.frame(mat)
  counts = mat[[count_col]]
  if(!is.numeric(counts) || anyNA(counts) ||
     any(!is.finite(counts)) || any(counts < 0)){
    stop("Counts must be finite, non-negative numbers.")
  }
  # Keep the existing rounding rule, including rounding ties to even.
  mat = mat[round(counts) > 0,,drop = FALSE]
  if(nrow(mat) == 0L){
    return(ShortRead::ShortReadQ(sread = Biostrings::DNAStringSet(),
                               quality = Biostrings::BStringSet(),
                               id = Biostrings::BStringSet()))
  }

  qname = as.data.frame(mat[,name_col,drop = FALSE])
  qname = qname %>% mutate(across(everything(), ~ paste0(cur_column(), "=", .)))
  qname = do.call(paste, c(qname,sep = "|"))

  seq = extractSeq(genome,chr,unlist(mat[,isoform_col]),strand,...)

  read_qual = strrep(quality, Biostrings::width(seq))


  count = unlist(round(mat[,count_col]))

  qname = rep(qname,count)
  qname = paste(qname,format(read_offset + seq_along(qname),scientific=FALSE,trim=TRUE),sep = "_")

  seq = rep(seq,count)

  read_qual = rep(read_qual,count)

  qname = BStringSet(qname)
  read_qual = BStringSet(read_qual)
  fastq = ShortReadQ(sread = seq, quality = read_qual, id = qname)

  return(fastq)
}


#' @title isoformCount2Reads
#' @description Transform the isoform quantification into corresponding read sequences for multiple genes.
#' @inheritParams isoformCount2ReadsForGene
#' @param gene_bed The gene bed annotation
#' @param filename The name of the output fastq file
#' @param gene_col,chr_col,strand_col The name of the column which stores the gene name/chromosome/gene strand.
#' @param return_reads Return all reads in memory for compatibility. Set FALSE to
#' write one gene at a time without retaining the full FASTQ in memory.
#' @param annotation_file Optional CSV path for read-name and cell annotations,
#' written alongside each gene's reads.
#' @param batch_reads,batch_bases Export batch limits. A single longer read may exceed the base limit.
#' @param qc_file Optional TSV reporting exported reads, bases and reference N bases per gene.
#' @importFrom dplyr filter_at
#' @importFrom ShortRead writeFastq append
#' @return A ShortReadQ object when return_reads is TRUE; otherwise the number
#' of written reads, invisibly. The FASTQ is written in both modes.
#' @export
isoformCount2Reads = function(mat,genome,gene_bed,filename,
                        gene_col = "gene",
                        chr_col = "chr",strand_col = "strand",
                        isoform_col = "isoform",count_col = "count",
                        name_col = c("cell","gene","polyA"),
                        quality = "~",return_reads = TRUE,annotation_file = NULL,
                        batch_reads = 10000,batch_bases = 10000000,qc_file = NULL,...){
  if(anyNA(c(batch_reads,batch_bases)) || any(!is.finite(c(batch_reads,batch_bases))) ||
     any(c(batch_reads,batch_bases) < 1)) stop("Export batch limits must be positive finite numbers.")
  batch_reads = floor(batch_reads)
  batch_bases = floor(batch_bases)
  mat = as.data.frame(mat)
  if(!is.numeric(mat[[count_col]]) || anyNA(mat[[count_col]]) ||
     any(!is.finite(mat[[count_col]])) || any(mat[[count_col]] < 0))
    stop("Counts must be finite, non-negative numbers.")
  gene_bed = as.data.frame(gene_bed)
  gene_uniq = unique(unlist(mat[,gene_col]))
  if(anyNA(gene_uniq)) stop("Gene identifiers must not be missing.")
  if(!is.null(annotation_file) && !"cell" %in% name_col){
    stop("name_col must include cell when writing annotation_file.")
  }
  if(!is.null(annotation_file) &&
     identical(normalizePath(filename, mustWork = FALSE),
               normalizePath(annotation_file, mustWork = FALSE))){
    stop("FASTQ and annotation_file must be different files.")
  }
  gene_rows = split(seq_len(nrow(mat)), match(mat[[gene_col]], gene_uniq))
  bed_rows = split(seq_len(nrow(gene_bed)), match(gene_bed[[gene_col]], gene_uniq))

  # Stage files beside their destinations; failed sequence construction leaves
  # any previous complete outputs intact.
  fastq_tmp = tempfile(".longcell-", tmpdir = dirname(filename), fileext = ".fq.gz")
  annot_tmp = if(!is.null(annotation_file))
    tempfile(".longcell-", tmpdir = dirname(annotation_file), fileext = ".csv")
  on.exit(unlink(c(fastq_tmp, annot_tmp)), add = TRUE)
  if(!is.null(annotation_file)){
    write.table(data.frame(qname = character(), annot = character()), annot_tmp,
                sep = ",", quote = FALSE, row.names = FALSE)
  }
  reads = list()
  qc = list()
  nreads = 0
  for(i in seq_along(gene_uniq)){
    key = as.character(i)
    sub_mat = mat[gene_rows[[key]],,drop = FALSE]
    sub_bed = gene_bed[bed_rows[[key]],,drop = FALSE]
    chr = unique(sub_bed[[chr_col]])
    strand = unique(sub_bed[[strand_col]])
    if(length(chr) != 1L || length(strand) != 1L ||
       is.na(chr) || is.na(strand) || !strand %in% c("+", "-")){
      stop("Gene ", gene_uniq[i], " must have exactly one chromosome and strand.")
    }
    if(!chr %in% reference_names(genome)){
      stop("Gene ", gene_uniq[i], ": chromosome ",chr," is absent from the reference genome.")
    }
    remaining = round(sub_mat[[count_col]])
    lengths = vapply(sub_mat[[isoform_col]],function(iso){
      bins = read2bins(iso,...)
      sum(bins[,2]-bins[,1]+1)
    },numeric(1))
    if(any(!is.finite(lengths)) || any(lengths < 1)) stop("Invalid isoform lengths.")
    offset = ambiguous = bases = 0
    while(any(remaining > 0)) {
      chunk_counts = numeric(length(remaining))
      nr = nb = 0
      for(j in which(remaining > 0)) {
        take = min(remaining[j],batch_reads-nr,floor((batch_bases-nb)/lengths[j]))
        # A single read longer than the base budget must still be exportable.
        if(take == 0 && nr == 0) take = 1
        if(take == 0) break
        chunk_counts[j] = take
        nr = nr + take
        nb = nb + take*lengths[j]
        if(nr >= batch_reads || nb >= batch_bases || take < remaining[j]) break
      }
      chunk = sub_mat[chunk_counts > 0,,drop=FALSE]
      chunk[[count_col]] = chunk_counts[chunk_counts > 0]
      gene_reads = isoformCount2ReadsForGene(chunk,genome,chr,strand,
                                          isoform_col = isoform_col,count_col = count_col,
                                          name_col = name_col,quality = quality,read_offset=offset,...)
    ShortRead::writeFastq(gene_reads, fastq_tmp,
                         mode = if(nreads == 0) "w" else "a", compress = TRUE)
    if(!is.null(annotation_file)){
      qname = as.character(ShortRead::id(gene_reads))
      annot = extractAnnotFromQname(qname, "cell")
      write.table(data.frame(qname, annot), annot_tmp, sep = ",", quote = FALSE,
                  row.names = FALSE, col.names = FALSE, append = TRUE)
    }
    nreads = nreads + length(gene_reads)
    offset = offset + length(gene_reads)
    bases = bases + sum(Biostrings::width(ShortRead::sread(gene_reads)))
    ambiguous = ambiguous + sum(Biostrings::letterFrequency(ShortRead::sread(gene_reads),"N"))
    remaining = remaining - chunk_counts
    if(return_reads) reads[[length(reads) + 1L]] = gene_reads
    rm(gene_reads)
    }
    qc[[i]] = data.frame(gene=gene_uniq[i],reads=offset,bases=bases,reference_N_bases=ambiguous)
  }
  if(nreads == 0) stop("No reads remain for FASTQ export after count/reference filtering.")
  staged = c(fastq_tmp,annot_tmp)
  destinations = c(filename,annotation_file)
  if(!is.null(qc_file)) {
    qc_tmp = tempfile(".longcell-qc-",tmpdir=dirname(qc_file))
    on.exit(unlink(qc_tmp),add=TRUE)
    write.table(dplyr::bind_rows(qc),qc_tmp,sep="\t",quote=FALSE,row.names=FALSE)
    staged = c(staged,qc_tmp)
    destinations = c(destinations,qc_file)
  }
  if(anyDuplicated(normalizePath(destinations,mustWork=FALSE))) stop("Export destination paths must be distinct.")
  publish_files(staged,destinations)
  if(return_reads) return(Reduce(ShortRead::append, reads))
  invisible(nreads)
}

#' @title extractAnnotFromQname
#' @description Extract the annotation for each read from the read name.
#' @param qname A vector of string which stores the read name.
#' @param annot The prefix of the annotation to be extracted
#' @importFrom stringr str_extract
#' @return A vector of string which stores the annotations.
extractAnnotFromQname = function(qname,annot = "cell"){
  regrex = paste(c("(?<=",annot,"=)[^|]+"),collapse = "")
  out = str_extract(qname, regrex)
  return(out)
}
