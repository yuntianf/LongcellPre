#' @title splice_site_table
#'
#' @description extract splice sites from isoform sequence.
#' @details extract splice sites from isoform sequence and filter out artifacts.
#'
#' @param isoform The string which represents the isoform
#' @param split The character to seperate different exons in the isoform representation.
#' @param sep The character to seperate the start and end position of an exon in the isoform representation.
#' @param splice_site_thresh The minimum number of occurance for a splice site to be preserved.
#' @importFrom Longcellsrc splice_site_table_cpp
#' @return A dataframe, the first and the last column stores the start and end position of each read,
#' the middle columns store if each read has the splicing sites
splice_site_table <- function(isoforms,
                              split = "|", sep = ",",
                              splice_site_thresh = 10) {
  out <- splice_site_table_cpp(isoforms, split, sep, splice_site_thresh)

  # Ensure numeric & data.frame once
  out$start <- as.numeric(out$start)
  out$end   <- as.numeric(out$end)
  out = as.data.frame(do.call(cbind,out))

  # return(out)
  if (ncol(out) > 3) {
    # Identify, numeric-order, and keep mid-site columns
    mid_names <- setdiff(colnames(out), c("id", "start", "end"))
    ord <- order(suppressWarnings(as.numeric(mid_names)))
    mid_cols <- mid_names[ord]

    # Row-wise bounds computed once
    end_thresh = 0
    lower <- out$start + end_thresh
    upper <- out$end   - end_thresh

    # Update each column in place to avoid big temporary matrices
    for (nm in mid_cols) {
      v <- suppressWarnings(as.numeric(out[[nm]]))  # numeric, NAs preserved
      # One logical vector per column; no large matrices
      bad <- (as.numeric(nm) < lower) | (as.numeric(nm) > upper)
      # Fast path: if any finite values exist, mask out-of-bounds
      if (any(bad, na.rm = TRUE)) v[bad] <- NA_real_
      out[[nm]] <- v
    }

    # Reassemble columns in desired order without cbind copying chains
    out <- out[c("id", "start", mid_cols, "end")]
  }

  rownames(out) <- NULL
  out
}

#' @title mid_len
#'
#' @description count the number of middle splicing sites covered by a read
#'
#' @param mid A vector of strings, each string denotes the existence of the middle splicing sites
#' @param sep The character to seperate splicing sites in the string
#' @return A dataframe with two columns, the first column denotes the middle splicing sites,
#' the second denotes the number of sites which covered.
mid_len = function(mid,sep = ","){
  mat = strsplit(mid,split = sep,fixed = TRUE)
  mat = do.call(rbind,mat)
  suppressWarnings(storage.mode(mat) <- "numeric")

  out = as.data.frame(cbind(mid,rowSums(!is.na(mat))))
  colnames(out) = c("mid","size")
  out = out %>% mutate(size = as.numeric(size))
  return(out)
}

#' @title mid_group
#'
#' @description group middle splicing sites which can transform to each other by truncations.
#'
#' @inheritParams mid_len
#' @importFrom Matrix Matrix
#' @importFrom Longcellsrc matrix_xor
#' @return A dataframe with two columns, the first column denotes the middle splicing site string,
#' the second denotes the splicing site string that can transform to this splicing site string by
#' truncations.
#'
mid_group = function (mid, sep = ",")
{
  mat = strsplit(mid, split = sep, fixed = TRUE)
  mat = do.call(rbind, mat)
  suppressWarnings(storage.mode(mat) <- "numeric")
  mid = mid[order(rowSums(is.na(mat)), decreasing = TRUE)]
  mat = mat[order(rowSums(is.na(mat)), decreasing = TRUE),
            , drop = FALSE]
  if (nrow(mat) == 1) {
    result = as.data.frame(cbind(mid, mid))
    colnames(result) = c("c", "p")
    return(result)
  }
  mask = as.matrix(is.na(mat))

  NA_flag = (mask %*% t(mask) - rowSums(mask)) == 0
  ones_flag = matrix_xor(mat)
  result = Matrix(NA_flag & ones_flag, sparse = TRUE)
  result = as.data.frame(summary(as(result, "generalMatrix")))
  if (length(result) > 0) {
    result = result %>% filter(i != j) %>% dplyr::select(-x)
    result = result[, c("j", "i")]
    colnames(result) = c("c", "p")
    result = result %>% mutate(c = mid[c], p = mid[p])
    orphan = setdiff(mid, result$c)
    if (length(orphan) > 0) {
      orphan = as.data.frame(cbind(orphan, NA))
      colnames(orphan) = c("c", "p")
      result = rbind(result, orphan)
    }
  }
  else {
    result = as.data.frame(cbind(mid, NA))
    colnames(result) = c("c", "p")
  }
  return(result)
}



#' @title disagree_sites
#'
#' @description Identify the different sites between different isoforms within the same
#' UMI cluster
#'
#' @param from,to A string vector of middle splicing sites. From is the original middle splicing
#' sites, to is the middle splicing sites the from mapped to within the same UMI cluster
#' @return A dataframe, recording the number of splicing sites existing in the correctly mapped
#' and wrongly mapped reads.
disagree_sites = function(from,to){
  if(length(from) != length(to)){
    stop("There should be a one to one correspondence!")
  }
  from_table = do.call(rbind,strsplit(from,split = ","))
  to_table = do.call(rbind,strsplit(to,split = ","))

  from_table <- suppressWarnings(matrix(as.numeric(from_table),ncol = ncol(from_table)))
  to_table <- suppressWarnings(matrix(as.numeric(to_table),ncol = ncol(to_table)))
  #return(list(from_table,to_table))
  disagree = xor(from_table,to_table) & (from_table == 1)

  filter = as.data.frame(cbind(colSums(disagree,na.rm = TRUE),
                               colSums(from_table,na.rm = TRUE),
                               colSums(to_table,na.rm = TRUE)))
  colnames(filter) = c("disagree","wrong","correct")

  return(filter)
}

#' @title site_correct
#'
#' @description correct wrongly mapped splicing sites.
#'
#' @inheritParams disagree_sites
#' @param sep The character to seperate splicing sites in the middle splicing site string.
#' @return A datafrom recording the mapping between wrongly mapped isoform and correct isoform.
site_correct = function(from,to,sep = ","){
  if(length(from) != length(to)){
    stop("The length of from and to should match!")
  }
  out = sapply(1:length(from),function(i){
    if(is.na(to[i])){
      return(to[i])
    }
    if(from[i] == to[i]){
      return(to[i])
    }
    else{
      from_vec = suppressWarnings(as.numeric(unlist(strsplit(from[i],split = sep,fixed = TRUE))))
      to_vec = suppressWarnings(as.numeric(unlist(strsplit(to[i],split = sep,fixed = TRUE))))

      from_vec[!is.na(from_vec) & !is.na(to_vec)] = to_vec[!is.na(from_vec) & !is.na(to_vec)]
      return(paste(from_vec,collapse = sep))
    }
  })
  return(out)

}

#' @title cells_mid_filter
#'
#' @description Filter out middle splice sites pattern caused by wrong mapping.
#' @details Filter out middle splice sites pattern caused by wrong mapping given the information from UMI cluster as
#' reads within the same UMI cluster should come from the same transcript.
#'
#' @param cells The string vector indicating the cell
#' @param cluster The factor vector indicating the UMI cluster
#' @param isoform_mid The string vector indicating middle splicing sites.
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr reframe
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr select
#' @importFrom tidyr pivot_wider
#' @return A string vector indicating preserved patterns of middle splicing sites

mid_correct_input = function(cells,cluster,gene_isoform){
  if(nrow(gene_isoform) != length(cells)){
    stop("The size of isoforms and cells don't match!")
  }
  if(length(cluster) != length(cells)){
    stop("The size of clusters and cells don't match!")
  }
  mid = apply(gene_isoform[,2:(ncol(gene_isoform)-1)],
              1,function(x) paste(x,collapse = ","))
  out = as.data.frame(cbind(cells,cluster,
                            gene_isoform$start,mid,gene_isoform$end))
  colnames(out) = c("cell","cluster","start","mid","end")
  return(out)
}

#' @title mid_coexist_fast
#' Collapse per-(cell, cluster) mids into a consensus and coexistence map (fast)
#'
#' @description
#' Given per-cell assignments of clone/lineage IDs (`mid`) within clusters, this
#' computes (1) a **consensus mid** for each `(cell, cluster)` and
#' (2) an aggregated **coexistence map** counting how often a source mid
#' maps to a (possibly collapsed) target mid. Ties are broken by within-group
#' count, then by global mid size.
#'
#' @param data data.frame or data.table with at least these columns:
#' \itemize{
#'   \item \code{cell}: cell/barcode identifier
#'   \item \code{cluster}: cluster/partition label for the cell
#'   \item \code{mid}: clone/lineage identifier
#' }
#'
#' @details
#' This function follows a two-step strategy per \code{(cell, cluster)}:
#'
#' \enumerate{
#'   \item Compute counts per \code{mid}; rank mids by descending count, then by
#'         a global \emph{mid size} (from \code{mid_len()}); take the top 1–2.
#'   \item If two mids are present and exactly one is a parent of the other
#'         (according to \code{mid_group()} edges \code{c -> p}), collapse both to
#'         the dominant direction's \emph{mode mid}; otherwise, map each observed
#'         mid to the group's top-1 \emph{mode mid}.
#' }
#'
#' Expectations for helpers:
#' \itemize{
#'   \item \code{mid_len(total)} returns a two-column object with columns
#'         \code{mid}, \code{size}.
#'   \item \code{mid_group(total)} returns a two-column object with columns
#'         \code{c} (child), \code{p} (parent). It may be empty; in that case,
#'         no directed collapsing occurs.
#' }
#'
#' Notes:
#' \itemize{
#'   \item If a \code{mid} has \code{size = NA}, it is treated as \code{-Inf} for ordering.
#'   \item The consensus column is named \code{concensus} (intentional spelling to
#'         match downstream code).
#' }
#'
#' @return A list of two \code{data.table}s:
#' \describe{
#'   \item{\code{concensus}}{Columns: \code{cell}, \code{cluster}, \code{concensus} (the mode mid).}
#'   \item{\code{coexist}}{Columns: \code{from}, \code{to}, \code{count}; aggregated counts of
#'         mappings after the collapse/mapping rules above.}
#' }
#'
#' @import data.table
#'
#' @examples
#' \dontrun{
#' # data: data.frame/data.table with columns cell, cluster, mid
#' res <- mid_coexist_fast(data)
#' concensus_dt <- res[[1]]
#' coexist_dt   <- res[[2]]
#' }
#'
#' @export

mid_coexist_fast <- function(data) {
  # Precompute once (same as your pipeline)
  total  <- names(sort(table(data$mid), decreasing = TRUE))
  len    <- mid_len(total)        # expects cols: mid, size
  parent <- mid_group(total)      # expects cols: c (child), p (parent)

  DT  <- as.data.table(data)
  LEN <- as.data.table(len);    setnames(LEN, c("mid","size"))
  P   <- as.data.table(parent); setnames(P,   c("c","p"))
  if (nrow(P)) setkey(P, c, p)

  # 1) counts per (cell, cluster, mid)
  CNT <- DT[, .N, by = .(cell, cluster, mid)]
  CNT <- LEN[CNT, on = "mid"]
  CNT[is.na(size), size := -Inf]                 # safe ordering if size NA
  setorder(CNT, cell, cluster, -N, -size)

  # 2) top-2 rows per group
  TOP2 <- CNT[, head(.SD, 2L), by = .(cell, cluster)]

  # split groups with 1 vs 2 mids
  NG      <- TOP2[, .N, by = .(cell, cluster)]
  groups1 <- TOP2[NG[N == 1L], on = .(cell, cluster)]
  groups2 <- TOP2[NG[N == 2L], on = .(cell, cluster)]

  # 3) groups with one mid: (mid -> mid, count)
  res_one  <- groups1[, .(from = mid, to = mid, count = N)]
  cons_one <- groups1[, .(cell, cluster, concensus = mid)]

  # 4) groups with two mids: check parent relation both directions
  groups2[, other_mid := mid[.N:1L], by = .(cell, cluster)]
  if (nrow(P)) {
    groups2[, is_child_of_other := P[.SD, on = .(c = mid, p = other_mid), .N, by = .EACHI]$N > 0]
  } else {
    groups2[, is_child_of_other := FALSE]
  }

  # >>> FIX: compute locals inside {...} and use them <<<
  CONS2 <- groups2[, {
    hd <- sum(is_child_of_other) == 1L
    mm <- if (hd) other_mid[is_child_of_other][1L] else mid[1L]
    sc <- sum(N)
    .(has_dir = hd, mode_mid = mm, sum_count = sc)
  }, by = .(cell, cluster)]

  # collapse if directed; otherwise map both to top1
  res_two_collapse <- CONS2[has_dir == TRUE,
                            .(from = mode_mid, to = mode_mid, count = sum_count)]
  res_two_map <- groups2[CONS2[has_dir == FALSE], on = .(cell, cluster)][
    , .(from = mid, to = mode_mid, count = N)]

  coexist <- rbindlist(list(res_one, res_two_collapse, res_two_map), use.names = TRUE, fill = TRUE)
  coexist <- coexist[, .(count = sum(as.numeric(count))), by = .(from, to)]

  concensus <- rbindlist(
    list(cons_one, CONS2[, .(cell, cluster, concensus = mode_mid)]),
    use.names = TRUE, fill = TRUE
  )

  list(concensus, coexist)
}

#' @title isoform_corres
#' Infer isoform correspondence from a coexistence count matrix
#'
#' @description
#' Builds a row-wise mapping from each isoform (row) to its most likely
#' corresponding isoform(s) (column) using a square coexistence matrix of
#' pairwise counts, with optional dispute resolution to drop unreliable
#' target columns before finalizing the map.
#'
#' @details
#' Workflow:
#' \enumerate{
#'   \item Identify a set of \emph{trustworthy targets}: keep columns whose
#'         diagonal count is strictly greater than twice their total off-diagonal
#'         mass, i.e. \eqn{M_{ii} > 2 \sum_{j \neq i} M_{ji}}.
#'   \item Initial mapping: for each row, pick the target column among the
#'         trustworthy set with the largest count (ties broken by first match);
#'         rows with all-zero counts map to \code{NA}.
#'   \item Dispute resolution: if any row maps to a different column than itself,
#'         call \code{disagree_sites(from, to)} and, for columns where
#'         \code{wrong > correct} and \code{disagree * 3 > wrong}, drop those
#'         columns from the trustworthy set and recompute the mapping restricted
#'         to the remaining targets.
#' }
#'
#' @param coexist_matrix A square numeric matrix of coexistence counts with
#'   matching row and column names (isoform identifiers). Diagonal entries
#'   represent self-coexistence.
#'
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#'
#' @return A two-column \code{data.frame} with row names set to \code{from}:
#'   \describe{
#'     \item{\code{from}}{Row isoform identifier.}
#'     \item{\code{to}}{Mapped target isoform identifier (or \code{NA} if no evidence).}
#'   }
#'
isoform_corres = function(coexist_matrix){
  isoform_coexist_filter <- cbind(diag(coexist_matrix),
                                  rowSums(coexist_matrix)-diag(coexist_matrix))
  isoform_coexist_filter <- isoform_coexist_filter[order(isoform_coexist_filter[,1],
                                                         decreasing = T),]
  isoform_coexist_filter <- names(which(isoform_coexist_filter[,1] >
                                          2*isoform_coexist_filter[,2]))

  corres = as.data.frame(coexist_matrix[,isoform_coexist_filter])
  corres = sapply(1:nrow(corres),function(i){
    x = unlist(corres[i,])
    if(sum(x) == 0){
      return(NA)
    }
    return(isoform_coexist_filter[which(x == max(x))[1]])
  })

  corres = as.data.frame(cbind(rownames(coexist_matrix),corres))
  colnames(corres) = c("from","to")
  rownames(corres) = corres$from

  disagree = corres %>% filter(from != to,!is.na(to))
  if(nrow(disagree) == 0){
    return(corres)
  }
  ds = disagree_sites(disagree$from,disagree$to)
  id = which(ds$wrong > ds$correct & ds$disagree*3 > ds$wrong)

  if(length(id) == 0){
    return(corres)
  }

  correct_iso = na.omit(unique(corres$to))
  to_table = do.call(rbind,strsplit(correct_iso,split = ","))
  to_table <- suppressWarnings(matrix(as.numeric(to_table),ncol = ncol(to_table)))
  to_table[is.na(to_table)] = 0
  if(length(id) == 1){
    correct_iso = correct_iso[to_table[,id] == 0]
  }
  else if(nrow(to_table) == 1){
    correct_iso = correct_iso[sum(to_table[,id]) == 0]
  }
  else{
    correct_iso = correct_iso[rowSums(to_table[,id]) == 0]
  }

  corres_new = as.data.frame(coexist_matrix[,correct_iso])
  corres_new = sapply(1:nrow(corres),function(i){
    x = unlist(corres_new[i,])
    if(sum(x) == 0){
      return(NA)
    }
    return(correct_iso[which(x == max(x))[1]])
  })

  corres_new = as.data.frame(cbind(rownames(coexist_matrix),corres_new))
  colnames(corres_new) = c("from","to")
  rownames(corres_new) = corres_new$from
  return(corres_new)
}

#' @title cluster_isoform_correct
#'
#' @description Corret for wrong mapping and truncation within each UMI cluster
#' @details Corret for wrong mapping and truncation within each UMI cluster
#'
#' @param start A numeric vector indicating start position of reads within the UMI cluster
#' @param mid A string vector indicating patterns of middle splicing sites
#' @param end A numeric vector indicating end position of reads within the UMI cluster
#' @param polyA A boolean vector indicating the polyA existence
#' @param preserve_mid A string vector indicating preserved patterns of middle splice sites
#' @return A vector including the corrected start, middle splicing sites, end position, polyA existence
#' and the number of reads within this UMI cluster
cluster_isoform_correct <- function(start,mid,end,concensus,polyA,preserve_mid){
    if(length(start) != length(mid) || length(start) != length(end)){
      stop("The size of isoforms representation don't match!")
    }
    if(length(start) != length(polyA)){
        stop("The size of isoforms and polyA don't match!")
    }
    cluster_size = length(start)
    mode_isoform = NA
    mode_start = NA
    mode_end = NA
    mode_polyA = NA

    if(length(preserve_mid) > 0 && nrow(preserve_mid) > 0){
      mid_corres = na.omit(preserve_mid[concensus,"to"])
      if(length(mid_corres) > 0){
        mode_isoform = mid_corres
        mode_preserve = which(mid == mode_isoform)
        if(length(mode_preserve) == 0){
          mode_start = -1
          mode_end = -1
        }
        else{
          mode_start = min(start[mode_preserve])
          #mode_start = table(start[mode_preserve])
          #mode_start = min(as.numeric(names(mode_start[which(mode_start == max(mode_start))])))

          mode_end = max(end[mode_preserve])
          #mode_end = table(end[start == mode_start])
          #mode_end = max(as.numeric(names(mode_end[which(mode_end == max(mode_end))])))
        }
        polyA = as.logical(polyA)
        mode_polyA = mean(as.numeric(polyA[mid == concensus]))
      }
    }

    return(c(mode_start,mode_isoform,mode_end,cluster_size,mode_polyA))
}

#' @title isoform_correct
#'
#' @description Corret for wrong mapping and truncation for each UMI cluster in each gene in each cell
#' @details Corret for wrong mapping and truncation for each UMI cluster in each gene in each cell
#'
#' @param gene_isoform A data frame containing cell, gene,and isoform informaion including
#' start, middle splicing sites, end position and polyA existence
#' @param preserve_mid A string vector indicating preserved patterns of middle splice sites
#' @importFrom tidyr replace_na
#' @importFrom dplyr mutate_at
#' @importFrom dplyr group_by
#' @return A dataframe including the corrected start, middle splicing sites, end position, polyA existence
#' and the number of reads of each UMI cluster for each gene in each cell
isoform_correct <- function(gene_isoform,preserve_mid){
  gene_isoform_adjust = gene_isoform %>% group_by(cell,cluster) %>%
                        summarise(adjust = list(cluster_isoform_correct(start,mid,end,unique(concensus),
                                                                        polyA,preserve_mid)),.groups = "drop")
  adjust = as.data.frame(do.call(rbind,gene_isoform_adjust$adjust))
  colnames(adjust) = c("start","mid","end","size","polyA")

  adjust = adjust %>% mutate_at(c("start","end","size","polyA"),as.numeric)
  #adjust = adjust %>% group_by(mid) %>%
  #  mutate_at(c("start","end"),~tidyr::replace_na(., median(., na.rm=TRUE)))

  out = as.data.frame(cbind(gene_isoform_adjust[,c("cell")],adjust))
  return(out[!is.na(out$start),])
}

#' @title cells_mid_correct
#' Correct per-cell isoform assignments using mid-level consensus and coexistence
#'
#' @description
#' Constructs per-(cell, cluster) inputs, computes consensus mids and their
#' coexistence using \code{mid_coexist_fast()}, derives an isoform correspondence
#' map from the coexistence matrix, and applies isoform-level corrections.
#'
#' @details
#' Steps:
#' \enumerate{
#'   \item Build input with \code{mid_correct_input(cells, cluster, gene_isoform)}.
#'   \item Run \code{mid_coexist_fast()} to obtain a consensus mid per group and
#'         a long table of mid-to-mid coexistence counts.
#'   \item If multiple mid pairs exist, convert to a square matrix with
#'         \code{long2square(..., symmetric = FALSE)}, infer the correspondence
#'         via \code{isoform_corres()}.
#'   \item Join consensus to the input, attach \code{polyA}, and apply
#'         \code{isoform_correct(data, corres)} to produce corrected isoform calls.
#' }
#'
#' @param cells A vector of cell/barcode identifiers.
#' @param cluster A vector of cluster labels aligned with \code{cells}.
#' @param gene_isoform Isoform annotation aligned with \code{cells} (format as
#'   required by \code{mid_correct_input()}).
#' @param polyA A numeric (or logical) vector aligned with \code{cells} giving
#'   poly(A) evidence to carry through to the output.
#'
#' @importFrom dplyr left_join select
#'
#' @return A \code{data.frame} of per-(cell, cluster) corrected isoform data
#'   after applying the consensus/coexistence-derived correspondence map.
#'
cells_mid_correct <- function(cells,cluster,gene_isoform,polyA){
  cat("Start to build the middle sites input\n")
  data = mid_correct_input(cells,cluster,gene_isoform)

  cat("Start to build coexistence table for middle sites\n")
  out = mid_coexist_fast(data)
  cat("Finish the coexistence table for middle sites\n")
  concensus = as.data.frame(out[[1]])
  coexist = as.data.frame(out[[2]])

  if(nrow(coexist) == 1){
    corres = as.data.frame(coexist %>% dplyr::select(-count))
    rownames(corres) = corres$from
  }
  # return(coexist)
  else{
    coexist_matrix = long2square(coexist,"from","to","count",symmetric = FALSE)
    corres = isoform_corres(coexist_matrix)
  }

  data = left_join(data,concensus,by = c("cell","cluster"))
  data$polyA = polyA
  # return(data)
  cat("Start to correct isoforms\n")
  data = isoform_correct(data,corres)
  return(data)
}

#' @title cells_nomid_correct
#' Summarize per-(cell, cluster) intervals without mid information
#'
#' @description
#' Produces a compact per-(cell, cluster) summary when mid assignments are
#' unavailable: the minimum \code{start}, maximum \code{end}, the number of
#' contributing records (\code{size}), and the mean \code{polyA}.
#'
#' @details
#' The inputs are combined and then grouped by \code{cell, cluster}. For each
#' group, it returns: \code{start = min(start)}, \code{end = max(end)},
#' \code{size = n()}, and \code{polyA = mean(polyA)}.
#'
#' @param cells A vector of cell/barcode identifiers.
#' @param cluster A vector of cluster labels aligned with \code{cells}.
#' @param gene_isoform A two-column object (matrix/data.frame) giving
#'   \code{start} and \code{end} coordinates aligned with \code{cells}.
#' @param polyA A numeric (or logical) vector aligned with \code{cells} used to
#'   compute per-group mean poly(A) support.
#'
#' @importFrom dplyr group_by summarise
#' @importFrom magrittr %>%
#'
#' @return A \code{data.frame} with columns \code{cell}, \code{cluster},
#'   \code{start}, \code{end}, \code{size}, and \code{polyA}.
#'
cells_nomid_correct <- function(cells,cluster,gene_isoform,polyA){
  data = as.data.frame(cbind(cells,cluster,gene_isoform,polyA))
  colnames(data) = c("cell","cluster","start","end","polyA")

  data = data %>% group_by(cell,cluster) %>%
    summarise(start = min(start),end = max(end),
              size = n(),
              polyA = mean(polyA),.groups = "drop")
  return(data)
}

#' @title site_recover
#'
#' @description transform the middle splicing sites pattern from boolean string into original sites
#' @details transform the middle splicing sites pattern from boolean string into original sites
#'
#' @param start A numeric value indicating start position of reads within the UMI cluster
#' @param mid A string indicating patterns of middle splicing sites
#' @param end A numeric value indicating end position of reads within the UMI cluster
#' @param sites A string vector indicating middle splicing sites
#' @param sep The character to seperate the start and end for each exon bin in the isoform representation
#' @param split The character to seperate each exon bin in the isoform representation
#' @return A string representing isoform
site_recover <- function(start,end,mid = NULL,sites = NULL,flank = 5,sep = ",",split = "|"){
    if(is.null(sites)){
        splice_sites = paste(start,end,sep = sep)
    }
    else{
      mid = suppressWarnings(as.logical(as.numeric(unlist(strsplit(mid,split = ",")))))
      mid[is.na(mid)] = FALSE
      if(length(mid) != length(sites)){
        stop("The size of splicing sites and binary indicator don't match!")
      }
      mid_sites = sites[mid]
      if(as.numeric(start) == -1){
        start = as.numeric(mid_sites[1])-flank
      }
      if(as.numeric(end) == -1){
        end = as.numeric(mid_sites[length(mid_sites)])+flank
      }
      if(is.na(start)){
        return(NA)
      }
      splice_sites = c(start,sites[mid],end)
      splice_sites = t(matrix(splice_sites,nrow = 2))
      splice_sites = apply(splice_sites,1,function(x) paste(x,collapse = sep))
      splice_sites = paste(splice_sites,collapse = split)
    }
    splice_sites = gsub(" ","",splice_sites,fixed = TRUE)
    return(splice_sites)
}

#' @title cells_build_isoform_dt
#' Build isoform strings per (row) with exact site-pairing semantics (fast, data.table)
#'
#' @description
#' Adds an `isoform` column to a table of `{start, end, mid}` using splice-site
#' structure when available. The implementation matches the original
#' `site_recover()` semantics while avoiding `apply(..., 1, ...)` and repeated
#' parsing for speed.
#'
#' @details
#' Behavior by case:
#' \itemize{
#'   \item If \code{sites} is \code{NULL}, each row's isoform is simply
#'         \code{"start,end"} (after trimming spaces).
#'   \item If \code{sites} is provided, \code{mid} is parsed as
#'         \code{as.logical(as.numeric(strsplit(mid, \",\")))}, so
#'         \code{0 -> FALSE}, non-zero \code{-> TRUE}, and \code{NA -> NA}.
#'         For \strong{boundary filling}, indices with non-NA \code{mid} are used
#'         as a mask: if \code{start == -1}, set \code{start = min(sites[mask]) - flank};
#'         if \code{end == -1}, set \code{end = max(sites[mask]) + flank}.
#'         If \code{start} is \code{NA} in this mode, the result is \code{NA}.
#'         For \strong{interior selection}, \code{NA} entries in \code{mid} are set to
#'         \code{FALSE} and only \code{TRUE} positions are kept: the sequence
#'         \code{c(start, sites[TRUE], end)} is then paired \emph{in chunks of two}
#'         (i.e., \code{(1,2)}, \code{(3,4)}, \code{(5,6)}, ...) to form
#'         \code{"a,b|c,d|..."} using \code{sep} between numbers and \code{split}
#'         between pairs. A length mismatch between parsed \code{mid} and
#'         \code{sites} raises an error.
#' }
#'
#' @param data A \code{data.frame} or \code{data.table} with columns:
#'   \code{start}, \code{end}, and (when \code{sites} is not \code{NULL}) \code{mid}.
#' @param sites \code{NULL} or a vector of splice-site coordinates (numeric or
#'   coercible to numeric) whose length must match the number of entries parsed
#'   from each \code{mid}.
#' @param flank Integer; number of bases to extend when filling \code{start == -1}
#'   and/or \code{end == -1}. Default \code{5}.
#' @param sep String used between endpoints within a pair (default \code{","}).
#' @param split String used between consecutive pairs (default \code{"|"}).
#'
#' @importFrom data.table as.data.table
#'
#' @return The input as a \code{data.table} with an added \code{isoform} column.
#'   Assignment is done by reference when possible; the function returns
#'   \code{dt[]} for convenience.
cells_build_isoform_dt <- function(data, sites = NULL, flank = 5L,
                                   sep = ",", split = "|") {
  dt <- as.data.table(data)

  if (is.null(sites)) {
    s_num <- suppressWarnings(as.numeric(dt$start))
    e_num <- suppressWarnings(as.numeric(dt$end))
    dt[, isoform := gsub(" ", "", paste(s_num, e_num, sep = sep), fixed = TRUE)]
    return(dt[])
  }

  sites_num <- suppressWarnings(as.numeric(sites))
  if (anyNA(sites_num)) stop("`sites` must be coercible to numeric.")

  s_num <- suppressWarnings(as.numeric(dt$start))
  e_num <- suppressWarnings(as.numeric(dt$end))

  mid_chr <- as.character(dt$mid); mid_chr[is.na(mid_chr)] <- ""
  mids_list <- strsplit(mid_chr, ",", fixed = TRUE)

  mids_num_list <- lapply(mids_list, function(x)
    if (length(x)) suppressWarnings(as.numeric(x)) else numeric(0L))
  mids_log_list <- lapply(mids_num_list, function(x) suppressWarnings(as.logical(x)))

  if (any(vapply(mids_log_list, length, 1L) != length(sites_num))) {
    stop("The size of splicing sites and binary indicator don't match!")
  }

  # Boundary mask: non-NA
  masks <- lapply(mids_log_list, function(m) !is.na(m))
  first_idx <- vapply(masks, function(m) if (any(m)) which(m)[1L] else NA_integer_, 1L)
  last_idx  <- vapply(masks, function(m) if (any(m)) tail(which(m), 1L) else NA_integer_, 1L)

  repl_s <- !is.na(s_num) & s_num == -1 & !is.na(first_idx)
  repl_e <- !is.na(e_num) & e_num == -1 & !is.na(last_idx)
  if (any(repl_s)) s_num[repl_s] <- sites_num[first_idx[repl_s]] - flank
  if (any(repl_e)) e_num[repl_e] <- sites_num[last_idx[repl_e]] + flank

  dt[, isoform := vapply(seq_len(.N), function(i) {
    si <- s_num[i]
    if (is.na(si)) return(NA_character_)
    ei <- e_num[i]
    m  <- mids_log_list[[i]]
    if (length(m)) m[is.na(m)] <- FALSE
    seg <- c(si, sites_num[if (length(m)) m else FALSE], ei)

    # Pair in chunks of two (matrix(nrow=2) semantics)
    k <- length(seg) %/% 2L
    if (k == 0L) {
      paste(si, ei, sep = sep)
    } else {
      left  <- seg[2L * seq_len(k) - 1L]
      right <- seg[2L * seq_len(k)]
      paste(paste(left, right, sep = sep), collapse = split)
    }
  }, FUN.VALUE = character(1L))]

  dt[, isoform := gsub(" ", "", isoform, fixed = TRUE)]
  dt[]
}


#' @title cells_isoform_correct
#' Correct and assign isoforms per (cell, cluster)
#'
#' @description
#' Generates corrected isoform calls for each \code{(cell, cluster)} by combining
#' mid-level correction (when splice-site information is available) or by
#' summarizing intervals (when no mid information exists), followed by recovery
#' of isoform identifiers from start/end positions and splice sites.
#'
#' @details
#' Workflow:
#' \enumerate{
#'   \item If \code{gene_isoform} has >2 columns (implying splice-site structure):
#'         \itemize{
#'           \item Extract splice-site names from column 2 to (n-1).
#'           \item Run \code{cells_mid_correct()} to obtain consensus mids and corrected
#'                 per-(cell, cluster) records.
#'           \item For each row, recover isoform identity using
#'                 \code{site_recover(start, end, mid, splice_sites)}.
#'         }
#'   \item Otherwise:
#'         \itemize{
#'           \item Summarize per-group intervals via \code{cells_nomid_correct()}.
#'           \item Recover isoforms using \code{site_recover(start, end)}.
#'           \item Assign a placeholder mid of \code{"null"}.
#'         }
#'   \item Drop empty results and remove \code{start, end} columns after recovery.
#' }
#'
#' @param cells A vector of cell/barcode identifiers.
#' @param cluster A vector of cluster labels aligned with \code{cells}.
#' @param gene_isoform A matrix/data.frame giving isoform structure. If it has
#'   more than two columns, internal splice-site columns are used for isoform
#'   recovery.
#' @param polyA A numeric (or logical) vector aligned with \code{cells} used to
#'   track poly(A) evidence.
#'
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#'
#' @return A \code{data.frame} of corrected isoform calls per (cell, cluster),
#'   including:
#'   \describe{
#'     \item{\code{cell}}{Cell identifier.}
#'     \item{\code{cluster}}{Cluster label.}
#'     \item{\code{mid}}{Consensus mid (or \code{"null"} if unavailable).}
#'     \item{\code{size, polyA}}{Group size and mean poly(A) support.}
#'     \item{\code{isoform}}{Recovered isoform identifier.}
#'   }
#'
cells_isoform_correct <- function(cells,cluster,gene_isoform,polyA){
  if(ncol(gene_isoform) > 2){
    cat("Start to prepare input\n")
    splice_sites = colnames(gene_isoform)[2:(ncol(gene_isoform)-1)]
    cat("Start to correct the middle sites\n")
    data = cells_mid_correct(cells,cluster,gene_isoform,polyA)
    if(nrow(data) == 0){
      return(NULL)
    }
    #data$isoform <- apply(data,1,function(x){
    #  site_recover(x["start"],x["end"],x["mid"],splice_sites)
    #})
    cat("Start to build isoform table\n")
    data = cells_build_isoform_dt(data,sites = splice_sites, flank = 5L,sep = ",", split = "|")
  }
  else{
    splice_sites = NULL
    data = cells_nomid_correct(cells,cluster,gene_isoform,polyA)
    if(nrow(data) == 0){
      return(NULL)
    }
    #data$isoform <- apply(data,1,function(x){
    #  site_recover(x["start"],x["end"])
    #})
    data = cells_build_isoform_dt(data, flank = 5L,sep = ",", split = "|")
    data$mid = "null"
  }

  data = as.data.frame(data)
  data <- na.omit(data) %>% dplyr::select(-c(start,end))
  return(data)
}
