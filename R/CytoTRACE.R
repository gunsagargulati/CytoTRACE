#' @title cytoTRACE
#'
#' @description This function generates single-cell predictions of differentiation status from scRNA-seq data.
#' It takes in a matrix of gene expression values where columns are cells and rows are genes.\cr\cr
#' (Optional) Technical batches can be corrected with ComBat by inputting a character vector of batch
#' annotations matching the number of columns (cells) in the matrix. Furthermore, for the analysis of
#' large datasets (>5,000 cells), users can increase the speed performance by enabling a sub-sampling approach for calculation
#' and using multiple cores
#'
#' @param mat matrix of gene expression values where columns are cells and rows are genes
#' @param batch character vector of length equal to the number of columns (cells) in the matrix
#' @param enableFast boolean indicating whether or not to run CytoTRACE in 'fast mode'
#' @param ncores integer indicating the number of cores to utilize when enableFast = TRUE
#' @param subsamplesize integer indicating the number of cells to subsample when enableFast = TRUE
#'
#' @return a list containing
#' \itemize{
#' \item CytoTRACE: a numeric vector of the predicted ordering of single cells by differentiation status
#' \item GCS: a numeric vector of the gene counts signature (geometric mean of the top 200 genes associated with gene counts)
#' \item GCSgenes: a numeric vector of the Pearson correlation of each gene with gene counts, ordered from highest to lowest
#' \item Counts: a numeric vector of the number of genes expressed per single cell (gene counts)
#' \item filteredCells = a character vector of the names of single cells (columns) that were filtered due to poor quality.
#' }
#'
#' @author Gunsagar Gulati <cytotrace@gmail.com>
#'
#' @seealso https://cytotrace.stanford.edu
#'
#' @references https://doi.org/10.1101/649848
#'
#' @examples
#' #Use the bone marrow 10x scRNA-seq dataset to run CytoTRACE
#' results <- cytoTRACE(marrow_10x_expr)
#'
#' #Run this dataset on fast mode using 8 cores
#' results <- cytoTRACE(marrow_10x_expr, enableFast = TRUE, ncores = 8)
#'
#' @export

cytoTRACE <- function(mat, batch = NULL, enableFast = FALSE,
                      ncores = 1,subsamplesize = 1000){

  #inputs
  a1 <- mat
  a2 <- batch

  #Checkpoint: NAs and poor quality genes
  pqgenes <- is.na(rowSums(mat>0)) | apply(mat, 1, var) == 0
  num_pqgenes <- length(which(pqgenes == TRUE))
  mat <- mat[!pqgenes,]
  if(num_pqgenes>0){
    warning(paste(num_pqgenes, "genes have zero expression in the matrix and were filtered"))
  }

  #Subsample routine
  if(enableFast == FALSE){
    size <- ncol(mat)
  } else if (enableFast == TRUE & subsamplesize < ncol(mat)){
    size <- subsamplesize
    message("CytoTRACE will be run on fast mode.")
  } else if (enableFast == TRUE & subsamplesize >= ncol(mat)){
    stop("Please choose a subsample size less than the number of cells in dataset.")
  }

  chunk <- round(ncol(mat)/size)
  subsamples <- split(1:ncol(mat), sample(factor(1:ncol(mat) %% chunk)))
  message(paste("CytoTRACE will be run on", chunk, "sub-sample(s) of approximately", round(mean(unlist(lapply(subsamples, length)))), "cells each using", min(chunk, ncores), "core(s)"))

  message(paste("Pre-processing data and generating similarity matrix..."))
  batches <- parallel::mclapply(subsamples, mc.cores = min(chunk, ncores), function(subsample){
    #Checkpoint: log2-normalization
    mat <- mat[,subsample]
    batch <- batch[subsample]

    if(max(mat)<50){
      mat <- 2^mat - 1
    }

    #Checkpoint: ERCC standards
    if(length(grep("ERCC-", rownames(mat)))>0){
      mat <- mat[-grep("ERCC-", rownames(mat)),]
    }

    #Checkpoint: Sequencing depth normalization
    mat <- t(t(mat)/apply(mat, 2, sum))*1000000

    #Checkpoint: NAs and poor quality cells
    pqcells <- is.na(apply(mat>0, 2, sum)) | apply(mat>0, 2, sum) <= 10
    num_pqcells <- length(which(pqcells == TRUE))
    mat <- mat[,!pqcells]

    #Checkpoint: log2-normalize
    mat <- log(mat+1,2)
    mat <- data.matrix(mat)

    #Calculate pre-batch corrected gene counts
    counts <- apply(mat>0, 2, sum)
    #Checkpoint: Batch correction
    if(ncol(a1) == length(a2)){
      #filter poor quality cells from batch vector
      batch <- batch[!pqcells]

      #Run Combat
      suppressMessages(mat <- sva::ComBat(mat, batch, c()))
      mat <- data.matrix(mat)

      #Replace negative values after batch correction with zeroes for compatibility with downstream steps
      mat[which(mat<0)] <- 0
    }
    #Rescale each single cell with gene counts to convert relative transcript abundances to absolute RNA content prior to cell lysis (credit: Census, Qiu et al., 2017)
    census_normalize <- function(mat, counts) {
      xnl <- 2^data.matrix(mat) - 1
      rs <- apply(xnl, 2, sum)
      rnorm <- t(t(xnl) * counts/rs)
      A <- log(rnorm+1,2)
      return(A)
    }

    mat2 <- census_normalize(mat, counts)
    #Function to identify the most variable genes
    mvg <- function(matn) {
      A <- matn
      n_expr <- rowSums(A > 0);
      A_filt <- A[n_expr >= 0.05 * ncol(A),];
      vars <- apply(A_filt, 1, var);
      means <- apply(A_filt, 1, mean);
      disp <- vars / means;
      last_disp <- tail(sort(disp), 1000)[1];
      A_filt <- A_filt[disp >= last_disp,];

      return(A_filt)
    }

    #Filter out cells not expressing any of the 1000 most variable genes
    mat2.mvg <- mvg(mat2)
    rm1 <- colSums(mat2.mvg) == 0
    mat2 <- mat2[, !rm1]
    counts <- counts[!rm1]

    #Calculate similarity matrix
    similarity_matrix_cleaned <- function(similarity_matrix){
      D <- similarity_matrix
      cutoff <- mean(as.vector(D))
      diag(D) <- 0;
      D[which(D < 0)] <- 0;
      D[which(D <= cutoff)] <- 0;
      Ds <- D
      D <- D / rowSums(D);
      D[which(rowSums(Ds)==0),] <- 0
      return(D)
    }
    D <- similarity_matrix_cleaned(HiClimR::fastCor(mvg(mat2)))

    return(list(mat2 = mat2,counts = counts, D = D))
  }
  )
  #Prepare for downstream steps
  mat2 <- do.call(cbind, lapply(batches, function(x) x$mat2))
  counts <- do.call(c, lapply(batches, function(x) x$counts))
  filter <- colnames(a1)[-which(colnames(a1) %in% colnames(mat2))]
  if(length(filter)>0){
    warning(paste(length(filter), "poor quality cells were filtered based on low or no expression. See 'filteredCells' in returned object for names of filtered cells."))
  }
  #Calculate gene counts signature (GCS) or the genes most correlated with gene counts
  message("Calculating gene counts signature...")
  ds2 <- sapply(1:nrow(mat2), function(x) ccaPP::corPearson(mat2[x,],counts))
  names(ds2) <- rownames(mat2)
  gcs <- apply(mat2[which(rownames(mat2) %in% names(rev(sort(ds2))[1:200])),],2,mean)

  samplesize <- unlist(lapply(lapply(batches, function(x) x$counts), length))
  gcs2 <- split(gcs, as.numeric(rep(names(samplesize), samplesize)))
  D2 <- lapply(batches, function(x) x$D)

  #Regress gene counts signature (GCS) onto similarity matrix
  regressed <- function(similarity_matrix_cleaned, score){
    out <- nnls::nnls(similarity_matrix_cleaned,score)
    score_regressed <- similarity_matrix_cleaned %*% out$x
    return(score_regressed)
  }

  #Apply diffusion to regressed GCS using similarity matrix
  diffused <- function(similarity_matrix_cleaned, score, ALPHA = 0.9){
    vals <- score
    v_prev <- rep(vals);
    v_curr <- rep(vals);

    for(i in 1:10000) {
      v_prev <- rep(v_curr);
      v_curr <- ALPHA * (similarity_matrix_cleaned %*% v_curr) + (1 - ALPHA) * vals;

      diff <- mean(abs(v_curr - v_prev));
      if(diff <= 1e-6) {
        break;
      }
    }
    return(v_curr)
  }

  message("Smoothing values with NNLS regression and diffusion...")
  cytotrace <- parallel::mclapply(1:length(D2), mc.cores = ncores, function(i) {
    gcs_regressed <- regressed(D2[[i]], gcs2[[i]])
    gcs_diffused <- diffused(D2[[i]], gcs_regressed)
    cytotrace <- rank(gcs_diffused)
  }
  )

  cytotrace <- unlist(cytotrace)
  names(cytotrace) <- names(gcs) <- names(counts) <- colnames(mat2)
  message("Done")
  return(list(CytoTRACE = cytotrace, GCS = gcs, GCSgenes = sort(ds2, decreasing = T), Counts = counts, filteredCells = filter))
}







