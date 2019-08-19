#' @title integrated CytoTRACE (iCytoTRACE)
#'
#' @description This function generates single-cell predictions of differentiation status
#' across multiple, heterogeneous scRNA-seq batches/datasets. This implementation leverages
#' mutual nearest neighbor and Gaussian kernel normalization techniques from Scanorama to merge datasets.
#' It takes in a list of gene expression matrices where columns are cells and rows are genes and outputs
#' integrated gene counts, gene counts signature (GCS), CytoTRACE, and the correct gene expression matrix.
#'
#'
#' @param datasets list of gene expression matrices where columns are cells and rows are genes
#'
#' @return a list containing
#' \itemize{
#' \item data: a matrix (genes by single cells) of the corrected gene expression values from Scanorama
#' \item CytoTRACE: a numeric vector of the predicted ordering of merged single cells by differentiation status
#' \item GCS: a numeric vector of the merged gene counts signature (geometric mean of the top 200 genes associated with gene counts)
#' \item Counts: a numeric vector of the corrected number of genes expressed per single cell (gene counts)
#' \item coord: a matrix containing the coordinates for the merged low-dimensional embedding (number of cells x 100 <i>t</i>-SNE components)
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
#'
#' #Create a list containing two bone marrow scRNA-seq datasets profiled on different platforms, 10x and Smart-seq2
#' datasets <- list(marrow_10x_expr, marrow_plate_expr)
#' \n
#' #Run iCytoTRACE
#' \dontrun{
#' results <- iCytoTRACE(datasets)
#' }
#' @export



iCytoTRACE <- function(datasets) {
  if(!have_scanoramaCT | !have_numpy){
    stop("The necessary python modules are not accessible. This function is disabled. Please follow the instructions in https://github.com/gunsagargulati/CytoTRACE to install the Python packages for this application.")
  }
  #import python libraries
  scanoramaCT<- reticulate::import('scanoramaCT')
  np <- reticulate::import('numpy')
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}

  if(class(datasets) != "list" | length(datasets) < 2){
    print("Please input a list of at least 2 gene by cell matrices")
  }

  processedMats <- lapply(datasets, function(mat){
    #Checkpoint: log2-normalization
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

    #Checkpoint: NAs and poor quality genes
    pqgenes <- is.na(rowSums(mat>0)) | apply(mat, 1, var) == 0
    num_pqgenes <- length(which(pqgenes == TRUE))
    mat <- mat[!pqgenes,]

    #Checkpoint: log2-normalize
    mat <- log(mat+1,2)
    mat <- data.matrix(mat)
    return(mat)
  }
  )
  #Calculate and min-max normalize gene counts
  countsOrig <- lapply(processedMats, function(x) colSums(x>0))
  countsNorm <- lapply(countsOrig, range01)
  countsNorm <- lapply(countsNorm, function(x) cbind(x, x))

  #Run modified scanormaCT to batch correct matrix and gene counts
  #Create genes_list
  genes_list <- lapply(processedMats, rownames)
  #transpose matrices
  processedMats <- lapply(processedMats, t)

  merge <- scanormaCT$merge_datasets(processedMats, genes_list)
  process <- scanormaCT$process_data(merge[[1]], merge[[2]], hvg = 0L, dimred = 100L)
  gcalign <- scanormaCT$find_alignments_gc(process[[1]], countsNorm)
  integrated.corrected.data <- scanormaCT$correct(processedMats, genes_list, return_dimred=TRUE, return_dense = T)
  countsNormCorrected <- unlist(lapply(gcalign, function(x) x[,1]))
  matCorrected <- t(do.call(rbind, integrated.corrected.data[[2]]))
  rownames(matCorrected) <- integrated.corrected.data[[3]]
  mat2 <- matCorrected

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
  countsNormCorrected <- countsNormCorrected[!rm1]

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
  #Calculate gene counts signature (GCS) or the genes most correlated with gene counts
  ds2 <- sapply(1:nrow(mat2), function(x) ccaPP::corPearson(mat2[x,],countsNormCorrected))
  names(ds2) <- rownames(mat2)
  gcs2 <- apply(mat2[which(rownames(mat2) %in% names(rev(sort(ds2))[1:200])),],2,mean)

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

  gcs_regressed <- regressed(D, gcs2)
  gcs_diffused <- gcs_regressed
  cytotrace <- rank(gcs_diffused)

  #Getting plotting coordinates
  m <- do.call(rbind, integrated.corrected.data[[1]])
  m <- m[!rm1,]

  #filter
  filteredCells <- !(unlist(lapply(datasets, colnames)) %in% unlist(lapply(processedMats, rownames))[!rm1])

  return(list(data = mat2, CytoTRACE = cytotrace, GCS= gcs2, Counts = countsNormCorrected, coord = m, filteredCells = filter))
}
