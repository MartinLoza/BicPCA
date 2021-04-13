
#' FindNumberPC
#'
#' @param object object for which to find the number of PCs. (options: Seurat, pca)
#' @param dim2Test Number of dimensions to test
#' @param returnBIC Return the BIC values used to select the number of PCs.
#'
#' @return The selected number of PCs.
#' @export
#'
#' @rdname FindNumberPC
FindNumberPC <- function(object = NULL, dim2Test = NULL, returnBIC = FALSE, ...){
  UseMethod("FindNumberPC")
}

#' @rdname FindNumberPC
#' @export
FindNumberPC.Seurat <- function(object = NULL, dim2Test = NULL, returnBIC = FALSE, ...){

  features <- VariableFeatures(object)

  # We use the Seurat scaling to get comparable results. If we manually scale the data we could get similar but different results.
  scaData <- t(object[["RNA"]]@scale.data)

  if(length(scaData) == 0)
    stop("The Seurat object must be scaled.")

  if(!exists("eigVal"))
    eigVal <- eigen(cov(scaData))$values

  if(is.null(dim2Test))
    dim2Test <- min(dim(scaData))

  bicPCA <- -getBicPca(eig = eigVal, N = ncol(scaData), K = dim2Test)
  selPCA <- which.min(bicPCA)

  if(returnBIC)
    output <- list(BIC = bicPCA, selPCA = selPCA)
  else
    output <- selPCA

  return(output)
}
