#' getBic
#'
#' @description Get the the BIC for a defined number of PCs
#'
#' @param eig Eigenvalues of the covariance matrix
#' @param N Number of observations
#' @param d Original number of dimensions
#' @param k Selected dimensions to calculate the BIC
#'
#' @return
#' @export
#'
#' @examples
getBic <- function(eig = NULL, N = NULL, d = NULL, k = NULL){

  # A
  A = -(N*0.5)*sum(log(eig[1:k]))

  # B
  v <-  sum(eig[(k+1):length(eig)])/(d-k)
  B <- -(N*(d-k)*0.5)*log(v)

  # C
  m = (d*k)-(k*(k+1)/2)
  C = -((m+k)*0.5)*log(N)

  return(A + B + C)
}

#' getBicPca
#'
#' @description Get the BIC for PCA
#'
#' @param pcaData PCA object
#' @param K Number of PCs to test
#'
#' @return
#' @export
#'
#' @examples
getBicPca <- function(pcaData = NULL, K = NULL){

  eig <- (pcaData$sdev)^2
  N = nrow(pcaData$x)
  d = length(eig)
  logL <- rep(NA, length(K))

  for(k in seq_len(K)){
    logL[k] <- getBic(eig = eig, N = N, d = d, k = k)
  }
  return(logL)
}
