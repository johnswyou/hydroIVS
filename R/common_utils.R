standardize_X <- function(X,...){

  X = as.matrix(X)
  Xmean = apply(X, MARGIN = 2,FUN = function(X) (mean(X,...))) # mean of X
  Xsd = apply(X, MARGIN = 2,FUN = function(X) (stats::sd(X,...))) # standard deviation of X

  Xs = t(t(sweep(X, 2, Xmean)) / Xsd) # standardized X

  # check for any columns that are completely NA and replace with Os
  test.nona = colSums(is.na(Xs)) == 0

  if (any(test.nona == FALSE)){
    Xs[,as.integer(which(test.nona == FALSE))] = 0
  }

  return(Xs)
}
# ------------------------------------------------------------------------------
MAD <- function(Z) {
  # To compute median absolute deviation
  #--------------------------------------

  # Compute median CMI/PMI values
  median.Z <- stats::median(Z)

  # Compute median absolute deviation of CMI/PMI values
  MAD.Z <- abs(Z - median.Z)

  # Compute S
  S <- 1.4826 * stats::median(MAD.Z)

  # Compute MAD Z.score
  Z.score <- MAD.Z / S

  return(Z.score)

}
