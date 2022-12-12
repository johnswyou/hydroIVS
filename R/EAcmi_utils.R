#-------------------------------------------------------------------------------
whiten_E0 <- function(X,...){

  X = as.matrix(X)
  Xmean = apply(X, MARGIN = 2,FUN = function(X) (mean(X,...)))
  Xw = t(t(sweep(X, 2, Xmean))) # transfom X to zero mean

  # check for any columns that are completely NA and replace with Os
  test.nona = colSums(is.na(Xw)) == 0

  if (any(test.nona == FALSE)){
    Xw[,as.integer(which(test.nona == FALSE))] = 0
  }

  return(Xw)
}
#-------------------------------------------------------------------------------
std_X <- function(X,...){

  X = as.matrix(X)
  Xsd = apply(X, MARGIN = 2,FUN = function(X) (stats::sd(X,...))) # standard deviation of X

  return(Xsd)
}
#-------------------------------------------------------------------------------
H_whiten <- function(s) {

  # take  scaling into account via the entropy transformation rule [ H(Wz) = H(z)+log(|det(W)|) ]
  H_w = log(prod(s))

  return(H_w)
}
#-------------------------------------------------------------------------------
H_normal <- function(X) {

  d = ncol(X) # number of variables

  # Shannon entropy of a normal variable with stats::cov(X) covariance
  H_w = log(det(stats::cov(X)))/2 + (d/2 * log(2*pi)) + d/2

  return(H_w)
}
#-------------------------------------------------------------------------------
Edgeworth_t1_t2_t3 <- function(X) {

  # Computes the three kappa_ijk := E[x_i x_j x_k] based terms (t1,t2,t3) in the
  # Edgeworth expansion based entropy estimator

  X = as.matrix(X)
  d = ncol(X) # number of variables

  # t1
  t1 = 0
  for(i in 1:d){ # d terms
    kappa_iii = mean(X[,i]^3)
    t1 = t1 + kappa_iii^2
  }

  # t2
  t2 = 0

  if (d > 1){ # 2*nchoosek(d,2) terms
    for(i in 1:(d-1)){ # for case where i<j
      for(j in (i+1):d){
        kappa_iij = mean((X[,i]^2) * X[,j])
        t2 = t2 + kappa_iij^2
      }
    }

    for(j in 1:(d-1)){ # for case where j<i
      for(i in (j+1):d){
        kappa_iij = mean((X[,i]^2) * X[,j])
        t2 = t2 + kappa_iij^2
      }
    }
  }
  t2 = 3 * t2

  # t3
  t3 = 0
  if(d > 2){ # nchoosek(d,3) terms
    for(i in 1:(d-2)){ # i<j<k
      for(j in (i+1):(d-1)){
        for(k in (j+1):d){
          kappa_ijk = mean(X[,i] * X[,j] * X[,k])
          t3 = t3 + kappa_ijk^2
        }
      }
    }
  }
  t3 = t3 / 6

  return(list(
    t1 = t1,
    t2 = t2,
    t3 = t3
  ))

}
#-------------------------------------------------------------------------------
H_EdgeworthApprox <- function(X) {

  # Computes the Edgeworth Approximation (EA)-based differential Shannon entropy

  Xw = whiten_E0(X) # whiten dataset (i.e., make each column mean of 0), this step does not change the Shannon entropy of the variable
  s = std_X(Xw) # get standard deviation of whitened dataset
  Y = standardize_X(X) # standardize dataset (i.e., make each column of unit standard deviation)
  H_w = H_whiten(s) # get log-product of dataset standard deviation (see H_whiten)
  H_n = H_normal(Y) # get Shannon entropy of a normal variable with stats::cov(Y) covariance
  ea = Edgeworth_t1_t2_t3(Y) # compute the three kappa_ijk := E[x_i x_j x_k] based terms (t1,t2,t3)
  H_ea = (H_n - (ea[["t1"]] + ea[["t2"]] + ea[["t3"]])/12) + H_w # EA differential Shannon entropy

  # version 2: using det(A) = prod(eigenvalues(A)) $ seems less accurate then method above (tested
  # on a 5-, 10-, and 20-dimensional normal distribution, version 2 is always results in lower
  # accuracy against the theoretical entropy)

  # H_ea = (H_n - (ea[["t1"]] + ea[["t2"]] + ea[["t3"]])/12) + log(prod(eigen(stats::cov(X))$values))

  return(H_ea)
}
#-------------------------------------------------------------------------------
MI_EdgeworthApprox <- function(Y,X) {

  # Computes the Edgeworth Approximation (EA)-based Shannon Mutual Information
  MI_ea = H_EdgeworthApprox(Y) + H_EdgeworthApprox(X) - H_EdgeworthApprox(cbind(Y,X)) # EA Shannon mutual information

  MI_ea = max(0,Re(MI_ea))
  return(MI_ea)
}
#-------------------------------------------------------------------------------
CMI_EdgeworthApprox <- function(Y,X,Z) {

  # I(Y;X|Z) = I(Y;(X,Z)) - I(Y;Z) = H(Y,Z) + H(X,Z) - H(Z) - H(Y,X,Z)

  # Computes the Edgeworth Approximation (EA)-based Shannon Conditional Mutual Information
  CMI_ea = (H_EdgeworthApprox(cbind(Y,Z)) + H_EdgeworthApprox(cbind(X,Z))
            - H_EdgeworthApprox(Z)- H_EdgeworthApprox(cbind(Y,X,Z))) # EA Shannon conditional mutual information

  CMI_ea = max(0,Re(CMI_ea))
  return(CMI_ea)
}
