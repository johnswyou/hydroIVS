# library(FNN)
# library(RANN)

# Reference(s):
#
# Quilty, J., J. Adamowski, B. Khalil, and M. Rathinasamy (2016), Bootstrap rank-
# ordered conditional mutual information (broCMI): A nonlinear input variable
# selection method for water resources modeling, Water Resour. Res., 52,
# doi:10.1002/2015WR016959.
#
# # Gao, W., O. Sewoong, and P. Viswanath (2018), Demystifying Fixed k-Nearest
# Neighbor Information Estimators, IEEE Trans. Inform. Theory, doi:
# 10.1109/TIT.2018.2807481.
#
# Kugiumtzis, D. (2013), Direct-coupling measure for nonuniform embedding,
# Physical Review E, 87, 062918.
#
# Tsimpiris, A., I. Vlachos, and D. Kugiumtzis (2012), Nearest neighbour
# estimation of conditional mutual information in feature selection, Expert
# Syst. Appl., 39, 697-708.

#-------------------------------------------------------------------------------
# %function [RM,ecC] = PMIME(allM,Lmax,T,nnei,A,showtxt)
# % function [RM,ecC] = PMIME(allM,Lmax,T,nnei,A,showtxt)
# % PMIME (Partial Mutual Information on Mixed Embedding)
# % computes the measure R_{X->Y|Z} for all combinations of X and Y time
# % series from the multivariate time series given in matrix 'allM', of size
# % N x K, where Z contains the rest K-2 time series.
# % The components of X,Y, and Z, are found from a mixed embedding aiming at
# % explaining Y. The mixed embedding is formed by using the progressive
# % embedding algorithm based on conditional mutual information (CMI).
# % CMI is estimated by the method of nearest neighbors (Kraskov's method).
# % The function is the same as PMIMEsig.m but defines the stopping criterion
# % differently, using a fixed rather than adjusted threshold. Specifically,
# % the algorithm terminates if the contribution of the selected lagged
# % variable in explaining the future response state is small enough, as
# % compared to a threshold 'A'. Concretely, the algorithm terminates if
# %        I(x^F; w| wemb) / I(x^F; w,wemb) <= A
# % where I(x^F; w| wemb) is the CMI of the selected lagged variable w and
# % the future response state x^F given the current mixed embedding vector,
# % and I(x^F; w,wemb) is the MI between x^F and the augmented mixed
# % embedding vector [wemb w].
# % We experienced that in rare cases the termination condition is not
# % satisfied and the algorithm does not terminate. Therefore we included a
# % second condition for termination of the algorithm when the ratio
# % I(x^F; w| wemb) / I(x^F; w,wemb) increases in the last two embedding
# % cycles.
# % The derived R measure indicates the information flow of time series X to
# % time series Y conditioned on the rest time series in Z. The measure
# % values are stored in a K x K matrix 'RM' and given to the output, where
# % the value at position (i,j) indicates the effect from i to j (row to
# % col), and the (i,i) components are zero.
# % INPUTS
# % - allM : the N x K matrix of the K time series of length N.
# % - Lmax : the maximum delay to search for X and Y components for the mixed
# %          embedding vector [default is 5].
# % - T    : T steps ahead that the mixed embedding vector has to explain.
# %          Note that if T>1 the future vector is of length T and contains
# %          the samples at times t+1,..,t+T [dafault is 1].
# % - nnei : number of nearest neighbors for density estimation [default is 5]
# % - A    : the threshold for the ratio of CMI over MI of the lagged variables
# %          for the termination criterion.
# % - showtxt : if 0 or negative do not print out anything,
# %             if 1 print out the response variable index at each run,
# %             if 2 or larger print also info for each embedding cycle [default is 1].
# % OUTPUTS
# % - RM   : A K x K matrix containing the R values computed by PMIME using
# %          the surrogates for setting the stopping criterion.
# % - ecC  : cell array of K components, where each component is a matrix of
# %          size E x 5, and E is the number of embedding cycles. For each
# %          embedding cycle the following 5 results are stored:
# %          1. variable index, 2. lag index, 3. CMI of the selected lagged
# %          variable w and the future response state x^F given the current
# %          mixed embedding vector, I(x^F; w| wemb). 4. MI between x^F and
# %          the augmented mixed embedding vector [wemb w], I(x^F; w,wemb).
# %          5. The ration of 3. and 4.: I(x^F; w| wemb)/I(x^F; w,wemb)
# %
# %     Copyright (C) 2015 Dimitris Kugiumtzis
# %
# %     This program is free software: you can redistribute it and/or modify
# %     it under the terms of the GNU General Public License as published by
# %     the Free Software Foundation, either version 3 of the License, or
# %     (at your option) any later version.
# %
# %     This program is distributed in the hope that it will be useful,
# %     but WITHOUT ANY WARRANTY; without even the implied warranty of
# %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# %     GNU General Public License for more details.
# %
# %     You should have received a copy of the GNU General Public License
# %     along with this program.  If not, see <http://www.gnu.org/licenses/>.
# %
# %=========================================================================
# % Reference : D. Kugiumtzis, "Direct coupling information measure from
# %             non-uniform embedding", Physical Review E, Vol 87, 062918,
# %             2013
# %             I. Vlachos, D. Kugiumtzis, "Non-uniform state space
# %             reconstruction and coupling detection", Physical Review E,
# %             Vol 82, 016207, 2010
# % Link      : http://users.auth.gr/dkugiu/
# %=========================================================================
#-------------------------------------------------------------------------------
vd <- function(d) {

  # Compute the volume of unit l_2 ball in d dimensional space
  vol2d = d*log(2*gamma(1 + 0.5)) - log(gamma(1 + 0.5*d))

  return(vol2d)
}
#-------------------------------------------------------------------------------
MI_BI_KNN <- function(Y,X,k,psi.nn,log.nsamp) {

  # Computes the K nearest neighbour (KNN)-based Shannon Mutual Information

  # get max distances
  YX = cbind(Y,X)
  d_YX = ncol(YX)
  maxdists = knn_max_dist(YX,k+1)

  # find points within max distance in projected space (add small noise
  # in case some discrete inputs)
  smallnoise = 10^-15
  npY = knn_radius(Y,maxdists + smallnoise)-1
  npX = knn_radius(X,maxdists + smallnoise)-1
  log.both = log(cbind(npY,npX))
  vd_xy = vd(d_YX)
  vd_x = vd( ncol(as.matrix(X)))
  vd_y = vd( ncol(as.matrix(Y)))
  vd_correction = vd_x + vd_y - vd_xy

  MI_nn = vd_correction + psi.nn + log.nsamp - mean(rowSums(log.both)) # MI (Eq. 40 in Gao et al. [2018])

  MI_nn = max(0,Re(MI_nn))
  return(MI_nn)
}
#-------------------------------------------------------------------------------
CMI_BI_KNN <- function(Y,X,Z,k,psi.nn,log.nsamp) {

  # I(Y;X|Z) = I(Y;(X,Z)) - I(Y;Z) = H(Y,Z) + H(X,Z) - H(Z) - H(Y,X,Z)

  # Computes the K nearest neighbour (KNN)-based Shannon Conditional Mutual Information

  # get max distances
  YXZ = cbind(Y,X,Z)
  d_YXZ = ncol(YXZ)
  maxdists = knn_max_dist(YXZ,k+1)

  # find points within max distance in projected space (add small noise
  # in case some discrete inputs)
  smallnoise = 10^-15
  XZ = cbind(X,Z)
  YZ = cbind(Y,Z)
  d_XZ = ncol(XZ)
  d_YZ = ncol(YZ)
  npZ = knn_radius(Z,maxdists + smallnoise)-1
  npXZ = knn_radius(XZ,maxdists + smallnoise)-1
  npYZ = knn_radius(YZ,maxdists + smallnoise)-1
  npY = knn_radius(Y,maxdists + smallnoise)-1
  log.npZ = log(npZ)
  log.npXZ = log(npXZ)
  log.npYZ = log(npYZ)
  log.npY = log(npY)
  log.three = cbind(log.npYZ,log.npXZ,-log.npZ)
  log.two = cbind(log.npY,log.npXZ)
  vd_xyz = vd(d_YXZ)
  vd_xz = vd(d_XZ)
  vd_yz = vd(d_YZ)
  vd_z = vd( ncol(as.matrix(Z)))
  vd_y = vd( ncol(as.matrix(Y)))
  vd_correction1 = vd_xz + vd_yz - vd_z - vd_xyz
  vd_correction2 = vd_y + vd_xz - vd_xyz

  CMI_nn = vd_correction1 + psi.nn - mean(rowSums(log.three)) # CMI (see Eq. 42 in Gao et al. [2018] and Eq. 10 in Vlachos and Kugiumtzis [2010])
  MIp_nn = vd_correction2 + psi.nn + log.nsamp - mean(rowSums(log.two)) # projected MI (see Eq. 10-12 in Vlachos and Kugiumtzis [2010])

  CMI_nn = max(0,Re(CMI_nn))
  MIp_nn = max(0,Re(MIp_nn))

  return(list(
    CMI = CMI_nn,
    MI = MIp_nn
  ))
}
#-------------------------------------------------------------------------------
KNNcmi_bi_framework_tol <- function(x, y, thresh,k=5, silent = FALSE) {
  # Stepwise KNNcmi_bi_tol
  # Use ratio of CMI over MI to identify significant inputs

  # Used Edgeworth Approximation-based conditional mutual information;
  # MI/CMI given value of 0 in case of negative values of MI/CMI
  # (which is due to bias)
  cpu.time <- proc.time()
  n.inputs <- ncol(x)
  n.data <- nrow(x)
  inp.names <- names(x)

  # record number of operations performed on data
  nfevals.2 <- 0 # number of cmi calls

  # Standardize the output data with zero mean and unit variance
  y.stand <- standardize_X(y)
  y <- y.stand

  # Standardize the input data with zero mean and unit variance
  z.stand <- standardize_X(x)
  x <- z.stand

  # Mock up an array to keep track of unselected inputs
  input.tracker <- 1:n.inputs
  counter <- n.inputs
  n.selected <- 0
  z.in <- vector()
  scores <- NULL
  max.iter <- ncol(x)+1

  # precompute values beforehand...
  psik = psigamma(k)
  logN = log(n.data)

  for(iter.1 in 1:max.iter) {
    if(n.selected > 0) {

      scores.0 <- data.frame(Input = input.tracker[tag],
                             CMI = signif(CMI[tag], 4),
                             MI = signif(MI.aug, 4),
                             CMI.MI.ratio = signif(CMI.MI.ratio, 4),
                             CMIevals = nfevals.2,
                             CPUtime = signif(as.numeric(proc.time()[2] - cpu.time[2]), 4),
                             ElapsedTime = signif(as.numeric(proc.time()[3] - cpu.time[3]), 4))

      scores <- rbind(scores, scores.0)

      # Output results to screen and file
      if(!silent) print(scores, quote = FALSE)

      # stopping criteria

      # max iterations
      if(iter.1 == max.iter) break # break if reached maximum iterations

      if(iter.1 > 2){

        if( is.na(CMI.MI.ratio <= thresh)
            || CMI.MI.ratio <= thresh) break} # if ratio less than threshold

      if(iter.1 > 3){
        if(
          (scores$CMI.MI.ratio[iter.1-1] > scores$CMI.MI.ratio[iter.1-2]) &
          (scores$CMI.MI.ratio[iter.1-2] > scores$CMI.MI.ratio[iter.1-3])
        ) break} # if ratio keeps increasing

      # Remove selected input from input data set and add to selected data set
      n.inputs <- n.inputs - 1
      input.tracker <- input.tracker[-tag]
      x <- as.matrix(x[,-tag])

      # Go through each input in turn to compute CMI
      # due to the already selected inputs
      Input <- vector()
      CMI <- vector()
      MIp <- vector()
      for(current.inp in 1:n.inputs) {
        z <- x[,current.inp]

        # Compute CMI
        CMI_MIp_nn = CMI_BI_KNN(y, z, z.in, k, psik, logN)
        CMI[current.inp] <- CMI_MIp_nn$CMI
        MIp[current.inp] <- CMI_MIp_nn$MI
        nfevals.2 <- nfevals.2 + 1
        Input[current.inp] <- input.tracker[current.inp]
      }
    } else {
      # Go through each input in turn to compute MI
      Input <- vector()
      CMI <- vector()
      for(current.inp in 1:n.inputs) {
        z <- x[,current.inp]
        CMI[current.inp] <- MI_BI_KNN(y, z, k, psik, logN)
        nfevals.2 <- nfevals.2 + 1
        Input[current.inp] <- input.tracker[current.inp]
      }
    }

    # Find input with highest CMI score
    best.CMI <- max(CMI,na.rm=TRUE)
    if (best.CMI <= 0){break}
    tag <- which.max(CMI)

    if(iter.1==1){best.MIp = best.CMI} else{best.MIp = MIp[tag]}


    if(!silent) cat("\nSelect input", inp.names[input.tracker[tag]], "....\n")

    z.in <- cbind(z.in, as.vector(x[,tag]))
    n.selected <- n.selected + 1

    # calculate ratio between CMI between target and current selected input
    # conditioned on previously selected inputs and MI of target using
    # agumented input set

    MI.aug = best.MIp
    CMI.MI.ratio = best.CMI / MI.aug

  }

  # The program will only reach this point upon termination
  cat("\nKNN_CMI_TOL ROUTINE COMPLETED\n")

  # print results to console
  print(scores, digits = 4,
        quote = FALSE, right = TRUE,
        row.names = TRUE)

  # optional write to file (comment/uncomment line below)
  utils::write.table(scores,"EA_CMI_ivs_results.txt",sep="\t",row.names=FALSE)


  return(scores[1:(iter.1-2),])  # first column is the set of selected inputs
}
#-------------------------------------------------------------------------------
