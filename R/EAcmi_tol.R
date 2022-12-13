# library(KernSmooth)
# library(nnet)

# Reference(s):
#
# Quilty, J., J. Adamowski, B. Khalil, and M. Rathinasamy (2016), Bootstrap rank-
# ordered conditional mutual information (broCMI): A nonlinear input variable
# selection method for water resources modeling, Water Resour. Res., 52,
# doi:10.1002/2015WR016959.
#
# Kugiumtzis, D. (2013), Direct-couplig measure for nonuniform embedding,
# Physical Review E, 87, 062918.
#
# Tsimpiris, A., I. Vlachos, and D. Kugiumtzis (2012), Nearest neighbour
# estimation of conditional mutual information in feature selection, Expert
# Syst. Appl., 39, 697-708.

#-------------------------------------------------------------------------------
# %function [H] = HShannon_Edgeworth_estimation(Y,co)
# %Estimates the Shannon entropy (H) of Y using Edgeworth expansion.
# %
# %We use the naming convention 'H<name>_estimation' to ease embedding new entropy estimation methods.
# %
# %INPUT:
#   %   Y: Y(:,t) is the t^th sample.
# %  co: entropy estimator object.
# %
# %REFERENCE:
#   %   Marc Van Hulle. Edgeworth approximation of multivariate differential entropy. Neural Computation, 17(9), 1903-1910, 2005.
#
# %Copyright (C) 2012- Zoltan Szabo ("http://www.gatsby.ucl.ac.uk/~szabo/", "zoltan (dot) szabo (at) gatsby (dot) ucl (dot) ac (dot) uk")
# %
# %This file is part of the ITE (Information Theoretical Estimators) toolbox.
# %
# %ITE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
# %the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# %
# %This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# %MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# %
# %You should have received a copy of the GNU General Public License along with ITE. If not, see <http://www.gnu.org/licenses/>.
#
# %co.mult:OK. The information theoretical quantity of interest can be (and is!) estimated exactly [co.mult=1]; the computational complexity of the estimation is essentially the same as that of the 'up to multiplicative constant' case [co.mult=0]. In other words, the estimation is carried out 'exactly' (instead of up to 'proportionality').

# %function [t1,t2,t3] = Edgeworth_t1_t2_t3(Y)
# %Computes the three kappa_ijk := E[x_i x_j x_k] based terms (t1,t2,t3) in the Edgeworth expansion based entropy estimator, see 'HShannon_Edgeworth_estimation.m'.
# %
# %INPUT:
#   %  Y: Y(:,t) is the t^th sample
#
# %Copyright (C) 2012- Zoltan Szabo ("http://www.gatsby.ucl.ac.uk/~szabo/", "zoltan (dot) szabo (at) gatsby (dot) ucl (dot) ac (dot) uk")
# %
# %This file is part of the ITE (Information Theoretical Estimators) toolbox.
# %
# %ITE is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
# %the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# %
# %This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# %MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# %
# %You should have received a copy of the GNU General Public License along with ITE. If not, see <http://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------------
EAcmi_framework_tol <- function(x, y,thresh, silent = FALSE) {
  # Stepwise EAcmi_tol
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

      # check stopping conditions

      # max iterations
      if(iter.1 == max.iter) break # break if reached maximum iterations

      # stopping criteria
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
      for(current.inp in 1:n.inputs) {
        z <- x[,current.inp]

        # Compute CMI
        CMI[current.inp] <- CMI_EdgeworthApprox(y, z, z.in)
        nfevals.2 <- nfevals.2 + 1
        Input[current.inp] <- input.tracker[current.inp]
      }
    } else {
      # Go through each input in turn to compute MI
      Input <- vector()
      CMI <- vector()
      for(current.inp in 1:n.inputs) {
        z <- x[,current.inp]
        CMI[current.inp] <- MI_EdgeworthApprox(y, z)
        nfevals.2 <- nfevals.2 + 1
        Input[current.inp] <- input.tracker[current.inp]
      }
    }

    # Find input with highest CMI score
    best.CMI <- max(CMI,na.rm=TRUE)
    if (best.CMI <= 0){break}
    tag <- which.max(CMI)
    if(!silent) cat("\nSelect input", inp.names[input.tracker[tag]], "....\n")

    z.in <- cbind(z.in, as.vector(x[,tag]))
    n.selected <- n.selected + 1

    # calculate ratio between CMI between target and current selected input
    # conditioned on previously selected inputs and MI of target using
    # agumented input set

    MI.aug = MI_EdgeworthApprox(y, z.in)
    CMI.MI.ratio = best.CMI / MI.aug

  }

  # The program will only reach this point upon termination
  cat("\nEA_CMI_TOL ROUTINE COMPLETED\n")

  # print results to console
  print(scores, digits = 4,
        quote = FALSE, right = TRUE,
        row.names = TRUE)

  # optional write to file (comment/uncomment line below)
  # utils::write.table(scores,"EA_CMI_ivs_results.txt",sep="\t",row.names=FALSE)


  return(scores[1:(iter.1-2),])  # first column is the set of selected inputs
}
#-------------------------------------------------------------------------------
