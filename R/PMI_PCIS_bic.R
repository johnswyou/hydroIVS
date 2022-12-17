# library(KernSmooth)
# library(nnet)
#-------------------------------------------------------------------------------
standardize_dat <- function(Z) {
  Z.mean <- mean(Z)
  Z.sd <- stats::sd(Z)
  Z.stand <- (Z - Z.mean) / Z.sd
  return(Z.stand)
}
#-------------------------------------------------------------------------------
dot.prod <- function(s, d) {
  s1 <- s[1:d]
  s2 <- s[(d + 1):(2*d)]
  return(sum(s1 * s2))
}
#-------------------------------------------------------------------------------
density.estimation.d1 <- function(Z) {

  dens <- 1 / (1 + exp(-Z))
  return(dens)
}
#-------------------------------------------------------------------------------
# get.norm.sigma <- function(Z) {
#
#   Z <- as.matrix(Z)
#   d <- ncol(Z)
#   sd <- sqrt(apply(Z, 2, stats::var))
#   sigma <- sd * (4/(d + 2))^(1/(d + 4)) * N^(-1/(d + 4))
# }
#-------------------------------------------------------------------------------
kernel.est.uvn <- function(Z) {

  N <- length(Z)
  d <- 1
  # compute sigma & constant
  sigma <- stats::bw.nrd(Z)
  constant <- sqrt(2*pi) * sigma * N

  # Commence main loop
  dens <- vector()
  for(h in 1:N) {
    dis.Z <- (Z - Z[h])^2
    exp.dis <- exp(-dis.Z / (2*sigma^2))
    dens[h] <- sum(exp.dis) / constant
  }
  return(dens)
}
#-------------------------------------------------------------------------------
density.estimation.d2 <- function(Z.1, Z.2) {

  d.1 <- density.estimation.d1(Z.1)
  d.2 <- density.estimation.d1(Z.2) # used to be density.estimation.c1

  dens <- d.1 * d.2

  return(dens)
}
#-------------------------------------------------------------------------------
kernel.est.mvn <- function(Z) {

  # Compute covariance and determinant of cov
  N <- nrow(Z)
  d <- ncol(Z)
  Cov <- stats::cov(Z)
  det.Cov <- det(Cov)

  # Compute sigma & constant
  sigma <- (4/(d + 2))^(1/(d + 4)) * N^(-1/(d + 4))
  constant <- (sqrt(2*pi)*sigma)^d * sqrt(det.Cov) * N

  # Commence main loop
  dens <- vector()
  for(h in 1:N) {
    dist.val <- stats::mahalanobis(Z, center = Z[h,], cov = Cov)
    exp.dis <- exp(-dist.val / (2*sigma^2))
    dens[h] <- sum(exp.dis) / constant
  }

  return(dens)
}
#-------------------------------------------------------------------------------
pmi.calc <- function(Y, X, Y.bin = FALSE) {
  # Use Gaussian kernel to compute marginal densities F(Y) & F(X) and
  # joint density F(X,Y)
  N <- length(Y)
  pdf.X <- kernel.est.uvn(X)
  if(Y.bin) {
    pdf.Y <- density.estimation.d1(Y)
    pdf.XY <- density.estimation.d2(Y, X)
  } else {
    pdf.Y <- kernel.est.uvn(Y)
    pdf.XY <- kernel.est.mvn(cbind(Y, X))
  }
  calc <- log(pdf.XY / (pdf.Y * pdf.X))
  return(sum(calc) / N)
}
#-------------------------------------------------------------------------------
enp <- function(Z, linear = FALSE) {
  # Compute effective number of parameters

  Z <- as.matrix(Z)
  if(linear) {
    A.mat <- Z %*% t(Z)
    A.inv <- chol2inv(A.mat)
    H.mat <- t(Z) %*% A.inv %*% Z
    return(sum(diag(H.mat)))
  } else {
    N <- nrow(Z)
    d <- ncol(Z)
    sigma <- (4/(d + 2))^(1/(d + 4)) * N^(-1/(d + 4))
    Cov <- stats::cov(Z)
    det.Cov <- det(Cov)
    inv.Cov <- chol2inv(Cov)
    constant <- (sqrt(2*pi)*sigma)^d * sqrt(det.Cov)

    H.mat <- array(dim = c(N, N))
    for(h in 1:N) {
      dis.Z <- vector()
      for(j in 1:d) {
        dis.Z <- cbind(dis.Z, Z[h,j] - Z[,j])
      }
      if(d == 1) {
        exp.dis <- exp(-(dis.Z)^2 / (2*sigma^2)) / constant
      } else {
        T.dis.Z <- apply(dis.Z, 1, "%*%", inv.Cov)
        M.dis.Z <- apply(cbind(t(T.dis.Z), dis.Z), 1, dot.prod, d = d)
        exp.dis <- exp(-M.dis.Z / (2*sigma^2)) / constant
      }
      H.mat[h,] <- exp.dis / sum(exp.dis)
    }
  }
  return(sum(diag(H.mat)))
}
#-------------------------------------------------------------------------------
GRNN.mah <- function(Z.out, Z.in) {

  Z.in <- as.matrix(Z.in)
  N <- nrow(Z.in)
  d <- ncol(Z.in)
  sigma <- (4 / (d + 2))^(1 / (d + 4)) * N^(-1 / (d + 4))
  Cov <- stats::cov(Z.in)
  det.Cov <- det(Cov)
  constant <- (sqrt(2*pi)*sigma)^d * sqrt(det.Cov)

  # Commence main loop
  Z.hat <- vector()
  for(i in 1:N) {
    dis <- vector()
    for(j in 1:d) {
      dis <- cbind(dis, Z.in[-i,j] - Z.in[i,j])
    }
    if(d == 1) {
      exp.dis <- exp(-(dis)^2 / (2*sigma^2)) / constant
    } else {
      dist.val <- stats::mahalanobis(Z.in[-i,], center = Z.in[i,], cov = Cov)
      exp.dis <- exp(-dist.val / (2*sigma^2)) / constant
    }


    dot.A <- sum(Z.out[-i] * exp.dis)      # Sum numerator
    dot.B <- sum(exp.dis)              # Sum denominator
    dot.B <- max(dot.B, 1e-6)
    Z.hat[i] <- dot.A / dot.B
  }
  return(Z.hat)
}
#-------------------------------------------------------------------------------
GRNN.gauss <- function(Z.out, Z.in) {

  Z.in <- as.matrix(Z.in)
  N <- nrow(Z.in)
  d <- ncol(Z.in)
  sigma <- (4 / (d + 2))^(1 / (d + 4)) * N^(-1 / (d + 4))
  Cov <- stats::cov(Z.in)
  det.Cov <- det(Cov)
  constant <- (sqrt(2*pi)*sigma)^d * sqrt(det.Cov)

  # Commence main loop
  Z.hat <- vector()
  for(i in 1:N) {
    dist.val <- vector()
    for(j in 1:d) {
      dist.val <- cbind(dist.val, (Z.in[-i,j] - Z.in[i,j])^2)
    }
    if(d > 1) {
      dist.val <- rowSums(dist.val)   # Sum Euclidean distances
    }
    exp.dis <- exp(-dist.val / (2*sigma^2)) / constant


    dot.A <- sum(Z.out[-i] * exp.dis)      # Sum numerator
    dot.B <- sum(exp.dis)              # Sum denominator
    dot.B <- max(dot.B, 1e-6)
    Z.hat[i] <- dot.A / dot.B
  }
  return(Z.hat)
}
#-------------------------------------------------------------------------------
PNN <- function(Z.out, Z.in) {

  Z.in <- as.matrix(Z.in)
  N <- nrow(Z.in)
  d <- ncol(Z.in)
  sigma <- (4/(d + 2))^(1 /(d + 4)) * N^(-1/(d + 4))

  # Commence main loop
  Z.hat <- vector()
  class.vals <- unique(Z.out)
  class.A <- which(Z.out == class.vals[1])
  class.B <- which(Z.out == class.vals[2])
  for(i in 1:N) {
    if(d > 1) {
      # print(Z.in[class.A[class.A != i],])
      # print(matrix(Z.in[i,], nrow = 1))
      dist.A <- rowSums((Z.in[class.A[class.A != i],] - as.vector(matrix(Z.in[i,], nrow = 1)))^2) # added as.vector() by JY October 2022
      dist.B <- rowSums((Z.in[class.B[class.B != i],] - as.vector(matrix(Z.in[i,], nrow = 1)))^2) # added as.vector() by JY October 2022
    } else {
      dist.A <- (Z.in[i,] - Z.in[class.A,])^2
      dist.B <- (Z.in[i,] - Z.in[class.B,])^2
    }
    dist.A <- exp(-dist.A / (2*sigma^2))
    dist.B <- exp(-dist.B / (2*sigma^2))
    dot.A <- sum(dist.A)/length(class.A[class.A != i])
    dot.B <- sum(dist.B)/length(class.B[class.B != i])
    Z.hat[i] <- ifelse(dot.A > dot.B, class.vals[1], class.vals[2])
  }
  return(Z.hat)
}
#-------------------------------------------------------------------------------
#' @title Partial Correlation Input Selection (PCIS)
#' @description An implementation of the partial correlation input selection algorithm (May et al., 2008)
#' with the Bayesian Information Criterion (BIC) used for the stopping criterion. This function is essentially the same as
#' `PMI_grnn_framework_bic`, except that the Pearson correlation coefficient and multiple linear regression are used instead of
#' partial mutual information (PMI) and the generalied regresssion neural network (GRNN), respectively.
#'
#' This function was written by Greer B. Humphrey (School of Civil, Environmental, and Mining Engineering,
#' University of Adelaide, SA, 5005, Australia). The original code repository for this function
#' can be found at (https://github.com/gbhumphrey1/PMIS_PCIS).
#' @param x A matrix of input variables as columns and observations as rows \[N x D\]
#' @param y A vector containing target values
#' @param silent Verbosity, Default: `FALSE`
#' @return A data frame with selected inputs (integer indices) in the first column, and additional information in the other columns.
#' @references
#' Galelli, S., Humphrey, G. B., Maier, H. R., Castelletti, A., Dandy, G. C., &amp; Gibbs, M. S. (2014).
#' An evaluation framework for input variable selection algorithms for environmental data-driven models.
#' Environmental Modelling &amp; Software, 62, 33–51. https://doi.org/10.1016/j.envsoft.2014.08.015
#'
#' May, R. J., Maier, H. R., Dandy, G. C., &amp; Fernando, T. M. K. G. (2008). Non-linear variable selection for
#' artificial neural networks using partial mutual information. Environmental Modelling &amp; Software, 23(10-11),
#' 1312–1326. https://doi.org/10.1016/j.envsoft.2008.03.007
#' @examples
#' \dontrun{
#' if(interactive()){
#'  X <- matrix(rnorm(10*1000), ncol = 10)
#'  y <- rnorm(1000)
#'  res <- PCIS_framework_bic(X, y)
#'  res$Input # selected inputs
#'  }
#' }
#' @seealso
#'  \code{\link[stats]{lm.fit}}, \code{\link[stats]{cor}}
#' @rdname PCIS_framework_bic
#' @export
#' @importFrom stats lm.fit cor
PCIS_framework_bic <- function(x, y, silent = FALSE) {

  #** Stepwise PCIS**
  # Author : Greer Humphrey
  # School of Civil and Environmental Engineering
  # University of Adelaide, SA, 5005, Australia

  cpu.time <- proc.time()
  n.inputs <- ncol(x)
  n.data <- nrow(x)
  inp.names <- names(x)

  # record number of operations performed on data
  nfevals.1 <- 0 # number of linear model calls
  nfevals.2 <- 0 # number of correlation calls

  # Standardize the output data with zero mean and unit variance
  y.stand <- standardize_dat(y)
  y <- y.stand

  # Standardize the input data with zero mean and unit variance
  z.stand <- apply(x, 2, standardize_dat)
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
      # Reset the array with the output in it back to the original values
      y <- y.stand
      # Compute output with already selected inputs and output residuals
      #y.hat <- lm(y ~ z.in)$fitted.values
      y.hat <- stats::lm.fit(cbind(1, z.in), y)$fitted.values
      nfevals.1 <- nfevals.1 + 1

      # Calculate residuals u  Note: the y array is reset back to
      # the original values at start of next iteration (see above)
      u <- y - y.hat

      # Standardize output residuals
      u.stand <- standardize_dat(u)

      mean.y <- mean(y)
      AIC.k <- n.data * log(sum(u^2) / n.data) + 2*(n.selected + 1)
      BIC.score <- n.data * log(sum(u^2) / n.data) + log(n.data)*(n.selected + 1)
      RMSE <- sqrt(sum(u^2) / n.data)
      R2 <- 1 - sum(u^2)/sum((y - mean.y)^2)

      scores.0 <- data.frame(Input = input.tracker[tag],
                             PC = signif(PC[tag], 4),
                             AIC_k = signif(AIC.k, 4),
                             BIC = signif(BIC.score, 4),
                             RMSE = signif(RMSE, 4),
                             R2 = signif(R2, 4),
                             Hampel = signif(Z.score[tag], 4),
                             ModEvals = nfevals.1,
                             PCEvals = nfevals.2,
                             CPUtime = signif(as.numeric(proc.time()[2] - cpu.time[2]), 4),
                             ElapsedTime = signif(as.numeric(proc.time()[3] - cpu.time[3]), 4))

      scores <- rbind(scores, scores.0)
      # if(iter.1 == 2) {
      #   utils::write.table(scores.0, file = "", row.names = FALSE, col.names = TRUE, sep = ",")
      # } else {
      #   utils::write.table(scores.0, file = "", row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
      # }

      # Output results to screen and file
      if(!silent) print(scores, quote = FALSE)
      if(iter.1 == max.iter) break # break if reached maximum iterations
      if(iter.1 > 2){if (is.na(scores$BIC[iter.1-1] > scores$BIC[iter.1-2])
                         || scores$BIC[iter.1-1] > scores$BIC[iter.1-2]) break} # break if BIC has increased

      # Remove selected input from input data set and add to selected data set
      n.inputs <- n.inputs - 1
      input.tracker <- input.tracker[-tag]
      x <- as.matrix(x[,-tag])

      # Go through each input in turn to compute PMI
      # due to the already selected inputs
      z.res <- vector()
      Input <- vector()
      PC <- vector()
      for(current.inp in 1:n.inputs) {
        z <- x[,current.inp]

        # Now compute input residuals
        # z.hat <- lm(z ~ z.in)$fitted.values
        z.hat <- stats::lm.fit(cbind(1, z.in), z)$fitted.values
        nfevals.1 <- nfevals.1 + 1
        v <- z - z.hat

        # Standardize input residuals
        v.stand <- standardize_dat(v)

        # Store input residuals to use for bootstrapping
        z.res <- cbind(z.res, v.stand)
        # Compute correlation between input and output
        PC[current.inp] <- stats::cor(u.stand, v.stand)^2
        nfevals.2 <- nfevals.2 + 1
        Input[current.inp] <- input.tracker[current.inp]
      }
    } else {
      # Go through each input in turn to compute MI
      Input <- vector()
      PC <- vector()
      for(current.inp in 1:n.inputs) {
        z <- x[,current.inp]
        # Compute correlation between input and output
        PC[current.inp] <- stats::cor(y, z)^2
        nfevals.2 <- nfevals.2 + 1
        Input[current.inp] <- input.tracker[current.inp]
      }
    }

    # if(n.selected == 0) {
    #   utils::write.table(data.frame(Input = Input,
    #                          PC = format(PC, digits = 4)), file = "",
    #               row.names = FALSE, quote = FALSE)
    # }

    # Find input with highest PC score
    best.PC <- max(PC,na.rm=TRUE)
    tag <- which.max(PC)

    # Output the input that has been selected
    # if(!silent) cat("\nSelect input", inp.names[input.tracker[tag]], "....\n")

    # Calculate various scores
    Z.score <- MAD(PC)

    n.selected <- n.selected + 1
    z.in <- cbind(z.in, as.vector(x[,tag]))
  }

  # The program will only reach this point upon termination
  # cat("\nPCIS_BIC ROUTINE COMPLETED\n")

  # print results to console
  # print(scores, digits = 4,
  #       quote = FALSE, right = TRUE,
  #       row.names = TRUE)

  # optional write to file (comment/uncomment line below)
  # utils::write.table(scores,"PCIS_results.txt",sep="\t",row.names=FALSE)

  return(scores[1:(iter.1-2),])  # first column is

}
#-------------------------------------------------------------------------------
#' @title Partial Mutual Information Input Selection using the Generalized Regression Neural Network
#' @description From Galelli S., Humphrey G.B., Maier H.R., Castelletti A., Dandy G.C. and Gibbs M.S. (2014) <DOI: 10.1016/j.envsoft.2014.08.015>:
#'
#' The PMIS algorithm adopts a forward selection strategy (i.e. one variable is selected at each iteration) based on the
#' estimation of the PMI, which measures the partial dependence between each input variable and the output, conditional
#' on the inputs that have already been selected. To estimate the PMI, the algorithm adopts a kernel density estimator,
#' whose accuracy depends on the value of a smoothing parameter (or bandwidth) lambda. Similarly to Sharma (2000) and
#' Bowden et al. (2005a), the Gaussian reference bandwidth (Scott, 1992) is adopted, because of its simplicity and computational efficiency.
#' The calculation of PMI also requires estimation of residual information in the input variables once the effect of the already
#' selected inputs has been taken in consideration: this is done through the identification of a General Regression Neural Network (GRNN),
#' which is a nonlinear and nonparametric regression method (Li et al., 2014). In addition to lambda, the other parameter to be set is the
#' stopping criterion, which is based on the coefficient of determination R2: the PMIS algorithm is stopped when the selection of a further
#' input variable leads to a decrease of R2 in the underlying GRNN being identified.
#'
#' This function was written by Greer B. Humphrey (School of Civil, Environmental, and Mining Engineering,
#' University of Adelaide, SA, 5005, Australia). The original code repository for this function
#' can be found at (https://github.com/gbhumphrey1/PMIS_PCIS).
#' @param x A matrix of input variables as columns and observations as rows \[N x D\]
#' @param y A vector containing target values
#' @param ybin Is `y` binary?, Default: FALSE
#' @param silent Verbosity, Default: `FALSE`
#' @return A data frame with selected inputs (integer indices) in the first column, and additional information in the other columns.
#' @references
#' Galelli, S., Humphrey, G. B., Maier, H. R., Castelletti, A., Dandy, G. C., &amp; Gibbs, M. S. (2014).
#' An evaluation framework for input variable selection algorithms for environmental data-driven models.
#' Environmental Modelling &amp; Software, 62, 33–51. https://doi.org/10.1016/j.envsoft.2014.08.015
#'
#' May, R. J., Maier, H. R., Dandy, G. C., &amp; Fernando, T. M. K. G. (2008). Non-linear variable selection for
#' artificial neural networks using partial mutual information. Environmental Modelling &amp; Software, 23(10-11),
#' 1312–1326. https://doi.org/10.1016/j.envsoft.2008.03.007
#'
#' Sharma, A. (2000). Seasonal to interannual rainfall probabilistic forecasts for improved water supply management:
#' Part 1 — a strategy for system predictor identification. Journal of Hydrology, 239(1-4), 232–239.
#' https://doi.org/10.1016/s0022-1694(00)00346-2
#'
#' Bowden, G. J., Dandy, G. C., &amp; Maier, H. R. (2005). Input determination for neural network models in water resources applications.
#' part 1—background and methodology. Journal of Hydrology, 301(1-4), 75–92. https://doi.org/10.1016/j.jhydrol.2004.06.021
#'
#' Scott, D. W. (2015). Multivariate density estimation: Theory, practice, and visualization. John Wiley &amp; Sons.
#'
#' Li, X., Zecchin, A. C., &amp; Maier, H. R. (2014). Selection of smoothing parameter estimators for general regression neural
#' networks – applications to hydrological and water resources modelling. Environmental Modelling &amp; Software, 59, 162–186.
#' https://doi.org/10.1016/j.envsoft.2014.05.010
#' @examples
#' \dontrun{
#' if(interactive()){
#'  X <- matrix(rnorm(10*1000), ncol = 10)
#'  y <- rnorm(1000)
#'  res <- PMI_grnn_framework_bic(X, y)
#'  res$Input # selected inputs
#'  }
#' }
#' @seealso
#'  \code{\link[stats]{nls}}, \code{\link[stats]{predict}}
#'  \code{\link[utils]{write.table}}
#' @rdname PMI_grnn_framework_bic
#' @export
#' @importFrom stats nls predict
#' @importFrom utils write.table
PMI_grnn_framework_bic <- function(x, y, ybin = FALSE, silent = FALSE) {
  # Stepwise PMI
  # Use tabulated critical values to identify significant inputs

  # Used Gaussian kernel for both PDF estimation & GRNN predictions
  # Average MI/PMI computed without considering effect of negative values of MI
  cpu.time <- proc.time()
  n.inputs <- ncol(x)
  n.data <- nrow(x)
  inp.names <- names(x)

  I.crit <- c(50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 220, 240, 260,
              280, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000,
              0.2224, 0.2031, 0.1879, 0.1756, 0.1657, 0.1572, 0.1434, 0.1321,
              0.1237, 0.1166, 0.1103, 0.1055, 0.1005, 0.0965, 0.0928, 0.0896,
              0.0775, 0.0689, 0.0627, 0.0578, 0.0539, 0.0507, 0.0481, 0.0333,
              0.0268, 0.0230, 0.0204)
  dim(I.crit) <- c(27, 2)
  I.crit <- as.data.frame(I.crit)
  names(I.crit) <- c("n", "val")

  # Fit model to tabulated values
  I.crit.mod <- stats::nls(val ~ a*n^b, data = I.crit,
                    start = list(a = 1, b = 1))
  new.data <- data.frame(n = n.data, val = NA)
  I.crit.val <- stats::predict(I.crit.mod, newdata = new.data)

  # record number of operations performed on data
  nfevals.1 <- 0 # number of grnn calls
  nfevals.2 <- 0 # number of pmi calls

  # Standardize the output data with zero mean and unit variance
  if(!ybin) {
    y.stand <- standardize_dat(y)
    y <- y.stand
  } else {
    y.stand <- y
  }

  # Standardize the input data with zero mean and unit variance
  z.stand <- apply(x, 2, standardize_dat)
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

      # Reset the array with the output in it back to the original values
      y <- y.stand
      # Compute output with already selected inputs and output residuals
      if(!ybin) {
        y.hat <- GRNN.gauss(y, z.in)
        nfevals.1 <- nfevals.1 + 1
      } else {
        y.hat <- PNN(y, z.in)
        nfevals.1 <- nfevals.1 + 1
      }
      # Calculate residuals u  Note: the y array is reset back to
      # the original values at start of next iteration (see above)
      u <- y - y.hat

      # Standardize output residuals
      if(!ybin) {
        u.stand <- standardize_dat(u)
      } else {
        u.stand <- u
      }

      mean.y <- mean(y)
      p <- enp(z.in)
      AIC.score <- n.data * log(sum(u^2) / n.data) + 2*p
      AIC.k <- n.data * log(sum(u^2) / n.data) + 2*(n.selected + 1)
      BIC.score <- n.data * log(sum(u^2) / n.data) + log(n.data)*p
      RMSE <- sqrt(sum(u^2) / n.data)
      R2 <- 1 - sum(u^2)/sum((y - mean.y)^2)

      scores.0 <- data.frame(Input = input.tracker[tag],
                             PMI = signif(PMI[tag], 4),
                             MC_I_95 = signif(I.crit.val, 4),
                             AIC_p = signif(AIC.score, 4),
                             AIC_k = signif(AIC.k, 4),
                             BIC = signif(BIC.score, 4),
                             RMSE = signif(RMSE, 4),
                             R2 = signif(R2, 4),
                             Hampel = signif(Z.score[tag], 4),
                             ModEvals = nfevals.1,
                             PMIEvals = nfevals.2,
                             CPUtime = signif(as.numeric(proc.time()[2] - cpu.time[2]), 4),
                             ElapsedTime = signif(as.numeric(proc.time()[3] - cpu.time[3]), 4))

      scores <- rbind(scores, scores.0)
      # if(iter.1 == 2) {
      #   utils::write.table(scores.0, file = outfile, row.names = FALSE, col.names = TRUE, sep = ",")
      # } else {
      #   utils::write.table(scores.0, file = outfile, row.names = FALSE, col.names = FALSE, sep = ",", append = TRUE)
      # }
      #      write("", file = outfile, append = TRUE)

      # Output results to screen and file
      if(!silent) print(scores, quote = FALSE)
      if(iter.1 == max.iter) break # break if reached maximum iterations
      if(iter.1 > 2){if (is.na(scores$BIC[iter.1-1] > scores$BIC[iter.1-2])
                         || scores$BIC[iter.1-1] > scores$BIC[iter.1-2]) break} # break if BIC has increased

      # Remove selected input from input data set and add to selected data set
      n.inputs <- n.inputs - 1
      input.tracker <- input.tracker[-tag]
      x <- as.matrix(x[,-tag])

      # Go through each input in turn to compute PMI
      # due to the already selected inputs
      Input <- vector()
      PMI <- vector()
      for(current.inp in 1:n.inputs) {
        # Compute and standardise input residuals
        z <- x[,current.inp]
        z.hat <- GRNN.gauss(z, z.in)
        nfevals.1 <- nfevals.1 + 1
        v <- z - z.hat
        v.stand <- standardize_dat(v)
        # Compute PMI
        PMI[current.inp] <- pmi.calc(u.stand, v.stand, ybin)
        nfevals.2 <- nfevals.2 + 1
        Input[current.inp] <- input.tracker[current.inp]
      }
    } else {
      # Go through each input in turn to compute MI
      Input <- vector()
      PMI <- vector()
      for(current.inp in 1:n.inputs) {
        z <- x[,current.inp]
        PMI[current.inp] <- pmi.calc(y, z, ybin)
        nfevals.2 <- nfevals.2 + 1
        Input[current.inp] <- input.tracker[current.inp]
      }
    }

    # if(!silent) print(cbind(inp.names[Input], format(PMI, digits = 4)), quote = FALSE)
    # if(n.selected == 0) {
    #   utils::write.table(data.frame(Input = format(inp.names[Input]),
    #                          PMI = format(PMI, digits = 4)), file = mifile,
    #               row.names = FALSE, quote = FALSE)
    # }

    # Find input with highest PMI score
    best.PMI <- max(PMI,na.rm=TRUE)
    if (best.PMI <= 0){break}
    tag <- which.max(PMI)
    if(!silent) cat("\nSelect input", inp.names[input.tracker[tag]], "....\n")

    # Calculate various scores
    Z.score <- MAD(PMI)

    z.in <- cbind(z.in, as.vector(x[,tag]))
    n.selected <- n.selected + 1
  }

  # The program will only reach this point upon termination
  cat("\nPMIS_BIC ROUTINE COMPLETED\n")

  # print results to console
  print(scores, digits = 4,
        quote = FALSE, right = TRUE,
        row.names = TRUE)

  # optional write to file (comment/uncomment line below)
  # utils::write.table(scores,"PMIS_results.txt",sep="\t",row.names=FALSE)

  return(scores[1:(iter.1-2),])  # first column is the set of selected inputs
}
#-------------------------------------------------------------------------------
