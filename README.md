
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hydroIVS <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
<!-- badges: end -->

This package implements several input variable selection methods.
Currently implemented methods are:

- `none`: No input variable selection
- `boruta`: Input Variable Selection using the Boruta algorithm.
- `rrf`: Input Variable Selection (IVS) using feature importance scores
  from fitting a regularized random forrest (rrf) regression model.
- `ea_cmi_htc`: Edgeworth Approximation (EA) based Shannon Conditional
  Mutual Information (CMI) Input Variable Selection (IVS) using the
  Hampel Test Criterion (HTC) to identify significant inputs.
- `ea_cmi_tol`: Edgeworth Approximation (EA) based Shannon Conditional
  Mutual Information (CMI) Input Variable Selection (IVS) using ratio of
  CMI over Mutual Information (MI) to identify significant inputs.
- `knn_cmi_tol`: K nearest neighbour (KNN) based Shannon Conditional
  Mutual Information (CMI) Input Variable Selection (IVS) using ratio of
  CMI over Mutual Information (MI) to identify significant inputs.
- `knn_cmi_bi_tol`: Bias Improved (BI) K nearest neighbour (KNN) based
  Shannon Conditional Mutual Information (CMI) Input Variable Selection
  (IVS) using ratio of CMI over Mutual Information (MI) to identify
  significant inputs.
- `pmis_bic`: Partial Mutual Information Selection (PMIS) using the
  Bayesian Information Criterion (BIC) to identify significant inputs.
- `pcis_bic`: Partial Correlation Input Selection (PCIS) using the
  Bayesian Information Criterion (BIC) to identify significant inputs.

## Installation

You can install the development version of `hydroIVS` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("johnswyou/hydroIVS")
```

## Example

Here is a basic example:

``` r
library(hydroIVS)
# library(MASS)

set.seed(1648)

# Toy data set with 40 input features

# NOTE: I generated my covariance matrix below based on Daryl's answer in
# https://math.stackexchange.com/questions/357980/how-to-generate-random-symmetric-positive-definite-matrices-using-matlab
# vcovar <- matrix(rnorm(41^2), ncol=41)
# vcovar <- 0.5*(t(vcovar) + vcovar) + 41*diag(nrow=41)
# mu_vec <- 41*rnorm(41)
# X <- mvrnorm(n=1e3, mu=mu_vec, Sigma=vcovar)
# X <- X[, 2:41]
# y <- c(X[, 1])

# Give input features basic names
# colnames(X) <- paste0("X", 1:40)

library(wooldridge)

data("hprice3")

hprice2$lprice <- NULL
hprice2$lnox <- NULL
hprice2$lproptax <- NULL

head(hprice2)
#>   price crime  nox rooms dist radial proptax stratio lowstat
#> 1 24000 0.006 5.38  6.57 4.09      1    29.6    15.3    4.98
#> 2 21599 0.027 4.69  6.42 4.97      2    24.2    17.8    9.14
#> 3 34700 0.027 4.69  7.18 4.97      2    24.2    17.8    4.03
#> 4 33400 0.032 4.58  7.00 6.06      3    22.2    18.7    2.94
#> 5 36199 0.069 4.58  7.15 6.06      3    22.2    18.7    5.33
#> 6 28701 0.030 4.58  6.43 6.06      3    22.2    18.7    5.21

y <- hprice2$price
X <- hprice2[, 2:ncol(hprice2)]
X <- as.matrix(X)

# *******************************************
# Partial Correlation Input Selection (PCIS)
# *******************************************

# Bayesian Information Criterion (BIC) used to identify significant inputs.
ivsIOData(y, X, ivsm = "pcis_bic")
#> $sel_inputs
#> [1] 8 3 7 4 2 1 5
#> 
#> $names_sel_inputs
#> [1] "lowstat" "rooms"   "stratio" "dist"    "nox"     "crime"   "radial" 
#> 
#> $scores
#> [1] 0.52760 0.21560 0.11020 0.03114 0.06044 0.01381 0.01370

# ********************************************************
# Edgeworth Approximation (EA) based Shannon Conditional 
# Mutual Information (CMI) Input Variable Selection (IVS) 
# ********************************************************

# ivs_param indicates ratio of CMI over Mutual Information (MI) 
# used to identify significant inputs.
ivsIOData(y, X, ivsm = "ea_cmi_tol", ivs_param = 0.1)
#> 
#> EA_CMI_TOL ROUTINE COMPLETED
#>   Input    CMI    MI CMI.MI.ratio CMIevals CPUtime ElapsedTime
#> 1     3 0.6670 0.667       1.0000        8    0.00        0.03
#> 2     1 0.4288 1.096       0.3913       15    0.00        0.05
#> 3     8 0.1540 1.250       0.1233       21    0.02        0.08
#> 4     7 0.1275 1.377       0.0926       26    0.02        0.11
#> $sel_inputs
#> [1] 3 1 8
#> 
#> $names_sel_inputs
#> [1] "rooms"   "crime"   "lowstat"
#> 
#> $scores
#> [1] 0.6670 0.4288 0.1540

# ********************************************************
# K nearest neighbour (KNN) based Shannon Conditional 
# Mutual Information (CMI) Input Variable Selection (IVS) 
# ********************************************************

# ivs_param[1] indicates ratio of CMI over Mutual Information (MI) 
# used to identify significant inputs.

# ivs_param[2] indicates number of nearest neighbors
ivsIOData(y, X, ivsm = "knn_cmi_tol", ivs_param = c(0.1, 5))
#> 
#> KNN_CMI_TOL ROUTINE COMPLETED
#>   Input      CMI     MI CMI.MI.ratio CMIevals CPUtime ElapsedTime
#> 1     8 0.422900 0.4229      1.00000        8    0.11        1.44
#> 2     3 0.008007 0.4439      0.01804       15    0.25        4.06
#> $sel_inputs
#> [1] 8
#> 
#> $names_sel_inputs
#> [1] "lowstat"
#> 
#> $scores
#> [1] 0.4229
```

## References

Quilty, J., Adamowski, J., Khalil, B., & Rathinasamy, M. (2016).
Bootstrap rank‐ordered conditional mutual information (BROCMI): A
nonlinear input variable selection method for Water Resources Modeling.
Water Resources Research, 52(3), 2299–2326.
<https://doi.org/10.1002/2015wr016959>

Galelli S., Humphrey G.B., Maier H.R., Castelletti A., Dandy G.C. and
Gibbs M.S. (2014) An evaluation framework for input variable selection
algorithms for environmental data-driven models, Environmental Modelling
and Software, 62, 33-51, DOI: 10.1016/j.envsoft.2014.08.015.

Kursa, M. B., & Rudnicki, W. R. (2010). Feature Selection with the
Boruta Package. Journal of Statistical Software, 36(11), 1–13.
<https://doi.org/10.18637/jss.v036.i11>

H. Deng(2013). Guided Random Forest in the RRF Package. arXiv:1306.0237.

H. Deng and G. Runger (2012). Feature Selection via Regularized Trees.
The 2012 International Joint Conference on Neural Networks (IJCNN).

Kugiumtzis, D. (2013), Direct-coupling measure for nonuniform embedding,
Physical Review E, 87, 062918.

Tsimpiris, A., I. Vlachos, and D. Kugiumtzis (2012), Nearest neighbour
estimation of conditional mutual information in feature selection,
Expert Syst. Appl., 39, 697-708.

Marc Van Hulle. Edgeworth approximation of multivariate differential
entropy. Neural Computation, 17(9), 1903-1910, 2005.

I. Vlachos, D. Kugiumtzis, “Non-uniform state space reconstruction and
coupling detection”, Physical Review E, Vol 82, 016207, 2010

May, R. J., H. R. Maier, G. C. Dandy, and T. Fernando (2008a),
Non-linear variable selection for artificial neural networks using
partial mutual information, Environ. Modell. Software, 23(10-11),
1312-1326.
