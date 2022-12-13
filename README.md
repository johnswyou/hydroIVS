
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
library(MASS)

# Toy data set with 40 input features
vcovar <- matrix(rnorm(40^2), ncol=40)
vcovar <- t(vcovar) %*% vcovar
mu_vec <- rnorm(40)
X <- mvrnorm(n=1e3, mu=mu_vec, Sigma=vcovar)
y <- rnorm(1e3)

# Give input features basic names
colnames(X) <- paste0("X", 1:40)

# *******************************************
# Partial Correlation Input Selection (PCIS)
# *******************************************

# Bayesian Information Criterion (BIC) used to identify significant inputs.
ivsIOData(y, X, ivsm = "pcis_bic")
#> $sel_inputs
#> [1] 5
#> 
#> $names_sel_inputs
#> [1] "X5"
#> 
#> $scores
#> [1] 0.005146

# ********************************************************
# Edgeworth Approximation (EA) based Shannon Conditional 
# Mutual Information (CMI) Input Variable Selection (IVS) 
# ********************************************************

# ivs_param indicates ratio of CMI over Mutual Information (MI) 
# used to identify significant inputs.
ivsIOData(y, X, ivsm = "ea_cmi_tol", ivs_param = 0.1)
#> 
#> EA_CMI_TOL ROUTINE COMPLETED
#>   Input      CMI       MI CMI.MI.ratio CMIevals CPUtime ElapsedTime
#> 1     8 0.004716 0.004716      1.00000       40    0.00        0.10
#> 2     5 0.004476 0.009191      0.48690       79    0.00        0.21
#> 3    23 0.004453 0.013640      0.32640      117    0.00        0.38
#> 4    37 0.003249 0.016890      0.19230      154    0.08        0.58
#> 5    30 0.002877 0.019770      0.14550      190    0.14        0.85
#> 6    10 0.002701 0.022470      0.12020      225    0.23        1.19
#> 7    36 0.002456 0.024930      0.09854      259    0.31        1.60
#> $sel_inputs
#> [1]  8  5 23 37 30 10
#> 
#> $names_sel_inputs
#> [1] "X8"  "X5"  "X23" "X37" "X30" "X10"
#> 
#> $scores
#> [1] 0.004716 0.004476 0.004453 0.003249 0.002877 0.002701

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
#> NULL
#> $sel_inputs
#> NULL
#> 
#> $names_sel_inputs
#> character(0)
#> 
#> $scores
#> NULL
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
