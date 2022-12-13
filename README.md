
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

# Toy data set with 100 input features
X <- matrix(rnorm(1e4*100), ncol=100)
y <- rnorm(1e4)

# Give input features basic names
colnames(X) <- paste0("X", 1:100)

# *******************************************
# Partial Correlation Input Selection (PCIS)
# *******************************************

# Bayesian Information Criterion (BIC) used to identify significant inputs.
ivsIOData(y, X, ivsm = "pcis_bic")
#> $sel_inputs
#> [1] 17
#> 
#> $names_sel_inputs
#> [1] "X17"
#> 
#> $scores
#> [1] 0.0008484

# ********************************************************
# Edgeworth Approximation (EA) based Shannon Conditional 
# Mutual Information (CMI) Input Variable Selection (IVS) 
# ********************************************************

# ivs_param indicates ratio of CMI over Mutual Information (MI) 
# used to identify significant inputs.
ivsIOData(y, X, ivsm = "ea_cmi_tol", ivs_param = 0.1)
#> 
#> EA_CMI_TOL ROUTINE COMPLETED
#>   Input       CMI        MI CMI.MI.ratio CMIevals CPUtime ElapsedTime
#> 1    48 0.0008005 0.0008005       1.0000      100    0.36        0.83
#> 2    17 0.0005960 0.0013960       0.4268      199    0.77        2.44
#> 3    79 0.0005430 0.0019390       0.2800      297    1.13        4.49
#> 4    43 0.0005314 0.0024710       0.2151      394    1.55        7.30
#> 5    24 0.0005176 0.0029880       0.1732      490    2.19       11.30
#> 6    70 0.0005028 0.0034910       0.1440      585    3.06       16.60
#> 7    58 0.0005064 0.0039980       0.1267      679    3.72       23.22
#> 8     6 0.0004638 0.0044610       0.1040      772    4.78       31.52
#> 9    37 0.0004305 0.0048920       0.0880      864    5.56       41.41
#> $sel_inputs
#> [1] 48 17 79 43 24 70 58  6
#> 
#> $names_sel_inputs
#> [1] "X48" "X17" "X79" "X43" "X24" "X70" "X58" "X6" 
#> 
#> $scores
#> [1] 0.0008005 0.0005960 0.0005430 0.0005314 0.0005176 0.0005028 0.0005064
#> [8] 0.0004638
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
