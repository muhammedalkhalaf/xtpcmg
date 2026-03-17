# xtpcmg

**Panel Cointegrating Polynomial Regressions via FM-OLS**

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/xtpcmg)](https://CRAN.R-project.org/package=xtpcmg)
<!-- badges: end -->

## Overview

`xtpcmg` estimates panel cointegrating polynomial regressions using:

- **Group-Mean FM-OLS** (Wagner & Reichold, 2023): panel-heterogeneous
  estimator allowing individual polynomial slopes.
- **Pooled FM-OLS** (de Jong & Wagner, 2022): common slope estimator.

Both support quadratic (EKC-type) and cubic polynomial forms with optional
additional I(1) controls and HAC long-run variance via Andrews (1991)
automatic bandwidth selection.

## Installation

```r
install.packages("xtpcmg")
```

## Usage

```r
library(xtpcmg)
dat <- grunfeld_cmg()

# Group-Mean FM-OLS (quadratic)
res <- xtpcmg(dat, y = "invest", x = "mvalue",
              panel_id = "firm", time_id = "year",
              model = "mg", q = 2L)
print(res)

# Pooled FM-OLS
res2 <- xtpcmg(dat, y = "invest", x = "mvalue",
               panel_id = "firm", time_id = "year",
               model = "pmg", q = 2L)
print(res2)
```

## References

- Wagner, M. & Reichold, K. (2023). Econometric Reviews, 42(9–10), 782–827.
  <https://doi.org/10.1080/07474938.2022.2070522>
- de Jong, R.M. & Wagner, M. (2022). Ann. Appl. Stat., 16(1), 416–442.
  <https://doi.org/10.1214/21-AOAS1536>
- Andrews, D.W.K. (1991). Econometrica, 59(3), 817–858.
  <https://doi.org/10.2307/2938229>

## Author

Muhammad Alkhalaf <muhammedalkhalaf@gmail.com>
