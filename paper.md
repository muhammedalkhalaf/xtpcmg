---
title: 'R Packages for Advanced Panel Data Econometrics'
tags:
  - R
  - econometrics
  - panel data
  - quantile regression
  - cross-sectional dependence
  - missing data
  - Granger causality
authors:
  - name: Muhammad Abdullah Alkhalaf
    orcid: 0009-0002-2677-9246
    corresponding: true
    email: muhammedalkhalaf@gmail.com
    affiliation: 1
affiliations:
  - name: Rufyq Elngeh for Academic and Business Services, Riyadh, Saudi Arabia
    index: 1
date: 18 March 2026
bibliography: paper.bib
---

# Summary

Seven R packages address gaps in panel data econometrics: `xtmispanel` for missing data diagnostics and imputation, `xtcsdq` for cross-sectional dependence testing in quantile regressions [@Demetrescu2023], `xtpcaus` for panel Granger causality, `xtpcmg` for panel cointegrating polynomial regressions via FM-OLS, `xtrec` for recursively detrended panel unit root tests [@Westerlund2015], `xtfifevd` for fixed effects with time-invariant regressors, and `xtqsh` for quantile slope homogeneity testing [@Galvao2017]. All are open-source under GPL-3.

# Statement of Need

R provides foundational panel tools through `plm` [@Croissant2008] and `fixest` [@Berge2021], yet several recent methodological advances lack R implementations: specialized panel missing data diagnostics, cross-sectional dependence tests for quantile regressions, panel Fourier Granger causality, cointegrating polynomial regressions with FM-OLS, recursively detrended unit root tests, estimation of time-invariant regressor effects, and quantile slope homogeneity tests.

# Packages

## xtmispanel

Detects, diagnoses, and imputes missing values in panel data. Provides per-variable, per-panel, and per-period summaries, approximate @Little1988 MCAR and MAR tests, and eleven imputation methods (mean, median, LOCF, NOCB, linear/spline interpolation, PMM, hot-deck, KNN, random forest, EM).

```r
library(xtmispanel)
result <- xtmispanel(data, index = c("country", "year"),
                     detect = TRUE, test = TRUE,
                     impute = "pmm", target = "gdp")
```

## xtcsdq

Tests for cross-sectional error dependence in panel quantile regressions following @Demetrescu2023. Provides T_τ and bias-corrected T̃_τ statistics, plus a portmanteau M_K statistic across quantiles.

```r
library(xtcsdq)
result <- xtcsdq(y ~ x1 + x2, data = panel_data,
                 index = c("country", "year"),
                 quantiles = c(0.25, 0.5, 0.75))
summary(result)
```

## xtpcaus

Panel Fourier Toda–Yamamoto (PFTY) and Panel Quantile Causality (PQC) tests with bootstrap p-values.

```r
library(xtpcaus)
result <- xtpcaus(data = panel_data, y = "gdp", x = "investment",
                  panel_id = "country", time_id = "year",
                  test = "pfty", kmax = 3, nboot = 499)
summary(result)
```

## xtpcmg

Group-Mean and Pooled FM-OLS for panel cointegrating polynomial regressions. Group-Mean follows @Wagner2023; Pooled follows @deJong2022. Supports quadratic and cubic polynomials with HAC long-run variance, turning point analysis, and @Swamy1970 slope homogeneity test.

```r
library(xtpcmg)
result <- xtpcmg(data = panel_data, y = "emissions",
                 x = "gdp", panel_id = "country",
                 time_id = "year", model = "mg", q = 2)
summary(result)
```

## xtrec

Recursively detrended panel unit root tests of @Westerlund2015. The basic t-REC assumes iid errors; the robust t-RREC accounts for serial correlation, cross-sectional dependence, and heteroskedasticity.

```r
library(xtrec)
result <- xtrec(data = panel_data, var = "gdp",
                panel_id = "country", time_id = "year",
                robust = TRUE)
summary(result)
```

## xtfifevd

FEVD [@Plumper2007], FEF, and FEF-IV [@PesaranZhou2016] estimators for time-invariant regressors in panel fixed effects models.

```r
library(xtfifevd)
result <- xtfifevd(y ~ x1 + x2 | z1 + z2, data = panel_data,
                   id = "country", time = "year",
                   method = "fef")
summary(result)
```

## xtqsh

Quantile regression slope homogeneity test of @Galvao2017. Provides Ŝ (chi-squared) and D̂ (standard normal) statistics for joint and marginal homogeneity.

```r
library(xtqsh)
result <- xtqsh(y ~ x1 + x2, data = panel_data,
                id = "country", time = "year",
                tau = c(0.25, 0.5, 0.75))
summary(result)
```

# Acknowledgements

The author acknowledges the econometric researchers whose original Stata, GAUSS, and MATLAB implementations provided the foundation for these R packages.

# References
