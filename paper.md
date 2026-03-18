---
title: 'R Packages for Advanced Panel Data Econometrics: Missing Data, Quantile Methods, Causality, and Estimation'
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

We present seven R packages for advanced panel data econometrics: `xtmispanel` for missing data diagnostics and imputation, `xtcsdq` for cross-sectional dependence testing in quantile regressions, `xtpcaus` for panel Granger causality with Fourier terms and quantile methods, `xtpcmg` for panel cointegrating polynomial regressions via FM-OLS, `xtrec` for recursively detrended panel unit root tests, `xtfifevd` for fixed effects estimation with time-invariant regressors, and `xtqsh` for quantile regression slope homogeneity testing. These packages address methodological gaps in the R ecosystem for panel data analysis, implementing recent advances that were previously available only in Stata, GAUSS, or MATLAB.

# Statement of Need

Panel data methods are central to modern empirical economics, political science, and sociology. While R provides foundational panel tools through `plm` [@Croissant2008] and `fixest` [@Berge2021], several important methodological advances lack R implementations:

1. **Missing data in panels** requires specialized diagnostics beyond cross-sectional methods, including per-panel and per-period analysis, MCAR/MAR testing, and panel-aware imputation [@Little1988].
2. **Cross-sectional dependence in quantile regressions** has only recently received formal testing procedures [@Demetrescu2023], with no existing software implementation.
3. **Panel Fourier Granger causality** [@Yilanci2020panel] and **panel quantile causality** [@Wang2022] extend standard panel VAR methods but lack R packages.
4. **Panel cointegrating polynomial regressions** with FM-OLS [@Wagner2023; @deJong2022] enable estimation of nonlinear long-run relationships (e.g., Environmental Kuznets Curves) but have no R implementation.
5. **Recursively detrended panel unit root tests** [@Westerlund2015] offer improved power through novel detrending but are unavailable in R.
6. **Fixed Effects Vector Decomposition** [@Plumper2007] and **Fixed Effects Filtered** [@PesaranZhou2016] methods address the long-standing problem of estimating time-invariant regressor effects in panel models.
7. **Quantile slope homogeneity tests** [@Galvao2017] extend mean-based homogeneity tests to distributional analysis but have no R implementation.

# Packages

## xtmispanel

Provides comprehensive tools for detecting, diagnosing, and imputing missing data in panel datasets. Includes per-variable, per-panel, and per-period summary tables, missing data pattern matrices, approximate @Little1988 MCAR test, logistic regression MAR test, and eleven imputation methods (mean, median, LOCF, NOCB, linear/spline interpolation, PMM, hot-deck, KNN, random forest, and EM).

```r
library(xtmispanel)
diag <- xtmisdiag(data, id = "country", time = "year")
summary(diag)
imputed <- xtmisimpute(data, id = "country", time = "year",
                       method = "pmm")
```

## xtcsdq

Implements tests of no cross-sectional error dependence (CSD) in panel quantile regressions following @Demetrescu2023. Provides the T_τ statistic and its bias-corrected version, plus a portmanteau statistic M_K aggregating evidence across multiple quantile levels. Supports pooled fixed-effects QR, individual unit QR, and post-estimation modes.

```r
library(xtcsdq)
result <- xtcsdq(y ~ x1 + x2, data = panel_data,
                 id = "country", time = "year",
                 taus = c(0.25, 0.5, 0.75))
summary(result)
```

## xtpcaus

Implements two panel Granger causality tests. The Panel Fourier Toda-Yamamoto (PFTY) test extends the @TodaYamamoto1995 framework with Fourier terms following @Yilanci2020panel and @Emirmahmutoglu2011. The Panel Quantile Causality (PQC) test examines causality across the conditional distribution following @Wang2022 and @Chuang2009. Bootstrap p-values account for cross-sectional dependence.

```r
library(xtpcaus)
# Panel Fourier Toda-Yamamoto
result <- pfty(y, x, data = panel_data, id = "country",
               time = "year", max_freq = 3, nboot = 1000)
summary(result)
```

## xtpcmg

Implements Group-Mean and Pooled Fully Modified OLS (FM-OLS) estimators for panel cointegrating polynomial regressions. The Group-Mean FM-OLS follows @Wagner2023, and the Pooled FM-OLS follows @deJong2022. Supports quadratic and cubic polynomial regressors with HAC long-run variance estimation, turning point analysis with delta-method confidence intervals, and @Swamy1970 slope homogeneity test.

```r
library(xtpcmg)
result <- xtpcmg(y ~ x + I(x^2), data = panel_data,
                 id = "country", time = "year",
                 estimator = "group_mean", degree = 2)
summary(result)
# Turning point analysis
tp <- turning_point(result)
```

## xtrec

Implements the recursively detrended panel unit root tests of @Westerlund2015. The basic t-REC test assumes iid errors; the robust t-RREC test accounts for serial correlation, cross-sectional dependence, and heteroskedasticity via defactoring and BIC-selected lag augmentation. Both tests have standard normal null distributions.

```r
library(xtrec)
result <- xtrec(y ~ 1, data = panel_data,
                id = "country", time = "year", robust = TRUE)
summary(result)
```

## xtfifevd

Implements the Fixed Effects Vector Decomposition (FEVD) estimator of @Plumper2007 and the Fixed Effects Filtered (FEF) and FEF-IV estimators of @PesaranZhou2016. Provides consistent standard errors for time-invariant regressors using the corrected variance estimator, with comparison of corrected versus raw standard errors.

```r
library(xtfifevd)
result <- xtfifevd(y ~ x1 + x2 | z1 + z2, data = panel_data,
                   id = "country", time = "year",
                   method = "fef", time_invariant = c("z1", "z2"))
summary(result)
```

## xtqsh

Implements the quantile regression slope homogeneity test of @Galvao2017. Provides Ŝ (chi-squared) and D̂ (standard normal) test statistics for both joint and marginal slope homogeneity. Returns minimum distance QR estimates with bandwidth selection following Bofinger (1975) or Hall and Sheather (1988).

```r
library(xtqsh)
result <- xtqsh(y ~ x1 + x2, data = panel_data,
                id = "country", time = "year",
                taus = c(0.25, 0.5, 0.75))
summary(result)
```

# Acknowledgements

The author acknowledges the econometric researchers whose original methodologies and software implementations (in Stata, GAUSS, and MATLAB) provided the foundation for these R packages.

# References
