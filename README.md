[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19839765.svg)](https://doi.org/10.5281/zenodo.19839765)
![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)
![Version](https://img.shields.io/badge/version-0.1.0-blue.svg)

# CSB 

<!-- badges: start -->
<!-- Add badges later: CRAN, R-CMD-check, DOI, etc. -->
<!-- badges: end -->

## Overview

**CSB** is an R package for modeling bounded data on the unit interval \((0,1)\) using the **Contracted Symmetric Beta (CSB)** distribution.

The package provides tools for:

- Density, distribution, and random generation
- Closed-form moments and shape measures
- Parameter estimation by:
  - Method of Moments
  - Maximum Likelihood
- Model comparison with:
  - Beta distribution
  - Kumaraswamy distribution
- Bootstrap goodness-of-fit procedures
- Reproducible workflows for simulation and applications

The CSB model is a flexible bounded distribution obtained through a random contraction mechanism applied to a symmetric beta baseline.

---

## Installation

### From local source

```r
remotes::install_github("YurIriarte/CSB")
