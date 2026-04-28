[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19839765.svg)](https://doi.org/10.5281/zenodo.19839765)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

# CSB 

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
- Reproducible workflows for simulations and applications

The CSB model is a flexible bounded distribution obtained through a random contraction mechanism applied to a symmetric beta baseline.

---

## Statistical Background

Let

$$
T = X U^{-1/q},
$$

where:

- $X \sim \text{Beta}(\alpha,\alpha)$
- $U \sim \text{Uniform}(1,2^q)$
- $\alpha > 0$, $q > 0$
- $X$ and $U$ are independent

Then $T \in (0,1)$ follows the Contracted Symmetric Beta distribution.

This construction introduces an additional shape parameter $q$, controlling the contraction intensity and generating a wide variety of bounded densities.

---

## Installation

### Option 1: Install from GitHub (development version)

```r
install.packages("remotes")
remotes::install_github("YurIriarte/CSB")
```

### Option 2: Install from Zenodo (archived release)

Download the latest source package from Zenodo:

👉 https://doi.org/10.5281/zenodo.19839765

Then install locally in R:

```r
install.packages("CSB_0.1.1.tar.gz", repos = NULL, type = "source")
```
