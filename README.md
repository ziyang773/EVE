# Equivariant Variance Estimation for Multiple Change-point Model

Software accompanying the paper *“Equivariant Variance Estimation for Multiple Change-point Model”* by  Hao, N., Niu, Y. S., Xiao, H. (2023). [arXiv:2108.09431] (https://arxiv.org/pdf/2108.09431).

---

## Overview

This repository provides the R implementation of the **Equivariant Variance Estimator (EVE)** proposed in the paper. The package includes:

- `eve()` — main function for estimating the variance under frequent mean shifts.
- Automatic, data-driven selection of the tuning parameter \(K\).
- Comparison with other variance estimators such as MS, MAD, DK, and Rice.
- Numerical analysis scripts used to reproduce the tables and figures in the paper.

---

## File Structure

| File / Folder                              | Description                                                    |
|-------------------------------------------|----------------------------------------------------------------|
| `R/core.R`                                | Core implementation of the EVE estimator.                      |
| `Numerical Analysis/Numerical Analysis.R` | Scripts to reproduce simulation tables and figures.            |
| `Numerical Analysis/trios.Rdata`          | SNP genotyping data set used in the numerical study.           |
| `Numerical Analysis/lp1.csv`              | Labor productivity data for Table 6 and Figure 1.              |
| `Numerical Analysis/lp2.csv`              | Labor productivity data for Table 6 and Figure 1.              |

---

## Installation and Usage

```r
# 1. Install the package from GitHub
install.packages("devtools")          # if not already installed
devtools::install_github("ziyang773/EVE")

# 2. Load the package
library(EVE)

# 3. Example usage
set.seed(1234)
n  <- 1000
mu <- rep(c(rep(5, 50), rep(0, 50)), length.out = n)
x  <- rnorm(n) + mu   # true standard deviation is 1

# EVE with automatic selection of K
eve(x)

# EVE with user–specified K
eve(x, K = 5)

# Compare several choices of K
eve(x, K = c(5, 10, 15, 20))

```r
---

## Contact
For questions, comments, or collaborations, please contact:

- **Ning Hao** — <nhao@math.arizona.edu>
- **Ziyang Liu** — <ziyang773@arizona.edu>

---
