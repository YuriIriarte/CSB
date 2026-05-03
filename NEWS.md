# CSB 0.2.0

## New features
- Added `fitCSB_mle_fixedq()` for maximum likelihood estimation of the CSB distribution when the contraction parameter `q` is fixed.
- Added support for the canonical case `q = 1` through `fitCSB_mle_fixedq(x, q = 1)`.

## Improvements
- Updated scripts in `inst/paper` for reproducibility of the manuscript analyses.
- Improved organization of computational routines related to fixed-parameter estimation.

# CSB 0.1.2

## Improvements
- Added CITATION file for automatic `citation()` support.
- Minor documentation improvements.

# CSB 0.1.0

## Initial release
- Added density function `dCSB()`.
- Added distribution function `pCSB()`.
- Added random generation `rCSB()`.
- Added moment functions:
  - `mCSB()`
  - `meanCSB()`
  - `varCSB()`
  - `cvCSB()`
  - `skewCSB()`
  - `kurtCSB()`
- Added parameter estimation:
  - `fitCSB_mom()`
  - `fitCSB_mle()`
- Added model comparison:
  - `compareCSBmodels()`
- Added bootstrap goodness-of-fit:
  - `gofCSB_boot()`
- Added S3 summary methods.
- Added package documentation.
