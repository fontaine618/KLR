# KLR

Kernel Logistic Regression

## Installation

Install the `devtools` package.

Install this package using:
```R
devtools::install_github("fontaine618/KLR")
```

## Example

Fit the model:
```R
library(KLR)
KLRobj = KLR(y, x, kernel="gaussian", lambda=0.001, sigma2=2.0)
```
Obtain predictions on the original data or on new data:
```R
predict(KLRobj)
predict(KLRobj, newx)
```
Obtain level curves (contours):
```R
contours(KLRobj, dims=1:2, res=100, levels=0.5)
```