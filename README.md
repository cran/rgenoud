[![R-CMD-check](https://github.com/JasjeetSekhon/rgenoud/actions/workflows/check-noncontainerized.yaml/badge.svg)](https://github.com/JasjeetSekhon/rgenoud/actions/workflows/check-noncontainerized.yaml)

## rgenoud: R-GENetic Optimization Using Derivatives (RGENOUD)

Walter R. Mebane, Jr. and Jasjeet S. Sekhon

## Introduction

Genoud is a function that combines evolutionary search algorithms with
derivative-based (Newton or quasi-Newton) methods to solve difficult
optimization problems. Genoud may also be used for optimization
problems for which derivatives do not exist.

## How to install

A version is on CRAN. The latest development version can be installed directly from Github
using [devtools](https://github.com/r-lib/devtools). 

```R
if (!require("devtools")) install.packages("devtools")
devtools::install_github("JasjeetSekhon/rgenoud")
```

The package contains compiled code, and you must have a development
environment to install the development version. (Use
`devtools::has_devel()` to check whether you do.) If no development
environment exists, Windows users download and install Rtools and macOS
users download and install Xcode.
