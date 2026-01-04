# Internal functions

Internal functions, not intended for use on their own.

## Usage

``` r
cmulti.fit0(Y, D, type=c("rem", "mix", "dis", "fmix"),
    interval=c(-25, 25), ...)
drop.scope.svisit(terms1, terms2, model = c("sta", "det"))
```

## Arguments

- Y:

  this contains the cell counts. See
  [`cmulti.fit`](https://peter.solymos.org/detec/reference/cmulti.md).

- D:

  design matrix, that describe the interval endpoints for the sampling
  methodology, dimensions must match dimensions of `Y`. See
  [`cmulti.fit`](https://peter.solymos.org/detec/reference/cmulti.md).

- type:

  character, one of `"rem"` (removal sampling, homogeneous singing
  rates), `"mix"` and `"fmix"` (removal sampling, heterogeneous singing
  rates), `"dis"` (distance sampling, half-normal detection function for
  point counts, circular area). See
  [`cmulti.fit`](https://peter.solymos.org/detec/reference/cmulti.md).

- interval:

  the interval used in
  [`optimize`](https://rdrr.io/r/stats/optimize.html).

- terms1:

  the terms or formula for the upper/lower scope. See
  [`drop.scope`](https://rdrr.io/r/stats/factor.scope.html).

- model:

  character, the type of model to be considered.

## Author

Peter Solymos
