# Hierarchical bootstrap indices

Generates hierarchical bootstrap indices.

## Usage

``` r
hbootindex(groups, strata, B = 199)
```

## Arguments

- groups:

  group membership vector.

- strata:

  strata, optional.

- B:

  number of bootstrap iterations.

## Details

Resampling with replacement with weights proportional to the number of
observations in each of the group level (unique values in `groups`).

Values of `groups` within levels (unique values) of `strata` are
resampled independently of other `strata` levels.

## Value

A matrix with bootstrapped indices, number of columns is `B` + 1. The
column is a resample without replacement (random subsets can be selected
without further reshuffling). Other elements contain indices according
to rules described in Details section (these also randomly reshuffled).

## Author

Peter Solymos

## Examples

``` r
## equal group sizes
groups <- rep(1:4, each=5)
strata <- rep(1:2, each=10)
hbootindex(groups, strata, 3)
#>       [,1] [,2] [,3] [,4]
#>  [1,]    1   13    6   14
#>  [2,]   19   18   16   18
#>  [3,]   18    4   11   17
#>  [4,]    5   13    4   19
#>  [5,]    6   18   20    9
#>  [6,]   11    4    1   10
#>  [7,]   17    6   18   10
#>  [8,]   16    9   18   15
#>  [9,]    4    9    1    5
#> [10,]   12   19   15   15
#> [11,]   13   13   14   10
#> [12,]   20   14    4   16
#> [13,]   10    1    2    7
#> [14,]    2    3    7   10
#> [15,]   14    4    1   17
#> [16,]    7   17   18    9
#> [17,]    9    4   18   19
#> [18,]    8   18    6   16
#> [19,]    3   13    8    2
#> [20,]   15    6   19    1

## unequal group sizes
groups <- groups[-c(5,9,10,11)]
strata <- strata[-c(5,9,10,11)]
hbootindex(groups, strata, 3)
#>       [,1] [,2] [,3] [,4]
#>  [1,]   16    3    6   16
#>  [2,]    3    3    7    1
#>  [3,]    5    4    9   16
#>  [4,]    2    4    3    2
#>  [5,]    9    7   14    6
#>  [6,]   14   15   15    1
#>  [7,]   13    1   12    9
#>  [8,]    7   14    8   12
#>  [9,]    1   10   16    2
#> [10,]   11   12    8    5
#> [11,]   10   10    4   11
#> [12,]    4   14    2   10
#> [13,]   12    9    7    6
#> [14,]    6    4   16   13
#> [15,]   15    9    4   12
#> [16,]    8   16   15   11
```
