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
#>  [1,]   19   14   19   17
#>  [2,]    6   19    1   12
#>  [3,]    7    6   15   14
#>  [4,]   18    4   18   16
#>  [5,]    1   17    7   17
#>  [6,]    3   13    4    7
#>  [7,]    2    9   11    5
#>  [8,]    5    9    8    6
#>  [9,]   20   19   18   10
#> [10,]   11    1    6   15
#> [11,]   13   13   14   10
#> [12,]   12   14   12   19
#> [13,]    4   11    1    1
#> [14,]    8    2    1    9
#> [15,]    9    6   20   15
#> [16,]   16    5    9    9
#> [17,]   14   16    6   18
#> [18,]   17   10   18   19
#> [19,]   15    4   16   10
#> [20,]   10   18    4   10

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
