# Ovenbird abundances

Ovenbird abundances from BBS

## Usage

``` r
data(oven)
```

## Format

A data frame with 891 observations on the following 11 variables.

- `count`:

  observations

- `route`:

  route id

- `stop`:

  stop id within route

- `pforest`:

  proportion of forest

- `pdecid`:

  proportion of deciduous forest

- `pagri`:

  proportion of agricultural areas

- `long`:

  longitude

- `lat`:

  latitude

- `observ`:

  observer, a factor with levels `ARS` `DW` `RDW` `SVW`

- `julian`:

  Julian day

- `timeday`:

  time of day

## Source

BBS, Erin Bayne (Univ. Alberta), unpublished data set used in Solymos et
al. 2012.

## References

Solymos, P., Lele, S. R. and Bayne, E. 2012. Conditional likelihood
approach for analyzing single visit abundance survey data in the
presence of zero inflation and detection error. *Environmetrics*,
**23**, 197–205. \<doi:10.1002/env.1149\>

## Examples

``` r
data(oven)
str(oven)
#> 'data.frame':    891 obs. of  11 variables:
#>  $ count  : int  1 0 0 1 0 0 0 0 0 0 ...
#>  $ route  : int  2 2 2 2 2 2 2 2 2 2 ...
#>  $ stop   : int  2 4 6 8 10 12 14 16 18 20 ...
#>  $ pforest: num  0.947 0.903 0.814 0.89 0.542 ...
#>  $ pdecid : num  0.575 0.562 0.549 0.679 0.344 ...
#>  $ pagri  : num  0 0 0 0 0.414 ...
#>  $ long   : num  609343 608556 607738 607680 607944 ...
#>  $ lat    : num  5949071 5947735 5946301 5944720 5943088 ...
#>  $ observ : Factor w/ 4 levels "ARS","DW","RDW",..: 4 4 4 4 4 4 4 4 4 4 ...
#>  $ julian : int  181 181 181 181 181 181 181 181 181 181 ...
#>  $ timeday: int  2 4 6 8 10 12 14 16 18 20 ...
```
