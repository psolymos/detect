# Conversion between truncated and unlimited effective detection distance (EDR)

Conversion between truncated and unlimited effective detection distance
(EDR).

## Usage

``` r
convertEDR(edr, r, truncated=FALSE)
```

## Arguments

- edr:

  effective detection distance. In same units as `r`.

- r:

  truncation distance (radius of point count). In same units as `edr`.

- truncated:

  logical, see Details.

## Details

`truncated = FALSE` means that `edr` is unlimited EDR, and the function
returns the truncated EDR given `r`.

`truncated = TRUE` means that `edr` is truncated EDR given `r`, and the
function returns the unlimited EDR.

## Value

A numeric vector with converted EDR values.

## References

Matsuoka, S. M., Bayne, E. M., Solymos, P., Fontaine, P., Cumming, S.
G., Schmiegelow, F. K. A., & Song, S. A., 2012. Using binomial
distance-sampling models to estimate the effective detection radius of
point-counts surveys across boreal Canada. *Auk*, **129**, 268–282.
\<doi:10.1525/auk.2012.11190\>

Solymos, P., Matsuoka, S. M., Bayne, E. M., Lele, S. R., Fontaine, P.,
Cumming, S. G., Stralberg, D., Schmiegelow, F. K. A. & Song, S. J.,
2013. Calibrating indices of avian density from non-standardized survey
data: making the most of a messy situation. *Methods in Ecology and
Evolution*, **4**, 1047–1058. \<doi:10.1111/2041-210X.12106\>

Supporting info, including a tutorial for the above paper:
<https://github.com/psolymos/QPAD/tree/master/inst/doc/v2>

## Author

Peter Solymos

## Examples

``` r
convertEDR(1, 0.5, truncated=FALSE)
#> [1] 0.4703182
## should be close to 1
convertEDR(convertEDR(1, 0.5, truncated=FALSE), 0.5, truncated=TRUE)
#> [1] 1.000002
```
