# Example point count data set with paired human and ARU data

Data set from Van Wilgenburg et al. 2017
\<doi:10.5751/ACE-00975-120113\>.

## Usage

``` r
data("paired")
```

## Format

A data frame with 11340 observations on the following 37 variables.

- `sort`:

  A numeric vector, sorting ID.

- `UniqueID`:

  Character, unique location ID.

- `SurveyDate`:

  Date.

- `Time`:

  Time.

- `Visit`:

  A numeric vector, visit to same location.

- `Observer`:

  Character with the observers' initials.

- `SurveyType`:

  Character with values `ARU` (automated recording unit) `HUM` (human
  observer)

- `NoiseLevel`:

  Character with values `Heavy` `Light` `Moderate` `None` `Unusable`.

- `SPECIES`:

  Character with values for species, 4-letter codes following AOU alpha
  codes.

- `Count`:

  A numeric vector, number of individuals counted.

- `TimeInterval`:

  Character, original interval the individual was detected in.

- `DISTANCE`:

  Character with values `>100 m`, `0-49 m`, `50-100 m`, and `ARU`.

- `FID`:

  A numeric vector for FID.

- `Strata`:

  Character, strata.

- `Plot`:

  Character, plot.

- `Station`:

  Character, station.

- `Latitude`:

  A numeric vector, latitude in decimal degrees.

- `Longitude`:

  A numeric vector, longitude in decimal degrees.

- `Join_Count`:

  A numeric vector.

- `YearLoss`:

  A numeric vector.

- `Class_Name`:

  Character, the `ModisLCC` description.

- `ModisLCC`:

  A numeric vector, Modis land cover class.

- `FIRENAME`:

  Character, fire name.

- `FIREYEAR`:

  A numeric vector, year of fire.

- `Disturbance`:

  Character with levels `Cutblock`, `Fire`, and `Undisturbed`

- `PKEY`:

  Character, point count ID.

- `PKEYm`:

  Character, `PKEY` with abbreviated `SurveyType`.

- `Noise`:

  A numeric (ordinal) version of `NoiseLevel`.

- `JULIAN`:

  A numeric vector, ordinal day of the year.

- `JDAY`:

  A numeric vector, normalized version of `JULIAN`.

- `srise`:

  A numeric vector, sunrise time in hours.

- `start_time`:

  A numeric vector, survey start time in hours.

- `TSSR`:

  A numeric vector, normalized time since sunrise (survey time and
  sunrise time difference as a fraction of a day).

- `Interval`:

  Character with levels `0-3 min`, `3-5 min`, `5-10 min`, and `UNK`.

- `SS`:

  Character, station ID.

- `RandomSel`:

  A numeric vector for training (1) and hold-out (0) observations.

- `SpeciesName`:

  Species common name.

## Source

Van Wilgenburg et al. 2016.

## References

Van Wilgenburg, S. L., P. Solymos, K. J. Kardynal, and M. D. Frey. 2017.
Paired sampling standardizes point count data from humans and acoustic
recorders. *Avian Conservation and Ecology* **12(1)**:13.
\<doi:10.5751/ACE-00975-120113\>

## Examples

``` r
data(paired)
str(paired)
#> 'data.frame':    11340 obs. of  37 variables:
#>  $ sort        : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ UniqueID    : chr  "05-041-01" "05-041-01" "05-041-01" "05-041-01" ...
#>  $ SurveyDate  : POSIXlt, format: "2015-06-22" "2015-06-22" ...
#>  $ Time        : POSIXlt, format: "2015-06-22 08:29:00" "2015-06-22 08:29:00" ...
#>  $ Visit       : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ Observer    : chr  "BFDO" "BFDO" "BFDO" "BFDO" ...
#>  $ SurveyType  : chr  "ARU" "ARU" "ARU" "ARU" ...
#>  $ NoiseLevel  : chr  "Moderate" "Moderate" "Moderate" "Moderate" ...
#>  $ SPECIES     : chr  "FOSP" "LISP" "SWTH" "YBFL" ...
#>  $ Count       : int  1 1 1 1 1 1 1 1 1 2 ...
#>  $ TimeInterval: chr  "0-1 min" "1-2 min" "0-1 min" "9-10 min" ...
#>  $ DISTANCE    : chr  "ARU" "ARU" "ARU" "ARU" ...
#>  $ FID         : int  44 44 44 44 44 44 44 44 44 44 ...
#>  $ Strata      : chr  "5" "5" "5" "5" ...
#>  $ Plot        : chr  "41" "41" "41" "41" ...
#>  $ Station     : chr  "1" "1" "1" "1" ...
#>  $ Latitude    : num  58.1 58.1 58.1 58.1 58.1 ...
#>  $ Longitude   : num  -109 -109 -109 -109 -109 ...
#>  $ Join_Count  : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ YearLoss    : int  0 0 0 0 0 0 0 0 0 0 ...
#>  $ Class_Name  : chr  "Sparse Coniferous" "Sparse Coniferous" "Sparse Coniferous" "Sparse Coniferous" ...
#>  $ ModisLCC    : int  8 8 8 8 8 8 8 8 8 8 ...
#>  $ FIRENAME    : chr  "ALTA" "ALTA" "ALTA" "ALTA" ...
#>  $ FIREYEAR    : int  1980 1980 1980 1980 1980 1980 1980 1980 1980 1980 ...
#>  $ Disturbance : chr  "Fire" "Fire" "Fire" "Fire" ...
#>  $ PKEY        : chr  "05-041-01_1" "05-041-01_1" "05-041-01_1" "05-041-01_1" ...
#>  $ PKEYm       : chr  "05-041-01_1_A" "05-041-01_1_A" "05-041-01_1_A" "05-041-01_1_A" ...
#>  $ Noise       : int  2 2 2 2 2 2 2 2 2 2 ...
#>  $ JULIAN      : int  172 172 172 172 172 172 172 172 172 172 ...
#>  $ JDAY        : num  0.471 0.471 0.471 0.471 0.471 ...
#>  $ srise       : num  4.22 4.22 4.22 4.22 4.22 ...
#>  $ start_time  : num  8.48 8.48 8.48 8.48 8.48 ...
#>  $ TSSR        : num  0.178 0.178 0.178 0.178 0.178 ...
#>  $ Interval    : chr  "0-3 min" "0-3 min" "0-3 min" "5-10 min" ...
#>  $ SS          : chr  "5.041" "5.041" "5.041" "5.041" ...
#>  $ RandomSel   : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ SpeciesName : chr  "Fox Sparrow" "Lincoln's Sparrow" "Swainson's Thrush" "Yellow-bellied Flycatcher" ...
```
