\name{paired}
\alias{paired}
\docType{data}
\title{
Example point count data set with paired human and ARU data
}
\description{
Data set from Van Wilgenburg et al. 2017 <doi:10.5751/ACE-00975-120113>.
}
\usage{data("paired")}
\format{
  A data frame with 11340 observations on the following 37 variables.
  \describe{
    \item{\code{sort}}{A numeric vector, sorting ID.}
    \item{\code{UniqueID}}{Character, unique location ID.}
    \item{\code{SurveyDate}}{Date.}
    \item{\code{Time}}{Time.}
    \item{\code{Visit}}{A numeric vector, visit to same location.}
    \item{\code{Observer}}{Character with the observers' initials.}
    \item{\code{SurveyType}}{Character with values
      \code{ARU} (automated recording unit) \code{HUM} (human observer)}
    \item{\code{NoiseLevel}}{Character with values \code{Heavy} \code{Light} \code{Moderate} \code{None} \code{Unusable}.}
    \item{\code{SPECIES}}{Character with values for species,
      4-letter codes following AOU alpha codes.}
    \item{\code{Count}}{A numeric vector, number of individuals counted.}
    \item{\code{TimeInterval}}{Character, original interval the
        individual was detected in.}
    \item{\code{DISTANCE}}{Character with values \code{>100 m}, \code{0-49 m}, \code{50-100 m}, and \code{ARU}.}
    \item{\code{FID}}{A numeric vector for FID.}
    \item{\code{Strata}}{Character, strata.}
    \item{\code{Plot}}{Character, plot.}
    \item{\code{Station}}{Character, station.}
    \item{\code{Latitude}}{A numeric vector, latitude in decimal degrees.}
    \item{\code{Longitude}}{A numeric vector, longitude in decimal degrees.}
    \item{\code{Join_Count}}{A numeric vector.}
    \item{\code{YearLoss}}{A numeric vector.}
    \item{\code{Class_Name}}{Character, the \code{ModisLCC} description.}
    \item{\code{ModisLCC}}{A numeric vector, Modis land cover class.}
    \item{\code{FIRENAME}}{Character, fire name.}
    \item{\code{FIREYEAR}}{A numeric vector, year of fire.}
    \item{\code{Disturbance}}{Character with levels \code{Cutblock}, \code{Fire}, and \code{Undisturbed}}
    \item{\code{PKEY}}{Character, point count ID.}
    \item{\code{PKEYm}}{Character, \code{PKEY} with abbreviated \code{SurveyType}.}
    \item{\code{Noise}}{A numeric (ordinal) version of \code{NoiseLevel}.}
    \item{\code{JULIAN}}{A numeric vector, ordinal day of the year.}
    \item{\code{JDAY}}{A numeric vector, normalized version of \code{JULIAN}.}
    \item{\code{srise}}{A numeric vector, sunrise time in hours.}
    \item{\code{start_time}}{A numeric vector, survey start time in hours.}
    \item{\code{TSSR}}{A numeric vector, normalized time since sunrise (survey time and sunrise time difference as a fraction of a day).}
    \item{\code{Interval}}{Character with levels
      \code{0-3 min}, \code{3-5 min}, \code{5-10 min}, and \code{UNK}.}
    \item{\code{SS}}{Character, station ID.}
    \item{\code{RandomSel}}{A numeric vector for training (1) and hold-out (0) observations.}
    \item{\code{SpeciesName}}{Species common name.}
  }
}
\source{
Van Wilgenburg et al. 2016.
}
\references{
Van Wilgenburg, S. L., P. Solymos, K. J. Kardynal, and M. D. Frey. 2017.
Paired sampling standardizes point count data from humans and acoustic recorders. 
\emph{Avian Conservation and Ecology} \bold{12(1)}:13.
<doi:10.5751/ACE-00975-120113>
}
\examples{
data(paired)
str(paired)
}
\keyword{datasets}
