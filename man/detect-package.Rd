\name{detect-package}
\alias{detect-package}
\alias{detect}
\docType{package}
\title{
Analyzing Wildlife Data with Detection Error
}
\description{
Models for analyzing site occupancy and count data models
  with detection error, including single-visit based models,
  conditional distance sampling and time-removal models.
  Package development was supported by the
  Alberta Biodiversity Monitoring Institute
  (\url{https://abmi.ca})
  and the Boreal Avian Modelling Project
  (\url{https://borealbirds.ualberta.ca/}).
}
\details{
\code{\link{svocc}}: single visit occupancy model
(Lele et al. 2011, Moreno et al. 2010).

\code{\link{svabu}}: single visit abundance model
based on conditional maximum likelihood
(Solymos et al. 2012, Solymos and Lele 2016, Denes et al. 2016).

\code{\link{cmulti}}: conditional multinomial maximum likelihood
estimation for removal and (point count) distance sampling,
efficient and flexible setup for varying methodologies
(Solymos et al. 2013, Solymos et al. 2018).
}
\author{
Peter Solymos, Monica Moreno, Subhash R Lele

Maintainer: Peter Solymos <solymos@ualberta.ca>
}
\references{
Denes, F., Solymos, P., Lele, S. R., Silveira, L. & Beissinger, S. 2017.
Biome scale signatures of land use change on raptor abundance:
insights from single-visit detection-based models.
\emph{Journal of Applied Ecology}, \bold{54}, 1268--1278.
<doi:10.1111/1365-2664.12818>

Lele, S.R., Moreno, M. and Bayne, E. 2011.
Dealing with detection error in site occupancy surveys:
What can we do with a single survey?
\emph{Journal of Plant Ecology}, \bold{5(1)}, 22--31.
<doi:10.1093/jpe/rtr042>

Moreno, M. and Lele, S. R. 2010.
Improved estimation of site occupancy using penalized likelihood.
\emph{Ecology}, \bold{91}, 341--346.
<doi:10.1890/09-1073.1>

Solymos, P., Lele, S. R. and Bayne, E. 2012.
Conditional likelihood approach for analyzing single visit
abundance survey data in the presence of zero inflation and
detection error.
\emph{Environmetrics}, \bold{23}, 197--205.
<doi:10.1002/env.1149>

Solymos, P., Matsuoka, S. M., Bayne, E. M., Lele, S. R., Fontaine, P.,
Cumming, S. G., Stralberg, D., Schmiegelow, F. K. A. & Song, S. J., 2013.
Calibrating indices of avian density from non-standardized survey data:
making the most of a messy situation.
\emph{Methods in Ecology and Evolution}, \bold{4}, 1047--1058.
<doi:10.1111/2041-210X.12106>

Solymos, P., Lele, S. R. 2016.
Revisiting resource selection probability functions and single-visit methods:
clarification and extensions.
\emph{Methods in Ecology and Evolution}, \bold{7}, 196--205.
<doi:10.1111/2041-210X.12432>

Solymos, P., Matsuoka, S. M., Cumming, S. G., Stralberg, D., Fontaine, P.,
Schmiegelow, F. K. A., Song, S. J., and Bayne, E. M., 2018.
Evaluating time-removal models for estimating availability of boreal birds
during point-count surveys: sample size requirements and model complexity.
\emph{Condor}, \bold{120}, 765--786.
<doi:10.1650/CONDOR-18-32.1>

Supporting info, including a tutorial for the QPAD method:
\url{https://github.com/psolymos/QPAD/tree/master/inst/doc/v2}
}
\keyword{ package }
