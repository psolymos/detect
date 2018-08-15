# detect: analyzing wildlife data with detection error

[![CRAN version](http://www.r-pkg.org/badges/version/detect)](http://cran.rstudio.com/web/packages/detect/index.html)
[![CRAN download stats](http://cranlogs.r-pkg.org/badges/grand-total/detect)](https://www.rdocumentation.org/packages/detect/)
[![Linux build Status](https://travis-ci.org/psolymos/detect.svg?branch=master)](https://travis-ci.org/psolymos/detect)
[![Windows build status](https://ci.appveyor.com/api/projects/status/5y5fwgv90f8i84ck?svg=true)](https://ci.appveyor.com/project/psolymos/detect)
[![Code coverage status](https://codecov.io/gh/psolymos/detect/branch/master/graph/badge.svg)](https://codecov.io/gh/psolymos/detect)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

The R package implements models to analyze
site occupancy and count data models with detection error.
The package development was supported by the
[Alberta Biodiversity Monitoring Institute (ABMI)](http://www.abmi.ca/)
and the [Boreal Avian Modelling (BAM) Project](http://www.borealbirds.ca/).

Main functions:

* `svocc`: single visit occupancy model (Lele et al. 2011, Moreno et al. 2010).
* `svabu`: single visit Poisson and Negative Binomial abundance model based on conditional maximum likelihood (Solymos et al. 2012, Denes et al. 2016, Solymos & Lele 2016).
* `cmulti`: conditional multinomial maximum likelihood estimation for removal and (point count) distance sampling, efficient and flexible setup for varying methodologies (Solymos et al. 2013, Solymos et al. 2018).

## Versions

Install CRAN release version (recommended):

```R
install.packages("detect")
```

Development version:

```R
library(devtools)
install_github("psolymos/detect")
```

User visible changes in the package are listed in the [NEWS](https://github.com/psolymos/detect/blob/master/NEWS.md) file.

Use the [issue tracker](https://github.com/psolymos/detect/issues)
to report a problem.

## References

Denes, F., Solymos, P., Lele, S. R., Silveira, L. & Beissinger, S. 2017.
Biome scale signatures of land use change on raptor abundance:
insights from single-visit detection-based models.
_Journal of Applied Ecology_, **54**, 1268--1278.
[DOI: 10.1111/1365-2664.12818](https://dx.doi.org/10.1111/1365-2664.12818)

Lele, S.R., Moreno, M. and Bayne, E. 2011.
Dealing with detection error in site occupancy surveys:
What can we do with a single survey?
_Journal of Plant Ecology_, **5(1)**, 22--31.
[DOI: 10.1093/jpe/rtr042](https://dx.doi.org/10.1093/jpe/rtr042)

Moreno, M. and Lele, S. R. 2010.
Improved estimation of site occupancy using penalized likelihood.
_Ecology_, **91**, 341--346.
[DOI: 10.1890/09-1073.1](https://dx.doi.org/10.1890/09-1073.1)

Solymos, P., Lele, S. R. and Bayne, E. 2012.
Conditional likelihood approach for analyzing single visit
abundance survey data in the presence of zero inflation and
detection error.
_Environmetrics_, **23**, 197--205.
[DOI: 10.1002/env.1149](https://dx.doi.org/10.1002/env.1149)

Solymos, P., Matsuoka, S. M., Bayne, E. M., Lele, S. R., Fontaine, P.,
Cumming, S. G., Stralberg, D., Schmiegelow, F. K. A. & Song, S. J., 2013.
Calibrating indices of avian density from non-standardized survey data:
making the most of a messy situation.
_Methods in Ecology and Evolution_, **4**, 1047--1058.
[DOI: 10.1111/2041-210X.12106](https://dx.doi.org/10.1111/2041-210X.12106)

Solymos, P., Lele, S. R. 2016.
Revisiting resource selection probability functions and single-visit methods:
clarification and extensions.
_Methods in Ecology and Evolution_, **7**, 196--205.
[DOI: 10.1111/2041-210X.12432](https://dx.doi.org/10.1111/2041-210X.12432)

Solymos, P., Matsuoka, S. M., Cumming, S. G., Stralberg, D., Fontaine, P.,
Schmiegelow, F. K. A., Song, S. J., and Bayne, E. M., 2018.
Evaluating time-removal models for estimating availability of boreal birds
during point-count surveys: sample size requirements and model complexity.
_Condor_, in press.

Supporting info, including a tutorial for the 
[QPAD method](https://github.com/psolymos/QPAD/tree/master/inst/doc/v2).
