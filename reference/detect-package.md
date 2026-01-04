# Analyzing Wildlife Data with Detection Error

Models for analyzing site occupancy and count data models with detection
error, including single-visit based models, conditional distance
sampling and time-removal models. Package development was supported by
the Alberta Biodiversity Monitoring Institute and the Boreal Avian
Modelling Project.

## Details

[`svocc`](https://peter.solymos.org/detec/reference/svocc.md): single
visit occupancy model (Lele et al. 2012, Moreno et al. 2010).

[`svabu`](https://peter.solymos.org/detec/reference/svabu.md): single
visit abundance model based on conditional maximum likelihood (Solymos
et al. 2012, Solymos and Lele 2016, Denes et al. 2016).

[`cmulti`](https://peter.solymos.org/detec/reference/cmulti.md):
conditional multinomial maximum likelihood estimation for removal and
(point count) distance sampling, efficient and flexible setup for
varying methodologies (Solymos et al. 2013, Solymos et al. 2018).

## Author

Peter Solymos, Monica Moreno, Subhash R Lele

Maintainer: Peter Solymos \<solymos@ualberta.ca\>

## References

Denes, F., Solymos, P., Lele, S. R., Silveira, L. & Beissinger, S. 2017.
Biome scale signatures of land use change on raptor abundance: insights
from single-visit detection-based models. *Journal of Applied Ecology*,
**54**, 1268–1278. \<doi:10.1111/1365-2664.12818\>

Lele, S.R., Moreno, M. and Bayne, E. 2012. Dealing with detection error
in site occupancy surveys: What can we do with a single survey? *Journal
of Plant Ecology*, **5(1)**, 22–31. \<doi:10.1093/jpe/rtr042\>

Moreno, M. and Lele, S. R. 2010. Improved estimation of site occupancy
using penalized likelihood. *Ecology*, **91**, 341–346.
\<doi:10.1890/09-1073.1\>

Solymos, P., Lele, S. R. and Bayne, E. 2012. Conditional likelihood
approach for analyzing single visit abundance survey data in the
presence of zero inflation and detection error. *Environmetrics*,
**23**, 197–205. \<doi:10.1002/env.1149\>

Solymos, P., Matsuoka, S. M., Bayne, E. M., Lele, S. R., Fontaine, P.,
Cumming, S. G., Stralberg, D., Schmiegelow, F. K. A. & Song, S. J.,
2013. Calibrating indices of avian density from non-standardized survey
data: making the most of a messy situation. *Methods in Ecology and
Evolution*, **4**, 1047–1058. \<doi:10.1111/2041-210X.12106\>

Solymos, P., Lele, S. R. 2016. Revisiting resource selection probability
functions and single-visit methods: clarification and extensions.
*Methods in Ecology and Evolution*, **7**, 196–205.
\<doi:10.1111/2041-210X.12432\>

Solymos, P., Matsuoka, S. M., Cumming, S. G., Stralberg, D., Fontaine,
P., Schmiegelow, F. K. A., Song, S. J., and Bayne, E. M., 2018.
Evaluating time-removal models for estimating availability of boreal
birds during point-count surveys: sample size requirements and model
complexity. *Condor*, **120**, 765–786. \<doi:10.1650/CONDOR-18-32.1\>

Supporting info, including a tutorial for the QPAD method:
<https://github.com/psolymos/QPAD/tree/master/inst/doc/v2>
