# Version 0.5-0, May 21, 2024

* Use `Authors@R` field in `DESCRIPTION`.
* Updated documentation.
* Added SQPAD approach: see `?sqpad`.

# Version 0.4-6, March 8, 2023

* CRAN alert fixed: internal `update.formula.svisit` renamed to `update_formula_svisit` (#9).
* Maintainer email changed to personal email.

# Version 0.4-5, November 7, 2022

* Documentation updates (#8).

# Version 0.4-4, August 12, 2020

* Using fully specified URLs in `DESCRIPTION` as per CRAN request.


# Version 0.4-3, August 11, 2020

* Namespace issue fixed for `scocc.fit` and `extractMLE` (#6).
* Adding `na.rm=TRUE` to `cmulti` to avoid issues with `optimize` interval.

# Version 0.4-2, August 29, 2018

* Updated documentation, references.

# Version 0.4-1, January 30, 2018

* Denes et al. 2017 paper added to documentation.
* `inst/ChangeLog` is removed and replaced by `NEWS.md`.
* `cmulti` models now have a `predict` method, see examples.
* `"fmix"`" model type with time varying rates added to `cmulti` models.
* Added new data set: `paired`.

# Version 0.4-0, March 2, 2016

* `load_BAM_QPAD` function is deprecated (migrated to **QPAD** package).
* `hbootindex` bugfix: `sample()` behaviour when
  `length(x) == 1` is recognized when resampling from a length 1 vector.
* **pbapply** package is now a dependency.
* Various minor fixes to satisfy R-devel checks.
* Negative Binomial option added to `svabu`.
* `cmulti` returns correct `nobs` (nonzero total counts).

# Version 0.3-2, May 15, 2014

* Bugfix in `cmulti2.fit`: n was undefined (reported by Julie Hart).
* Bugfix in `hbootindex`: n was not set when strata was missing.
* Documentation fixes (comma separated lists of variable names inside `\code`).

# Version 0.3-1, Sept 25, 2013

* `hbootindex`: routine improved, now uses weights proportional to groups sizes.
* Cleanup to satisfy R 3.0.2 check: remove `:::`.
* Bugfix: `svocc` returned k=1 object irrespective of the `n.clones`
  option settings for MLE with `method = "dc"`.
* Help page added for internal objects (e.g. `cmulti.fit0`).
* `drop.scope.svisit` is no longer treated as a method because
  `drop.scope` is not a generic.

# Version 0.3-0, July 29, 2013

* `jags.engine` removed, replaced by **dcmle** dependency.
* `dFormula` and `checkDesign` added with **Formula** dependency,
  this will replace the old `svisitFormula` interface
  to allow handling of more complex designs.
* `logdmultinom`: simplified version of `dmultinom`.
* `bymethod`: collapses (pools) observations by methodology.
* `cmulti`: new function for conditional multinomial estimation.
* `load_BAM_QPAD`: new function to load BAM QPAD parameter
  estimates and support functions.
* `hbootindex`: new function to sample indices with replacement
  for grouped/stratified data.
* `fitted.cmulti` added: get estimated singing rate or distance model
  parameter (not available for `type = "mix"`)
* `predict` methods: response (lhs) is dropped from formula
  when `newdata` is provided.

# Version 0.2-2, May 2, 2012

* `inst/COPYING` file removed (standard license).

# Version 0.2-1, November 25, 2011

* Help pages edited.

# Version 0.2-0, October 18, 2011

* `svabu` and related functions added.
* summary: nobs was returned as NULL, so AIC was not printed. Now fixed.
* R (>= 2.13.0) added to `DESCRIPTION` (reported by Uwe Ligges).

# Version 0.1-0, September 27, 2011

* First public release on CRAN.
