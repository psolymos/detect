# AUC ROC plot for fitted models

Area under the receiver-operator (ROC) curve (AUC), and ROC plot methods
for fitted models.

## Usage

``` r
AUC(object, ...)
rocplot(x, ...)
```

## Arguments

- object, x:

  a fitted model object

- ...:

  other arguments

## Value

`AUC` returns AUC value for a model, or a data frame with values for
more models.

`rocplot` returns the values used for the plot invisibly, and as a side
effect it draws a graph.

## Author

Peter Solymos and Monica Moreno
