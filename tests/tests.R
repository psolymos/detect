#devtools::install_github("psolymos/detect")
library(detect)

## --- run examples with \dontrun sections ---

help_pages <- c(#"bootstrap",
    "cmulti", "convertEDR", "databu", "datocc",
    "hbootindex", "oven", "svabu", "svocc")

for (i in help_pages) {
    cat("\n\n---------- detect example:", i, "----------\n\n")
    eval(parse(text=paste0("example(", i,
        ", package = 'detect', run.dontrun = TRUE)")))
}

if (FALSE) {

library(detect)
data(datocc)
m0 <- svocc(W ~ x1 | x1 + x3, datocc,
    method="dc",n.adapt=100,n.update=100,n.iter=100)
extractMLE(m0)

}
