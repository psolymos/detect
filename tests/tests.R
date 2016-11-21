#devtools::install_github("psolymos/detect")
library(detect)

## --- run examples with \dontrun sections ---

help_pages <- c("bootstrap", "cmulti", "convertEDR", "databu", "datocc",
"hbootindex", "oven", "svabu", "svocc")

for (i in help_pages) {
    cat("\n\n---------- detect example:", i, "----------\n\n")
    eval(parse(text=paste0("example(", i,
        ", package = 'detect', run.dontrun = TRUE)")))
}
