load("data-raw/lhreg/cor_matrix.rda")
load("data-raw/lhreg/lhreg_data.rda")

lhreg <- list(traits = lhreg_data, phylo_cor = cor_matrix)
save(lhreg, file="data/lhreg.rda", compress = "xz")
