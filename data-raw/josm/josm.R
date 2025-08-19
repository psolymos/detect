load("data-raw/josm/josm-data.rda")

f <- function(x, from, to) {
    levels(x)[levels(x) == from] <- to
    x
}

# BADO = Barred Owl=BDOW
spp_from <- "BADO"
spp_to <- "BDOW"
rownames(josm$species)[rownames(josm$species) == spp_from] <- spp_to
josm$species$SpeciesID <- f(josm$species$SpeciesID, spp_from, spp_to)
josm$counts$SpeciesID <- f(josm$counts$SpeciesID, spp_from, spp_to)

# "CANG" Canada Goose=CAGO
spp_from <- "CANG"
spp_to <- "CAGO"
rownames(josm$species)[rownames(josm$species) == spp_from] <- spp_to
josm$species$SpeciesID <- f(josm$species$SpeciesID, spp_from, spp_to)
josm$counts$SpeciesID <- f(josm$counts$SpeciesID, spp_from, spp_to)

# "GRAJ" Canada Jay=CAJA
spp_from <- "GRAY"
spp_to <- "CAJA"
rownames(josm$species)[rownames(josm$species) == spp_from] <- spp_to
josm$species$SpeciesID <- f(josm$species$SpeciesID, spp_from, spp_to)
josm$counts$SpeciesID <- f(josm$counts$SpeciesID, spp_from, spp_to)

# "YEWA" Yellow Warbler=YWAR
spp_from <- "YEWA"
spp_to <- "YWAR"
rownames(josm$species)[rownames(josm$species) == spp_from] <- spp_to
josm$species$SpeciesID <- f(josm$species$SpeciesID, spp_from, spp_to)
josm$counts$SpeciesID <- f(josm$counts$SpeciesID, spp_from, spp_to)

# "MYWA" Myrtle Warbler - considered Yellow-rumped Warbler=YRWA should we lump these 2?
spp_from <- "MYWA"
spp_to <- "YRWA"
josm$species <- josm$species[rownames(josm$species) != spp_from,]
josm$counts$SpeciesID <- f(josm$counts$SpeciesID, spp_from, spp_to)

save(josm, file="data/josm.rda", compress = "xz")

# check species between paired and josm
load("data/paired.rda")
load("data/josm.rda")
load("data-raw/josm/lhreg_data.rda")
lhreg_data$spp <- f(lhreg_data$spp, "GRAJ", "CAJA")

library(mefa4)

compare_sets(josm$species$SpeciesID, paired$SPECIES)
compare_sets(josm$species$SpeciesName, paired$SpeciesName)
setdiff(josm$species$SpeciesID, paired$SPECIES)

compare_sets(josm$species$SpeciesID, lhreg_data$spp)
setdiff(josm$species$SpeciesID, lhreg_data$spp)
compare_sets(paired$SPECIES, lhreg_data$spp)
setdiff(paired$SPECIES, lhreg_data$spp)


sort(setdiff(josm$species$SpeciesName, paired$SpeciesName))
sort(setdiff(paired$SpeciesName, josm$species$SpeciesName))


sort(setdiff(josm$species$SpeciesName, lhreg_data$common_name))
sort(setdiff(paired$SpeciesName, lhreg_data$common_name))

lh <- lhreg_data[,c("spp", "scientific_name", "common_name", "mass", 
    "MaxFreqkHz", "Migr", "Nesthm", "habitat", "food", "behavior")]
colnames(lh) <- c("id", "scientific_name", "common_name", "mass", 
    "max_freq_kHz", "migratory_status", "nest_height_m", "habitat", "food", "behavior")
head(lhreg_data)
head(lh)
write.csv(lh, row.names=F, file="data-raw/josm/traits.csv")
