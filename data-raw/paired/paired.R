if (FALSE) {

library(rJava)
library(tabulapdf)
library(pdftools)

z <- extract_tables("https://www.birdpop.org/docs/misc/Alpha_codes_eng.pdf",
    col_names = FALSE)
zz <- do.call(rbind, z)
# '+' before English name indicates a non-species taxon
colnames(zz) <- c("is_species", "english_name", "code_4_letter", "scientific_name", "code_6_letter")
table(zz$is_species,useNA="a")
zz$is_species <- is.na(zz$is_species)
zz <- zz[,c("code_4_letter", "english_name", "code_6_letter", "scientific_name", "is_species")]
zz <- zz[order(zz$code_4_letter),]

zz$not_1st_order_4_letter <- nchar(zz$code_4_letter) > 4
zz$not_1st_order_6_letter <- nchar(zz$code_6_letter) > 6

zz$code_4_letter <- gsub("\\*", "", zz$code_4_letter)
zz$code_6_letter <- gsub("\\*", "", zz$code_6_letter)

write.csv(zz, "data-raw/aou-alpha-codes-eng.csv", row.names = FALSE)

}

x <- read.csv("data-raw/paired/aru-join-table_2016-AUG-22.csv")
str(x)

x$TimeInterval[x$TimeInterval == ""] <- "UNK"
x$Interval[x$Interval == "0-1 min"] <- "0-3 min"
x$Interval[x$Interval == "2-3 min"] <- "0-3 min"
x$Interval[x$Interval == "8-9 min"] <- "5-10 min"
table(x$TimeInterval,x$Interval)

x$SurveyDate <- strptime(x$SurveyDate, "%d-%b-%y", tz = "Canada/Central")
x$Time <- strptime(paste0(x$SurveyDate, " ", x$Time)[1], "%Y-%m-%d %I:%M:%S %p", tz = "Canada/Central")

s <- read.csv("data-raw/paired/aou-codes.csv")

# "BADO" ???
# "BCFR" ???
# "WOFR" ???
# "WOSP" ???
# "CANG" Canada Goose=CAGO
x$SPECIES[x$SPECIES == "CANG"] <- "CAGO"
# "GRAJ" Canada Jay=CAJA
x$SPECIES[x$SPECIES == "GRAJ"] <- "CAJA"
# "YEWA" Yellow Warbler=YWAR
x$SPECIES[x$SPECIES == "YEWA"] <- "YWAR"
# "MYWA" Myrtle Warbler - considered Yellow-rumped Warbler=YRWA should we lump these 2?
x$SPECIES[x$SPECIES == "MYWA"] <- "YRWA"
# "UNBI" "UNBL" "UNDU" "UNHA" "UNKN" "UNSP" "UNWO" - will keep these as NA
x$SPECIES[x$SPECIES %in% c("UNBI", "UNBL", "UNDU", "UNHA", "UNKN", "UNSP", "UNWO")] <- "UNKN"
# "RESQ" Red Squirrel - not a bird

sort(setdiff(unique(x$SPECIES), s$id))

x$SpeciesName <- s$name[match(x$SPECIES, s$id)]

colnames(x)[colnames(x) == "YEAR"] <- "FIREYEAR"
paired <- x
save(paired, file="data/paired.rda", compress = "xz")
