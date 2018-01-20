#--- SETTINGS
library("curl")
library("taxize")
library("traits")


dir_prj <- "Prj04_GenetFragment_MetaAnalysis"
dir_in <- file.path(dir_prj, "2_Data", "4_DataExtracted")
dir_out <- file.path(dir_prj, "2_Data", "5_DataCleaned")


# File names
if (!exists("tag_newdate")) tag_newdate <- format(Sys.Date(), "%Y%m%d")
if (!exists("tag_usedate")) tag_usedate <- "20171107"

ftaxonomy <- file.path(dir_in, "20170130_taxonomy.rda")
ftaxonomy_new <- file.path(dir_in, paste0(tag_newdate, "_taxonomy.rda"))
feol <- file.path(dir_in, "20170130_eol-traitbank.rda")
feol_new <- file.path(dir_in, paste0(tag_newdate, "_eol-traitbank.rda"))

# fname_traits_plants <- file.path(dir_in, "SpeciesTraits_Plants16.3.17.xls")
# fname_traits_plants <- file.path(dir_in, "Plant-Traits-10.4.17.xls")
fname_traits_plants <- file.path(dir_in, "20170411_SpeciesTraits_Plants.xlsx")
fname_traits_animals <- file.path(dir_in, "20170411_SpeciesTraits_Animals.xlsx")
fname_age_animals <- file.path(dir_in, "20170411_SpeciesTraits_AnimalsAge.xlsx")
fname_habitatreclass <- file.path(dir_in, "20171105_HabitatReClassification.csv")

fin <- file.path(dir_out, paste0(tag_usedate, "_ExtractedDataCleaned2_wide.rds"))
foutw <- file.path(dir_out, paste0(tag_newdate, "_ExtractedDataCleaned3_wide.rds"))
foutl <- file.path(dir_out, paste0(tag_newdate, "_ExtractedDataCleaned3_long.rds"))


#--- READ DATA
print(paste(Sys.time(), "reading input", shQuote(basename(fin))))
dat0 <- readRDS(fin)

ftemp <- file.path("~", ".eol_API_key")
eol_API_key <- if (file.exists(ftemp)) {
    temp <- try(taxize::getkey(readLines(ftemp), "eolApiKey"))
    if (inherits(temp, "try-error")) NULL else temp
  } else NULL

age_categories <- levels(dat0$FragmentAgeCat_ord_yrs)
stopifnot(identical(age_categories, c("(0, 2] yrs", "(2-10] yrs", "(10-15] yrs",
  "(15-50] yrs", "(50-100] yrs", "(100-200] yrs", ">200 yrs")))

age_categories2 <- levels(dat0$EffectAgeCat2_ord_yrs)
stopifnot(identical(age_categories2, c("(0, 10] yrs", "(10-50] yrs", ">50 yrs")))


#------------------------------------
#--- TAXONOMY
do_taxonomy <- TRUE
test_taxonomy <- NULL

temp <- dat0$Species_accepted
itemp <- is.na(temp)
temp[itemp] <- dat0$Species_resolved[itemp]
specieslist2 <- unique(temp)

if (file.exists(ftaxonomy)) {
  test_taxonomy <- new.env(parent = emptyenv())
  load(ftaxonomy, envir = test_taxonomy)
  do_taxonomy <- !identical(specieslist2, test_taxonomy$specieslist2)
  do_taxonomy <- !all(specieslist2 %in% test_taxonomy$specieslist2)
}

if (!do_taxonomy) {
  ids <- match(specieslist2, test_taxonomy$specieslist2, nomatch = 0)

  if (identical(specieslist2, test_taxonomy$specieslist2[ids])) {
    tax2 <- get("tax2", envir = test_taxonomy)[ids]
  } else {
    do_taxonomy <- TRUE
  }
}

#--- Obtain taxonomy
db_pref <- "ncbi"
tax_ranks <- c("kingdom", "phylum", "class", "order", "family")

if (do_taxonomy &&
  (requireNamespace("curl") && curl::has_internet()) || !requireNamespace("curl")) {

  # Extract taxonomy from NCBI and ITIS
  ntemp <- c(1, 1)
  repeat {
    repeat {
      print(paste(Sys.time(), "attempt", ntemp[1], "to reach ITIS and NCBI"))
      pinged <- taxize::itis_ping() && taxize::ncbi_ping()

      if (pinged) {
        break
      }

      print("------ ITIS and/or NCBI are currently not accessible")

      if (ntemp[1] > 5) {
        stop("Too many attempts failed 1")
      }
      ntemp[1] <- ntemp[1] + 1

      print("------ code waits for 30 s before attempting again")
      Sys.sleep(30)
    }

    print(paste(Sys.time(), "attempt", ntemp[2], "to download data from ITIS and NCBI"))
    tax <- try(taxize::tax_name(query = specieslist2, get = tax_ranks, db = "both", pref = db_pref,
      accepted = TRUE, verbose = FALSE, ask = FALSE))

    if (!inherits(tax, "try-error")) break

    if (ntemp[2] > 5) {
      stop("Too many attempts failed 2")
    }
    ntemp[2] <- ntemp[2] + 1
  }
  print(warnings())

  # Clean taxonomy
  tax$kingdom[tax$kingdom %in% "Viridiplantae"] <- "Plantae"
  tax$kingdom[tax$kingdom %in% "Animalia"] <- "Metazoa"
  #   - our study doesn't consider other classes besides "Magnoliopsida" in phylum "Magnoliopsida"
  tax$phylum[is.na(tax$phylum) & tax$class %in% "Magnoliopsida"] <- "Streptophyta"
  tax$class[is.na(tax$class) & tax$phylum %in% "Streptophyta"] <- "Magnoliopsida"

  nonNA_taxonomy <- function(x) {
    if (NROW(x) > 1) {
      apply(x, 2, function(r) {
          inotna <- !is.na(r)
          if (any(inotna)) r[inotna][1] else NA
        })
    } else {
      unlist(x)
    }
  }

  tax2 <- data.frame(array(NA, dim = c(length(specieslist2), dim(tax)[2]),
    dimnames = list(NULL, c("Species_resolved", dimnames(tax)[[2]][-1]))))
  tax2[, "Species_resolved"] <- unique(dat0$Species_resolved)
  tax2[, "query"] <- specieslist2

  for (k in seq_along(specieslist2)) {
    irows_tax <- which(specieslist2[k] == tax[, "query"])
    temp1 <- NULL

    if (length(irows_tax) > 1) {
      temp <- tax[irows_tax, ]
      i_db_pref <- which(temp[, "db"] == db_pref)

      if (length(i_db_pref) == 1) {
        temp1 <- unlist(temp[i_db_pref, tax_ranks])
        temp2 <- nonNA_taxonomy(temp[-i_db_pref, tax_ranks])
        if (anyNA(temp1)) {
          temp1 <- ifelse(is.na(temp1), temp2, temp1)
        }
      } else {
        temp1 <- nonNA_taxonomy(temp[, tax_ranks])
      }

      tax2[k, tax_ranks] <- temp1
    }
  }

  #Fill in missing information
  id_addtax <- tax2$query == "Dianthus seguieri glaber"
  tax2[id_addtax, tax_ranks] <- c("Plantae", "Streptophyta", "Magnoliopsida", "Caryophyllales", "Caryophyllaceae")

  id_addtax <- tax2$query == "Swainsona recta"
  tax2[id_addtax, tax_ranks] <- c("Plantae", "Streptophyta", "Magnoliopsida", "Fabales", "Fabaceae")

  id_addtax <- tax2$query == "Medicago sativa falcata"
  tax2[id_addtax, tax_ranks] <- c("Plantae", "Streptophyta", "Magnoliopsida", "Fabales", "Fabaceae")

  id_addtax <- tax2$query == "Leucochrysum albicans tricolor"
  tax2[id_addtax, tax_ranks] <- c("Plantae", "Streptophyta", "Magnoliopsida", "Asterales", "Asteraceae")

  id_addtax <- tax2$query == "Ardisia crenata bicolor"
  tax2[id_addtax, tax_ranks] <- c("Plantae", "Streptophyta", "Magnoliopsida", "Primulales", "Myrsinaceae")

  id_addtax <- tax2$query == "Pheidole annexus"
  tax2[id_addtax, tax_ranks] <- c("Metazoa", "Arthropoda", "Insecta", "Hymenoptera", "Formicidae")

  id_addtax <- tax2$query == "Persoonia glaucescens"
  tax2[id_addtax, tax_ranks] <- c("Plantae", "Streptophyta", "Magnoliopsida", "Proteales", "Proteaceae")

  id_addtax <- tax2$query == "Gnypetoscincus queenslandiae"
  tax2[id_addtax, "class"] <- "Reptilia"

  id_addtax <- tax2$query == "Plestiodon skiltonianus"
  tax2[id_addtax, "class"] <- "Reptilia"

  for (k in rev(seq_along(tax_ranks))) {
    if (anyNA(tax2[, tax_ranks[k]])) {
      print(tax_ranks[k])
      print(tax2[is.na(tax2[, tax_ranks[k]]), ])
      stop("Fill in missing taxonomy before continuing...")
    }
  }

  save(specieslist2, tax, tax2, file = ftaxonomy_new)
  if (!identical(ftaxonomy, ftaxonomy_new))
    unlink(ftaxonomy)
}


#--- Transfer taxonomic data
# Init
tax_ranks2 <- stringr::str_to_title(tax_ranks)
temp <- rep(NA, dim(dat0)[1])
for (k in seq_along(tax_ranks2)) {
  dat0[, tax_ranks2[k]] <- temp
}

for (k in seq_along(specieslist2)) {
  irows_dat <- which(tax2[k, "Species_resolved"] == dat0$Species_resolved)

  if (length(irows_dat) > 0) {
    dat0[irows_dat, tax_ranks2] <- as.data.frame(matrix(unlist(tax2[k, tax_ranks]),
      nrow = length(irows_dat), ncol = length(tax_ranks), byrow = TRUE),
      stringsAsFactors = FALSE)
  }
}


#--- Extract genus name
dat0[, "Genus"] <- sapply(strsplit(dat0[, "Species_resolved"], split = "\\b\\W+\\b"),
  function(x) x[[1]])


#------------------------------------
#--- SPECIES TRAITS

#--- Prepare spreadsheet for trait extraction from publications
if (FALSE) {
  temp <- c("Species_resolved", "Class", "Order", "Family", "Organism_group1",
  "Paper_Nr", "ID_article")
  temp <- unique(dat0[, temp])

  write.csv(temp[order(temp[, "Class"]), ], file = file.path(dir_out, "20170130_SpeciesData.csv"))
}


#--- eol TraitBank: http://www.eol.org/info/traitbank

# Available traits/moderators
# TraitBank currently features over 11 million records related to more than 330 attributes for 1.7 million taxa obtained from over 50 data sources.

# from 67 data sources: http://eol.org/collections/97700

# Suggestions Dec 15, 2016:
eol_traits <- list(
  body_size = c("abdomen length", "body length", "body length (VT)",
    "Body length, nose to tail", "body mass", "plant height", "weight"),
  generation_time = c("age at first birth", "age at first reproduction",
    "age at maturity", "age at maturity of resprouts", "dispersal age",
    "onset of fertility"),
  lifespan = c("life span", "total life span"), #"annual", "ephemeral", "biennial", "perennial"
  dispersal = c("dispersal vector", "propagule", "seed mass"),
#   - trophic level
  reproductive_syndrome = c("Mating System", "reproduction"), #"sexual system", "asexual reproduction", "cooperative breeding"
  food_specialization = c("diet breadth"),
  habitat_specialization = c("habitat breadth", "geographic range (size of area)"),
#   - climate specialization --> latitudinal/elevational range
#   - rarity --> red list status; local abundance; threats in addition to fragmentation
  rarity = c("conservation status")
)

do_eol <- TRUE
test_eol <- NULL

if (file.exists(feol)) {
  test_eol <- new.env(parent = emptyenv())
  load(feol, envir = test_eol)
  do_eol <- !all(specieslist2 %in% test_eol$specieslist2)
}

if (!do_eol) {
  ids <- match(specieslist2, test_eol$specieslist2, nomatch = 0)

  if (identical(specieslist2, test_eol$specieslist2[ids])) {
    eol_traits <- get("eol_traits", envir = test_eol)[ids]
    eol_traits2 <- get("eol_traits2", envir = test_eol)[ids]
  } else {
    do_eol <- TRUE
  }
}


if (do_eol && taxize::eol_ping() &&
  (requireNamespace("curl") && curl::has_internet()) || !requireNamespace("curl")) {

  eol_tb_predicates <- unlist(eol_traits)

  eol_traits <- NULL
  eol_traits2_cols <- c("@id", "dwc:scientificname", "predicate", "value", "units")
  eol_traits2 <- as.data.frame(matrix(NA, nrow = 0, ncol = length(eol_traits2_cols),
    dimnames = list(NULL, eol_traits2_cols)))

  for (k in seq_along(specieslist2)) {
    temp <- taxize::eol_search(terms = specieslist2[k], key = eol_API_key, cache_ttl = 10)
    pageid <- if (NROW(temp) > 0 && "pageid" %in% names(temp)) temp[1, "pageid"] else NULL

    if (all(is.finite(pageid)) && !is.null(pageid)) {
        res <- traits::traitbank(pageid, cache_ttl = 10)[["graph"]]

        if (NROW(res) > 0) {
          eol_traits[[specieslist2[k]]] <- res

          add_data <- as.data.frame(matrix(NA, nrow = nrow(res),
            ncol = length(eol_traits2_cols), dimnames = list(NULL, eol_traits2_cols)))
          icol <- match(eol_traits2_cols, colnames(res), nomatch = 0)
          add_data[, eol_traits2_cols %in% colnames(res)] <- res[, icol]
          eol_traits2 <- rbind(eol_traits2, add_data)

        } else {
          print(paste(k, specieslist2[k], ": traitbank found no traits:",
            pageid))
        }

    } else {
      print(paste(k, specieslist2[k], ": has invalid eol-IDs:", paste(pageid,
        collapse = ", ")))
    }
  }

}


if (do_eol) {
  save(specieslist2, eol_traits, eol_traits2, file = feol_new)
  if (!identical(feol, feol_new))
    unlink(feol)
}

print(length(eol_traits)) # --> 99 species
print(sort(unique(eol_traits2[, "predicate"]))) # --> 195 traits
print(nrow(eol_traits2)) # --> 6659 trait records

### Selected traits
# Animals:
#   - body size
#   - generation time:
#   - dispersal: flying, non-flying
#   - sex: sex, no sex, mixed
#   - trophic group: herbivore, predator, parasite, parasitoid, omnivore, detrivore
#   - breeding system: e.g., monogamy, polygyny, polyandry, lekking

# Plants:
#   - Raunkiaer lifeform
#   - Age
#   - sex: sex, no sex, mixed
#   - Self compatibility: self-compatible, obligate outcrossing, self-compatible but mainly outcrossing
#   - Pollination syndrome: wind, wind and pollinators, mostly one pollinator, several pollinators, many pollinators
#   - Pollination vector: wind, insect, bird, bat
#   - Seed dispersal vector: wind, ants, bats, mammals, birds & mammals, gravity, ejection, water

used_predicates <- c("age at first birth", "age at first reproduction", "age at maturity",
  "body length (VT)", "body mass", "developmental mode", "dispersal vector", "feeds on",
  "inter-birth interval", "life span", "Mating System", "onset of fertility", "preys upon",
  "primary diet", "primary growth form", "propagule", "reproduction", "total life span",
  "weight")


eol_traits3 <- eol_traits2[eol_traits2$predicate %in% used_predicates, ]
eol_traits4 <- reshape2::dcast(reshape2::melt(eol_traits3),
  `dwc:scientificname` + `@id` ~ predicate + units)
eol_traits4 <- eol_traits4[apply(eol_traits4, 1, function(x) any(!is.na(x))), ]

if (FALSE)
  write.csv(eol_traits4, file = file.path(dir_in, "20170220_SpeciesTraits_EOLtraitbank.csv"))



#--- IUCN red list: www.iucnredlist.org

# Available traits/moderators
#   - geographic range
#   - habitat and ecology
#   - habitat suitability
#   - systems (e.g., terrestrial)


#--- Plants Trait Database: www.try-db.org

# Available traits
traits_tryDB <- read.table(file.path(dir_in, "tryDB_tde201612712946.txt"), header = TRUE, skip = 3, sep = "\t")

# traits that cover most plant species
traits_tryDB[order(traits_tryDB$AccSpecNum, decreasing = TRUE), ][1:20, ]
    #    TraitID                                                                          Trait ObsNum ObsGRNum PubNum AccSpecNum  X
    #685      38                                                                Plant woodiness 108763    32530  75967      54466 NA
    #605      42                                                              Plant growth form 294097   192688 126025      38164 NA
    #658      31                                                       Plant tolerance to frost  34060       NA  33805      29425 NA
    #608      18                                                                   Plant height 226560   130159 132824      26837 NA
    #763      26                                                                  Seed dry mass 183083     4803  19922      26095 NA
    #430      37                                                            Leaf phenology type  93015    47079  43678      20341 NA
    #580      22                                                         Photosynthetic pathway  46471     7675  28097      17919 NA
    #766      95                                 Seed germination rate (germination efficiency)  63202       NA     NA      15585 NA
    #236       1                                                                      Leaf area 129720   111720  55281      15238 NA
    #281      17                                                              Leaf compoundness  36138    16287  21850      13387 NA
    #240      11                          Leaf area per leaf dry mass (specific leaf area, SLA) 137112   107515  50922      11991 NA
    #614      59                                                     Plant lifespan (longevity)  43446    10324  36267      11609 NA
    #628       8                                               Plant nitrogen fixation capacity  36942    19551   9232      11569 NA
    #506      43                                                                      Leaf type  46256    31658  15530      11055 NA
    #790      98                                                         Seed storage behaviour  12492       NA     NA      10291 NA
    #613     343                                          Plant life form (Raunkiaer life form)  82710    60915  16635      10169 NA
    #399      14                                    Leaf nitrogen (N) content per leaf dry mass  64835    48534  14108      10149 NA
    #938       4 Stem dry mass per stem fresh volume (stem specific density, SSD, wood density)  35076    16190  22843      10147 NA
    #146      28                                                             Dispersal syndrome 354368      817 299601      10045 NA
    #478     154                                                                     Leaf shape  35406     6857  33028       9085 NA



#------------------------------------
#--- Species traits extracted from articles

#--- Plants
temp <- names(gdata::read.xls(fname_traits_plants, as.is = TRUE, nrows = 1))
dat_traits_plants <- gdata::read.xls(fname_traits_plants, as.is = TRUE, skip = 1, header = FALSE)
names(dat_traits_plants) <- temp
dat_traits_plants[dat_traits_plants == ""] <- NA

itemp <- match(dat0$ID_article, dat_traits_plants$Article_ID, nomatch = 0)
iuse <- itemp > 0
if (grepl("SpeciesTraits_Plants16.3.17", basename(fname_traits_plants))) {
  ctemp <- c("Lifeform", "Age..years.", "Sex", "Self_compatibility", "Mating_system",
    "Pollination_syndrome", "Pollination_vector", "Seed_dispersal_vector", "Remarks")
  ctemp_new <- paste0("Plants_", c("Lifeform", "Age_yr", "Sex", "Self_compatibility",
    "Mating_system", "Pollination_syndrome", "Pollination_vector",
    "Seed_dispersal_vector", "TraitsRemarks"))
} else if (grepl("20170411", basename(fname_traits_plants))) {
  ctemp <- c("Lifeform", "Age_years", "Age_class", "Sex", "Self_compatibility", "Mating_system",
    "Pollination_syndrome", "Pollination_vector", "Seed_dispersal_vector")
  ctemp_new <- paste0("Plants_", ctemp)
} else {
  stop("unknown plant trait file")
}
for (k in seq_along(ctemp))
  dat0[iuse, ctemp_new[k]] <- as.character(dat_traits_plants[itemp, ctemp[k]])

#--- Sort the plant trait factors
#   * Life span
unique(dat0$Plants_Age_years)

temp <- rep(NA, nrow(dat0))
temp[dat0$Plants_Age_years %in% c("1")] <- 1 # 1-2 yrs
temp[dat0$Plants_Age_years %in% c("<10", "<5", "3-10", "5-8")] <- 2 # 3-10 yrs
temp[dat0$Plants_Age_years %in% c(">10", "13-18")] <- 3 # 11-15 yrs
temp[dat0$Plants_Age_years %in% c("<50", ">15", ">25", "20-50", "30")] <- 4 # 16-50 yrs
temp[dat0$Plants_Age_years %in% c("<100", ">50", "50-100")] <- 5 # 51-100 yrs
temp[dat0$Plants_Age_years %in% c(">100")] <- 6 # 101-200 yrs
temp[dat0$Plants_Age_years %in% c(">200")] <- 7 # >200 yrs
dat0[, "Plants_Age_cat"] <- factor(temp, levels = seq_along(age_categories),
  labels = age_categories, ordered = TRUE)

table(dat0$Plants_Age_cat, dat0$Plants_Age_years)
temp <- sapply(list(dat0$Plants_Age_cat, dat0$Plants_Age_years), function(x) sum(!is.na(x)))
stopifnot(temp[1] >= temp[2])

#   * Plant life form (Raunkiaer 1934)
dat0[, "Plants_Lifeform"] <- factor(dat0[, "Plants_Lifeform"], ordered = FALSE)

#   * Mode of reproduction
#     1 = exclusive sexual reproduction
#     2 = both sexual and vegetative reproduction
dat0[, "Plants_Sex"] <- factor(dat0[, "Plants_Sex"], ordered = FALSE)

#   * Compatibility system
#     CS = self-compatibility
#     SI = self-incompatibility
dat0[, "Plants_Self_compatibility"] <- factor(dat0[, "Plants_Self_compatibility"], ordered = FALSE)

#   * Mating system
#     1 = mainly outcrossing
#     2 = mainly selfing
#     3 = mixed (both outcrossing and selfing)
dat0[, "Plants_Mating_system"] <- factor(dat0[, "Plants_Mating_system"], ordered = FALSE)

#   * Pollination syndrome
#     1 = wind
#     2 = various insects
#     3 = birds
dat0[, "Plants_Pollination_syndrome"] <- factor(dat0[, "Plants_Pollination_syndrome"], ordered = FALSE)

#   * Dominant pollen vector
#     1 = wind
#     2 = insects
#     3 = birds and/or bats
dat0[, "Plants_Pollination_vector"] <- factor(dat0[, "Plants_Pollination_vector"], ordered = FALSE)

#   * Seed dispersal
#     1 = gravity
#     2 = short distance dispersed seeds (ejection, ants)
#     3 = long distance dispersed seeds (wind, birds, bats)
dat0[, "Plants_Seed_dispersal_vector"] <- factor(dat0[, "Plants_Seed_dispersal_vector"], ordered = FALSE)


#--- Animals
temp <- names(gdata::read.xls(fname_traits_animals, sheet = "Data", as.is = TRUE, nrows = 1))
dat_traits_animals <- gdata::read.xls(fname_traits_animals, sheet = "Data", as.is = TRUE, skip = 2, header = FALSE)
names(dat_traits_animals) <- temp
dat_traits_animals[dat_traits_animals == ""] <- NA

itemp <- match(dat0$ID_article, dat_traits_animals$Article_ID, nomatch = 0)
iuse <- itemp > 0
ctemp <- c("Max_adult_lenght", "Max_adult_mass", "Generation_lenght", "Max_age",
  "Dispersal_fine", "Trophic_group")
ctemp_new <- paste0("Animals_", c("Length_cat", "Mass_cat", "Generation_class", "Age_class",
  "Dispersal_ability", "Trophic_group"))
for (k in seq_along(ctemp))
  dat0[iuse, ctemp_new[k]] <- as.character(dat_traits_animals[itemp, ctemp[k]])

#--- Sort the animal trait factors
#   * Body length based on maximum reported including tails
#     tiny: ≤ 2 cm; small: 2.01–10 cm; medium: 10.01–50 cm; large: > 50 cm
cats <- c("tiny", "small", "medium", "large")
temp <- sapply(dat0[, "Animals_Length_cat"], function(x)
  switch(x, tiny = 1, small = 2, medium = 3, large = 4, NA))
dat0[, "Animals_Length_cat"] <- factor(temp, levels = seq_along(cats), labels = cats, ordered = TRUE)

#   * Body mass based on maximum reported for adult
#     small: < 15 g; medium: 15–99.99; large: 100–999.99; huge: 1000
cats <- c("small", "medium", "large", "huge")
temp <- sapply(dat0[, "Animals_Mass_cat"], function(x)
  switch(x, small = 1, medium = 2, large = 3, huge = 4, NA))
dat0[, "Animals_Mass_cat"] <- factor(temp, levels = seq_along(cats), labels = cats, ordered = TRUE)

#   * Generation length typically based on age at sexual maturity thus giving an
#     indication of maximum number of generations that could have lived in the fragments
#     since fragmentation
#     short: < 2 years; medium: 2–4.99 years; long: ≥ 5 years
cats <- c("small", "medium", "long")
temp <- sapply(dat0[, "Animals_Generation_class"], function(x)
  switch(x, small = 1, medium = 2, long = 3, NA))
dat0[, "Animals_Generation_class"] <- factor(temp, levels = seq_along(cats), labels = cats, ordered = TRUE)

#   * Longevity based on maximum reported lifespan including captivity
#     30 yrs is max among this set of animal species
#     short: ≤ 3 years; medium: 3.11–9.99; long: 10–14.99; very_long: ≥ 15
cats <- c("short", "medium", "long", "very_long")
temp <- sapply(dat0[, "Animals_Age_class"], function(x)
  switch(x, short = 1, medium = 2, long = 3, very_long = 4, NA))
dat0[, "Animals_Age_class"] <- factor(temp, levels = seq_along(cats), labels = cats, ordered = TRUE)

#   * Dispersal ability
#      walk_low: < 100 m dispersal; includes species simply listed as not moving much or
#                not having dispersal stage and small home ranges
#      walk_medium: 100.01–1000 m dispersal
#      walk_long: > 1000 m dispersal; includes some very wide ranging species
#      flying: Capable of active flight during at least some stages or for a part of the
#               population, thus including ants because of their flying queens or males,
#               or arthropods where only a small part of the population is macropter,
#               however, does not include gliders, balloning would have been accepted but
#               did not occur in sample
cats <- c("walk_low", "walk_medium", "walk_long", "flying")
temp <- sapply(dat0[, "Animals_Dispersal_ability"], function(x)
  switch(x, walk_low = 1, walk_medium = 2, walk_long = 3, flying = 4, NA))
dat0[, "Animals_Dispersal_ability"] <- factor(temp, levels = seq_along(cats),
  labels = cats, ordered = TRUE)

#   * Trophic group based on main diet of adult
#      herbivore: mostly eats plant products like leaves, fruit, nectar, seeds;
#                 includes some species that also take some invertebrates but clearly
#                 mostly eat plant matter
#      omnivore: eats plant matter and animals
#      mycophage: eats mostly fungi, e.g. truffle specialist
#      predator: eats mostly animals that it overpowers; includes some species where
#                juveniles have a different diet, or where adults have seasonaly varying
#                diets (e.g. birds that eat seeds in winter but prey on insects in most
#                of the year
cats <- c("herbivore", "omnivore", "mycophage", "predator")
dat0[, "Animals_Trophic_group"] <- factor(dat0[, "Animals_Trophic_group"],
  labels = cats, ordered = FALSE)


#--- Animal age in years and generation lenght
temp <- names(gdata::read.xls(fname_age_animals, sheet = "Sheet1", as.is = TRUE, nrows = 1))
dat_age_animals <- gdata::read.xls(fname_age_animals, sheet = "Sheet1", as.is = TRUE, skip = 2, header = FALSE)
names(dat_age_animals) <- temp
dat_age_animals[dat_age_animals == ""] <- NA

stopifnot(all(dat_age_animals$Species %in% dat0$Species_resolved))
itemp <- match(dat0$Species_resolved, dat_age_animals$Species, nomatch = 0)
iuse <- itemp > 0
ctemp <- c("Age", "Generation_lenght")
ctemp_new <- paste0("Animals_", c("Age_years", "Generation_years"))
for (k in seq_along(ctemp))
  dat0[iuse, ctemp_new[k]] <- as.character(dat_age_animals[itemp, ctemp[k]])

# Check that all information between 'class' and 'years' variables is consistent
table(dat0$Animals_Age_class, dat0$Animals_Age_years)
temp <- sapply(list(dat0$Animals_Age_class, dat0$Animals_Age_years),
  function(x) sum(!is.na(x)))
stopifnot(temp[1] >= temp[2])

temp <- rep(NA, nrow(dat0))
temp_num <- as.numeric(dat0$Animals_Age_years)
temp[(!is.na(temp_num) & temp_num <= 2) |
    (is.na(temp_num) & dat0$Animals_Age_class %in% "short")] <- 1 # 1-2 yrs
temp[(!is.na(temp_num) & temp_num > 2 & temp_num <= 10) |
    (is.na(temp_num) & dat0$Animals_Age_class %in% "medium")] <- 2 # 3-10 yrs
temp[(!is.na(temp_num) & temp_num > 10 & temp_num <= 15) |
    (is.na(temp_num) & dat0$Animals_Age_class %in% "long")] <- 3 # 11-15 yrs
temp[(!is.na(temp_num) & temp_num > 15 & temp_num <= 50) |
    (is.na(temp_num) & dat0$Animals_Age_class %in% "very long")] <- 4 # 16-50 yrs
temp[!is.na(temp_num) & temp_num > 50 & temp_num <= 100] <- 5 # 51-100 yrs
temp[!is.na(temp_num) & temp_num > 100 & temp_num <= 200] <- 6 # 101-200 yrs
temp[!is.na(temp_num) & temp_num > 200] <- 7 # >200 yrs
dat0[, "Animals_Age_cat"] <- factor(temp, levels = seq_along(age_categories),
  labels = age_categories, ordered = TRUE)
table(dat0$Animals_Age_cat, dat0$Animals_Age_years)


table(dat0$Animals_Generation_class, dat0$Animals_Generation_years)
temp <- sapply(list(dat0$Animals_Generation_class, dat0$Animals_Generation_years),
  function(x) sum(!is.na(x)))
#     short: < 2 years; medium: 2–4.99 years; long: ≥ 5 years

if (temp[1] < temp[2]) {
  temp <- rep(NA, nrow(dat0))
  temp_num <- as.numeric(dat0$Animals_Generation_years)
  temp[dat0$Animals_Generation_years %in% c("<<1")] <- 1 # 1-2 yrs
  temp[(!is.na(temp_num) & temp_num <= 2) |
      (is.na(temp_num) & dat0$Animals_Generation_class %in% "short")] <- 1 # 1-2 yrs
  temp[(!is.na(temp_num) & temp_num > 2 & temp_num <= 10) |
      (is.na(temp_num) & dat0$Animals_Generation_class %in% c("medium", "long"))] <- 2 # 3-10 yrs
  temp[!is.na(temp_num) & temp_num > 10 & temp_num <= 15] <- 3 # 11-15 yrs
  temp[!is.na(temp_num) & temp_num > 15 & temp_num <= 50] <- 4 # 16-50 yrs
  temp[!is.na(temp_num) & temp_num > 50 & temp_num <= 100] <- 5 # 51-100 yrs
  temp[!is.na(temp_num) & temp_num > 100 & temp_num <= 200] <- 6 # 101-200 yrs
  temp[!is.na(temp_num) & temp_num > 200] <- 7 # >200 yrs
  dat0[, "Animals_Generation_cat"] <- factor(temp, levels = seq_along(age_categories),
    labels = age_categories, ordered = TRUE)

  table(dat0$Animals_Generation_cat, dat0$Animals_Generation_years)
  table(dat0$Animals_Generation_cat, dat0$Animals_Generation_class)
} else {
  stop()
}

#--- Species age
table(dat0$Plants_Age_cat)
table(dat0$Plants_Age_years)
table(dat0$Animals_Age_class)
table(dat0$Animals_Age_cat)
table(dat0$Animals_Age_years)
table(dat0$Animals_Generation_class)
table(dat0$Animals_Generation_cat)
table(dat0$Animals_Generation_years)

temp <- rep(NA, nrow(dat0))
for (k in seq_along(age_categories)) {
  ids <- dat0$Plants_Age_cat == age_categories[k] |
    dat0$Animals_Age_cat == age_categories[k]
  temp[ids] <- k
}
dat0[, "SpeciesAgeCat_ord_yrs"] <- factor(temp, levels = seq_along(age_categories),
  labels = age_categories, ordered = TRUE)

utemp <- unique(dat0[, c("ID_effect", "SpeciesAgeCat_ord_yrs")])
print(table(utemp[, "SpeciesAgeCat_ord_yrs"]))
print(table(dat0[, "Plants_Age_cat"], dat0[, "SpeciesAgeCat_ord_yrs"]))
print(table(dat0[, "Animals_Age_cat"], dat0[, "SpeciesAgeCat_ord_yrs"]))


# Aggregated age categories
temp <- rep(NA, nrow(dat0))
temp[dat0[, "SpeciesAgeCat_ord_yrs"] %in% c("(0, 2] yrs", "(2-10] yrs")] <- 1
temp[dat0[, "SpeciesAgeCat_ord_yrs"] %in% c("(10-15] yrs", "(15-50] yrs")] <- 2
temp[dat0[, "SpeciesAgeCat_ord_yrs"] %in% c("(50-100] yrs", "(100-200] yrs", ">200 yrs")] <- 3
dat0[, "SpeciesAgeCat2_ord_yrs"] <- factor(temp, levels = seq_along(age_categories2),
  labels = age_categories2, ordered = TRUE)


#------------------------------------
#--- Species lifespan relative to fragment age (~ minimum generations since fragmentation)
table(dat0[, "FragmentAgeCat_ord_yrs"], dat0[, "SpeciesAgeCat_ord_yrs"], useNA = "always",
  deparse.level = 2)
table(dat0[, "EffectAgeCat_ord_yrs"], dat0[, "SpeciesAgeCat_ord_yrs"], useNA = "always",
  deparse.level = 2)
#  ==> categories: <= 1 generation since fragmentation; 2-5; > 5
relage_categories <- c("(0, 1] Tf/Nl", "(1, 5] Tf/Nl", ">5 Tf/Nl") # Time(fragmentation) / N(lifespan)

temp <- rep(NA, nrow(dat0))

#   * Case 1: numeric fragment age && numeric animal age
num_relage1 <-  as.numeric(dat0$Fragment_age_yr) / as.numeric(dat0$Animals_Age_years)
num_relage1_cat <- cut(num_relage1, breaks = c(0, 1, 5, Inf), labels = relage_categories,
  ordered_result = TRUE)
ids <- !is.na(num_relage1_cat)
temp[ids] <- num_relage1_cat[ids]

#   * Case 2: ordinal fragment age && numeric animal age
num_hasnona <- !is.na(dat0$Animals_Age_years)
ids <- c(
  which(dat0$Fragment_age_category %in% ">22" & num_hasnona & dat0$Animals_Age_years <= 22),
  which(dat0$Fragment_age_category %in% "20-90" & num_hasnona & dat0$Animals_Age_years <= 20),
  which(dat0$Fragment_age_category %in% ">30" & num_hasnona & dat0$Animals_Age_years <= 30),
  which(dat0$Fragment_age_category %in% ">40" & num_hasnona & dat0$Animals_Age_years <= 40),
  which(dat0$Fragment_age_category %in% "<49" & num_hasnona & dat0$Animals_Age_years <= 49),
  which(dat0$Fragment_age_category %in% ">49" & num_hasnona & dat0$Animals_Age_years <= 49),
  which(dat0$Fragment_age_category %in% "<50" & num_hasnona & dat0$Animals_Age_years <= 50),
  which(dat0$Fragment_age_category %in% ">100" & num_hasnona & dat0$Animals_Age_years <= 100),
  which(dat0$Fragment_age_category %in% ">300" & num_hasnona & dat0$Animals_Age_years <= 300),
  which(dat0$Fragment_age_category %in% ">500" & num_hasnona & dat0$Animals_Age_years <= 500),
  which(dat0$Fragment_age_category %in% ">1000" & num_hasnona & dat0$Animals_Age_years <= 1000)
)
stopifnot(is.na(temp[ids]))
temp[ids] <- 1 # <= 1 generation

ids <- c(
  which(dat0$Fragment_age_category %in% ">22" &
    num_hasnona & dat0$Animals_Age_years > 22 & dat0$Animals_Age_years <= 5 * 22),
  which(dat0$Fragment_age_category %in% "20-90" &
    num_hasnona & dat0$Animals_Age_years > 20 & dat0$Animals_Age_years <= 5 * 90),
  which(dat0$Fragment_age_category %in% ">30" &
    num_hasnona & dat0$Animals_Age_years > 30 & dat0$Animals_Age_years <= 5 * 30),
  which(dat0$Fragment_age_category %in% ">40" &
    num_hasnona & dat0$Animals_Age_years > 40 & dat0$Animals_Age_years <= 5 * 40),
  which(dat0$Fragment_age_category %in% "<49" &
    num_hasnona & dat0$Animals_Age_years > 49 & dat0$Animals_Age_years <= 5 * 49),
  which(dat0$Fragment_age_category %in% ">49" &
    num_hasnona & dat0$Animals_Age_years > 49 & dat0$Animals_Age_years <= 5 * 49),
  which(dat0$Fragment_age_category %in% "<50" &
    num_hasnona & dat0$Animals_Age_years > 50 & dat0$Animals_Age_years <= 5 * 50),
  which(dat0$Fragment_age_category %in% ">100" &
    num_hasnona & dat0$Animals_Age_years > 100 & dat0$Animals_Age_years <= 5 * 100),
  which(dat0$Fragment_age_category %in% ">300" &
    num_hasnona & dat0$Animals_Age_years > 300 & dat0$Animals_Age_years <= 5 * 300),
  which(dat0$Fragment_age_category %in% ">500" &
    num_hasnona & dat0$Animals_Age_years > 500 & dat0$Animals_Age_years <= 5 * 500),
  which(dat0$Fragment_age_category %in% ">1000" &
    num_hasnona & dat0$Animals_Age_years > 1000 & dat0$Animals_Age_years <= 5 * 1000)
)
stopifnot(is.na(temp[ids]))
temp[ids] <- 2 # >1 to <= 5 generations

ids <- c(
  which(dat0$Fragment_age_category %in% ">22" & num_hasnona & dat0$Animals_Age_years > 5 * 22),
  which(dat0$Fragment_age_category %in% "20-90" & num_hasnona & dat0$Animals_Age_years > 5 * 90),
  which(dat0$Fragment_age_category %in% ">30" & num_hasnona & dat0$Animals_Age_years > 5 * 30),
  which(dat0$Fragment_age_category %in% ">40" & num_hasnona & dat0$Animals_Age_years > 5 * 40),
  which(dat0$Fragment_age_category %in% "<49" & num_hasnona & dat0$Animals_Age_years > 5 * 49),
  which(dat0$Fragment_age_category %in% ">49" & num_hasnona & dat0$Animals_Age_years > 5 * 49),
  which(dat0$Fragment_age_category %in% "<50" & num_hasnona & dat0$Animals_Age_years > 5 * 50),
  which(dat0$Fragment_age_category %in% ">100" & num_hasnona & dat0$Animals_Age_years > 5 * 100),
  which(dat0$Fragment_age_category %in% ">300" & num_hasnona & dat0$Animals_Age_years > 5 * 300),
  which(dat0$Fragment_age_category %in% ">500" & num_hasnona & dat0$Animals_Age_years > 5 * 500),
  which(dat0$Fragment_age_category %in% ">1000" & num_hasnona & dat0$Animals_Age_years > 5 * 1000)
)
stopifnot(is.na(temp[ids]))
temp[ids] <- 3 # > 5 generations

#   * Case 3: numeric fragment age && ordinal plant age
num_hasnona <- !is.na(dat0$Fragment_age_yr)
ids <- c(
  which(dat0$Plants_Age_years %in% "1" & num_hasnona & dat0$Fragment_age_yr <= 1),
  which(dat0$Plants_Age_years %in% "<5" & num_hasnona & dat0$Fragment_age_yr <= 5),
  which(dat0$Plants_Age_years %in% "3-10" & num_hasnona & dat0$Fragment_age_yr <= 3),
  which(dat0$Plants_Age_years %in% "5-8" & num_hasnona & dat0$Fragment_age_yr <= 5),
  which(dat0$Plants_Age_years %in% "<10" & num_hasnona & dat0$Fragment_age_yr <= 10),
  which(dat0$Plants_Age_years %in% ">10" & num_hasnona & dat0$Fragment_age_yr <= 10),
  which(dat0$Plants_Age_years %in% "13-18" & num_hasnona & dat0$Fragment_age_yr <= 13),
  which(dat0$Plants_Age_years %in% ">15" & num_hasnona & dat0$Fragment_age_yr <= 15),
  which(dat0$Plants_Age_years %in% "20-50" & num_hasnona & dat0$Fragment_age_yr <= 20),
  which(dat0$Plants_Age_years %in% ">25" & num_hasnona & dat0$Fragment_age_yr <= 25),
  which(dat0$Plants_Age_years %in% "30" & num_hasnona & dat0$Fragment_age_yr <= 30),
  which(dat0$Plants_Age_years %in% "<50" & num_hasnona & dat0$Fragment_age_yr <= 50),
  which(dat0$Plants_Age_years %in% ">50" & num_hasnona & dat0$Fragment_age_yr <= 50),
  which(dat0$Plants_Age_years %in% "50-100" & num_hasnona & dat0$Fragment_age_yr <= 50),
  which(dat0$Plants_Age_years %in% "<100" & num_hasnona & dat0$Fragment_age_yr <= 100),
  which(dat0$Plants_Age_years %in% ">100" & num_hasnona & dat0$Fragment_age_yr <= 100),
  which(dat0$Plants_Age_years %in% ">200" & num_hasnona & dat0$Fragment_age_yr <= 200)
)
stopifnot(is.na(temp[ids]))
temp[ids] <- 1 # <= 1 generation

ids <- c(
  which(dat0$Plants_Age_years %in% "1" &
    num_hasnona & dat0$Fragment_age_yr > 1 & dat0$Fragment_age_yr <= 5 * 1),
  which(dat0$Plants_Age_years %in% "<5" &
    num_hasnona & dat0$Fragment_age_yr > 5 & dat0$Fragment_age_yr <= 5 * 5),
  which(dat0$Plants_Age_years %in% "3-10" &
    num_hasnona & dat0$Fragment_age_yr > 3 & dat0$Fragment_age_yr <= 5 * 3),
  which(dat0$Plants_Age_years %in% "5-8" &
    num_hasnona & dat0$Fragment_age_yr > 5 & dat0$Fragment_age_yr <= 5 * 5),
  which(dat0$Plants_Age_years %in% "<10" &
    num_hasnona & dat0$Fragment_age_yr > 10 & dat0$Fragment_age_yr <= 5 * 10),
  which(dat0$Plants_Age_years %in% ">10" &
    num_hasnona & dat0$Fragment_age_yr > 10 & dat0$Fragment_age_yr <= 5 * 10),
  which(dat0$Plants_Age_years %in% "13-18" &
    num_hasnona & dat0$Fragment_age_yr > 13 & dat0$Fragment_age_yr <= 5 * 13),
  which(dat0$Plants_Age_years %in% ">15" &
    num_hasnona & dat0$Fragment_age_yr > 15 & dat0$Fragment_age_yr <= 5 * 15),
  which(dat0$Plants_Age_years %in% "20-50" &
    num_hasnona & dat0$Fragment_age_yr > 20 & dat0$Fragment_age_yr <= 5 * 20),
  which(dat0$Plants_Age_years %in% ">25" &
    num_hasnona & dat0$Fragment_age_yr > 25 & dat0$Fragment_age_yr <= 5 * 25),
  which(dat0$Plants_Age_years %in% "30" &
    num_hasnona & dat0$Fragment_age_yr > 30 & dat0$Fragment_age_yr <= 5 * 30),
  which(dat0$Plants_Age_years %in% "<50" &
    num_hasnona & dat0$Fragment_age_yr > 50 & dat0$Fragment_age_yr <= 5 * 50),
  which(dat0$Plants_Age_years %in% ">50" &
    num_hasnona & dat0$Fragment_age_yr > 50 & dat0$Fragment_age_yr <= 5 * 50),
  which(dat0$Plants_Age_years %in% "50-100" &
    num_hasnona & dat0$Fragment_age_yr > 50 & dat0$Fragment_age_yr <= 5 * 50),
  which(dat0$Plants_Age_years %in% "<100" &
    num_hasnona & dat0$Fragment_age_yr > 100 & dat0$Fragment_age_yr <= 5 * 100),
  which(dat0$Plants_Age_years %in% ">100" &
    num_hasnona & dat0$Fragment_age_yr > 100 & dat0$Fragment_age_yr <= 5 * 100),
  which(dat0$Plants_Age_years %in% ">200" &
    num_hasnona & dat0$Fragment_age_yr > 200 & dat0$Fragment_age_yr <= 5 * 200)
)
stopifnot(is.na(temp[ids]))
temp[ids] <- 2 # >1 to <= 5 generations

ids <- c(
  which(dat0$Plants_Age_years %in% "1" & num_hasnona & dat0$Fragment_age_yr > 5 * 1),
  which(dat0$Plants_Age_years %in% "<5" & num_hasnona & dat0$Fragment_age_yr > 5 * 5),
  which(dat0$Plants_Age_years %in% "3-10" & num_hasnona & dat0$Fragment_age_yr > 5 * 3),
  which(dat0$Plants_Age_years %in% "5-8" & num_hasnona & dat0$Fragment_age_yr > 5 * 5),
  which(dat0$Plants_Age_years %in% "<10" & num_hasnona & dat0$Fragment_age_yr > 5 * 10),
  which(dat0$Plants_Age_years %in% ">10" & num_hasnona & dat0$Fragment_age_yr > 5 * 10),
  which(dat0$Plants_Age_years %in% "13-18" & num_hasnona & dat0$Fragment_age_yr > 5 * 13),
  which(dat0$Plants_Age_years %in% ">15" & num_hasnona & dat0$Fragment_age_yr > 5 * 15),
  which(dat0$Plants_Age_years %in% "20-50" & num_hasnona & dat0$Fragment_age_yr > 5 * 20),
  which(dat0$Plants_Age_years %in% ">25" & num_hasnona & dat0$Fragment_age_yr > 5 * 25),
  which(dat0$Plants_Age_years %in% "30" & num_hasnona & dat0$Fragment_age_yr > 5 * 30),
  which(dat0$Plants_Age_years %in% "<50" & num_hasnona & dat0$Fragment_age_yr > 5 * 50),
  which(dat0$Plants_Age_years %in% ">50" & num_hasnona & dat0$Fragment_age_yr > 5 * 50),
  which(dat0$Plants_Age_years %in% "50-100" & num_hasnona & dat0$Fragment_age_yr > 5 * 50),
  which(dat0$Plants_Age_years %in% "<100" & num_hasnona & dat0$Fragment_age_yr > 5 * 100),
  which(dat0$Plants_Age_years %in% ">100" & num_hasnona & dat0$Fragment_age_yr > 5 * 100),
  which(dat0$Plants_Age_years %in% ">200" & num_hasnona & dat0$Fragment_age_yr > 5 * 200)
)
stopifnot(is.na(temp[ids]))
temp[ids] <- 3 # > 5 generations

#   * Case 4: ordinal fragment age && ordinal plant age
page_set1 <- c("1", "<10", "<5", "3-10", "5-8", ">10", "13-18", ">15")
page_set2 <- c(page_set1, ">25", "20-50", "30")
page_set3 <- c(page_set2, "<50")
page_set4 <- c(page_set3, ">50", "50-100", "<100")
page_set5 <- c(page_set4, ">100")
page_set6 <- c(page_set5, ">200")

ids <- c(
  which(dat0$Fragment_age_category %in% ">22" & dat0$Plants_Age_years %in% page_set1),
  which(dat0$Fragment_age_category %in% "20-90" & dat0$Plants_Age_years %in% page_set1),
  which(dat0$Fragment_age_category %in% ">30" & dat0$Plants_Age_years %in% page_set2),
  which(dat0$Fragment_age_category %in% ">40" & dat0$Plants_Age_years %in% page_set2),
  which(dat0$Fragment_age_category %in% "<49" & dat0$Plants_Age_years %in% page_set3),
  which(dat0$Fragment_age_category %in% ">49" & dat0$Plants_Age_years %in% page_set3),
  which(dat0$Fragment_age_category %in% "<50" & dat0$Plants_Age_years %in% page_set3),
  which(dat0$Fragment_age_category %in% ">100" & dat0$Plants_Age_years %in% page_set4),
  which(dat0$Fragment_age_category %in% ">300" & dat0$Plants_Age_years %in% page_set6),
  which(dat0$Fragment_age_category %in% ">500" & dat0$Plants_Age_years %in% page_set6),
  which(dat0$Fragment_age_category %in% ">1000" & dat0$Plants_Age_years %in% page_set6)
)
stopifnot(is.na(temp[ids]))
temp[ids] <- 1 # <= 1 generation

ids <- c(
  which(dat0$Fragment_age_category %in% ">22" & !(dat0$Plants_Age_years %in% page_set1) & dat0$Plants_Age_years %in% page_set4),
  which(dat0$Fragment_age_category %in% "20-90" & !(dat0$Plants_Age_years %in% page_set1) & dat0$Plants_Age_years %in% page_set6),
  which(dat0$Fragment_age_category %in% ">30" & !(dat0$Plants_Age_years %in% page_set2) & dat0$Plants_Age_years %in% page_set5),
  which(dat0$Fragment_age_category %in% ">40" & !(dat0$Plants_Age_years %in% page_set2) & dat0$Plants_Age_years %in% page_set6),
  which(dat0$Fragment_age_category %in% "<49" & !(dat0$Plants_Age_years %in% page_set3) & dat0$Plants_Age_years %in% page_set6),
  which(dat0$Fragment_age_category %in% ">49" & !(dat0$Plants_Age_years %in% page_set3) & dat0$Plants_Age_years %in% page_set6),
  which(dat0$Fragment_age_category %in% "<50" & !(dat0$Plants_Age_years %in% page_set3) & dat0$Plants_Age_years %in% page_set6),
  which(dat0$Fragment_age_category %in% ">100" & !(dat0$Plants_Age_years %in% page_set4) & dat0$Plants_Age_years %in% page_set6),
  which(dat0$Fragment_age_category %in% ">300" & !(dat0$Plants_Age_years %in% page_set6) & dat0$Plants_Age_years %in% page_set6),
  which(dat0$Fragment_age_category %in% ">500" & !(dat0$Plants_Age_years %in% page_set6) & dat0$Plants_Age_years %in% page_set6),
  which(dat0$Fragment_age_category %in% ">1000" & !(dat0$Plants_Age_years %in% page_set6) & dat0$Plants_Age_years %in% page_set6)
)
stopifnot(is.na(temp[ids]))
temp[ids] <- 2 # >1 to <= 5 generations

ids <- c(
  which(dat0$Fragment_age_category %in% ">22" & !(dat0$Plants_Age_years %in% page_set4) & dat0$Plants_Age_years %in% page_set6),
  which(dat0$Fragment_age_category %in% "20-90" & !(dat0$Plants_Age_years %in% page_set6) & dat0$Plants_Age_years %in% page_set6),
  which(dat0$Fragment_age_category %in% ">30" & !(dat0$Plants_Age_years %in% page_set5) & dat0$Plants_Age_years %in% page_set6),
  which(dat0$Fragment_age_category %in% ">40" & !(dat0$Plants_Age_years %in% page_set6) & dat0$Plants_Age_years %in% page_set6),
  which(dat0$Fragment_age_category %in% "<49" & !(dat0$Plants_Age_years %in% page_set6) & dat0$Plants_Age_years %in% page_set6),
  which(dat0$Fragment_age_category %in% ">49" & !(dat0$Plants_Age_years %in% page_set6) & dat0$Plants_Age_years %in% page_set6),
  which(dat0$Fragment_age_category %in% "<50" & !(dat0$Plants_Age_years %in% page_set6) & dat0$Plants_Age_years %in% page_set6),
  which(dat0$Fragment_age_category %in% ">100" & !(dat0$Plants_Age_years %in% page_set6) & dat0$Plants_Age_years %in% page_set6),
  which(dat0$Fragment_age_category %in% ">300" & !(dat0$Plants_Age_years %in% page_set6) & dat0$Plants_Age_years %in% page_set6),
  which(dat0$Fragment_age_category %in% ">500" & !(dat0$Plants_Age_years %in% page_set6) & dat0$Plants_Age_years %in% page_set6),
  which(dat0$Fragment_age_category %in% ">1000" & !(dat0$Plants_Age_years %in% page_set6) & dat0$Plants_Age_years %in% page_set6)
)
stopifnot(is.na(temp[ids]))
temp[ids] <- 3 # > 5 generations

# Put into data.frame
dat0[, "NlifeSinceFrag_ord"] <- factor(temp,
  levels = seq_along(relage_categories), labels = relage_categories, ordered = TRUE)

# stopifnot(identical(
#   sum(table(dat0[, "FragmentAgeCat_ord_yrs"], dat0[, "SpeciesAgeCat_ord_yrs"])),
#   sum(!is.na(dat0[, "NlifeSinceFrag_ord"]))))

utemp <- unique(dat0[, c("ID_effect", "NlifeSinceFrag_ord")])
print(table(utemp[, "NlifeSinceFrag_ord"]))
print(table(dat0[, "SpeciesAgeCat_ord_yrs"], dat0[, "NlifeSinceFrag_ord"]))
print(table(dat0[, "FragmentAgeCat_ord_yrs"], dat0[, "NlifeSinceFrag_ord"]))

#--- Effect level value
dat0_sorted <- dat0[order(dat0$ID_effect), ]
ntimes_effect <- table(dat0_sorted$ID_effect)

temp <- with(dat0_sorted, tapply(NlifeSinceFrag_ord, ID_effect, function(f) {
  # nearest even order statistic
  quantile(as.integer(f), probs = 0.5, type = 3, na.rm = TRUE)
}))
temp <- factor(temp, levels = seq_along(relage_categories),
  labels = relage_categories, ordered = TRUE)
stopifnot(names(temp) == names(ntimes_effect))
dat0_sorted[, "Effect_NlifeSinceFrag_ord"] <- rep(temp, times = ntimes_effect)

# Restore original order
temp <- dat0$ID_record_orig
dat0 <- dat0_sorted[order(dat0_sorted$ID_record_orig), ]
stopifnot(temp == dat0$ID_record_orig)



#------------------------------------
#--- Species generation time relative to fragment age
# (~ maximum generations since fragmentation)
# replace missing data with 'NlifeSinceFrag_ord'
table(dat0[, "FragmentAgeCat_ord_yrs"], dat0[, "Animals_Generation_cat"], useNA = "always",
  deparse.level = 2)
relgen_categories <- c("(0, 1] Tf/Ng", "(1, 5] Tf/Ng", ">5 Tf/Ng") # Time(fragmentation) / N(generations)

temp <- rep(NA, nrow(dat0))

num_Generation_years <- as.numeric(dat0$Animals_Generation_years)
cats_Generation_years <- as.character(na.exclude(
  unique(dat0$Animals_Generation_years[is.na(num_Generation_years)])))
# [1] "<<1"

#   * Case 1: numeric fragment age && numeric generation time
num_relage1 <-  as.numeric(dat0$Fragment_age_yr) / num_Generation_years
num_relage1_cat <- cut(num_relage1, breaks = c(0, 1, 5, Inf), labels = relage_categories,
  ordered_result = TRUE)
ids <- !is.na(num_relage1_cat)
temp[ids] <- num_relage1_cat[ids]

#   * Case 2: ordinal fragment age && numeric generation time
num_hasnona <- !is.na(num_Generation_years)
ids <- c(
  which(dat0$Fragment_age_category %in% ">22" & num_hasnona & num_Generation_years <= 22),
  which(dat0$Fragment_age_category %in% "20-90" & num_hasnona & num_Generation_years <= 20),
  which(dat0$Fragment_age_category %in% ">30" & num_hasnona & num_Generation_years <= 30),
  which(dat0$Fragment_age_category %in% ">40" & num_hasnona & num_Generation_years <= 40),
  which(dat0$Fragment_age_category %in% "<49" & num_hasnona & num_Generation_years <= 49),
  which(dat0$Fragment_age_category %in% ">49" & num_hasnona & num_Generation_years <= 49),
  which(dat0$Fragment_age_category %in% "<50" & num_hasnona & num_Generation_years <= 50),
  which(dat0$Fragment_age_category %in% ">100" & num_hasnona & num_Generation_years <= 100),
  which(dat0$Fragment_age_category %in% ">300" & num_hasnona & num_Generation_years <= 300),
  which(dat0$Fragment_age_category %in% ">500" & num_hasnona & num_Generation_years <= 500),
  which(dat0$Fragment_age_category %in% ">1000" & num_hasnona & num_Generation_years <= 1000)
)
temp[ids] <- 1 # <= 1 generation

ids <- c(
  which(dat0$Fragment_age_category %in% ">22" &
    num_hasnona & num_Generation_years > 22 & num_Generation_years <= 5 * 22),
  which(dat0$Fragment_age_category %in% "20-90" &
    num_hasnona & num_Generation_years > 20 & num_Generation_years <= 5 * 90),
  which(dat0$Fragment_age_category %in% ">30" &
    num_hasnona & num_Generation_years > 30 & num_Generation_years <= 5 * 30),
  which(dat0$Fragment_age_category %in% ">40" &
    num_hasnona & num_Generation_years > 40 & num_Generation_years <= 5 * 40),
  which(dat0$Fragment_age_category %in% "<49" &
    num_hasnona & num_Generation_years > 49 & num_Generation_years <= 5 * 49),
  which(dat0$Fragment_age_category %in% ">49" &
    num_hasnona & num_Generation_years > 49 & num_Generation_years <= 5 * 49),
  which(dat0$Fragment_age_category %in% "<50" &
    num_hasnona & num_Generation_years > 50 & num_Generation_years <= 5 * 50),
  which(dat0$Fragment_age_category %in% ">100" &
    num_hasnona & num_Generation_years > 100 & num_Generation_years <= 5 * 100),
  which(dat0$Fragment_age_category %in% ">300" &
    num_hasnona & num_Generation_years > 300 & num_Generation_years <= 5 * 300),
  which(dat0$Fragment_age_category %in% ">500" &
    num_hasnona & num_Generation_years > 500 & num_Generation_years <= 5 * 500),
  which(dat0$Fragment_age_category %in% ">1000" &
    num_hasnona & num_Generation_years > 1000 & num_Generation_years <= 5 * 1000)
)
temp[ids] <- 2 # >1 to <= 5 generations

ids <- c(
  which(dat0$Fragment_age_category %in% ">22" & num_hasnona & num_Generation_years > 5 * 22),
  which(dat0$Fragment_age_category %in% "20-90" & num_hasnona & num_Generation_years > 5 * 90),
  which(dat0$Fragment_age_category %in% ">30" & num_hasnona & num_Generation_years > 5 * 30),
  which(dat0$Fragment_age_category %in% ">40" & num_hasnona & num_Generation_years > 5 * 40),
  which(dat0$Fragment_age_category %in% "<49" & num_hasnona & num_Generation_years > 5 * 49),
  which(dat0$Fragment_age_category %in% ">49" & num_hasnona & num_Generation_years > 5 * 49),
  which(dat0$Fragment_age_category %in% "<50" & num_hasnona & num_Generation_years > 5 * 50),
  which(dat0$Fragment_age_category %in% ">100" & num_hasnona & num_Generation_years > 5 * 100),
  which(dat0$Fragment_age_category %in% ">300" & num_hasnona & num_Generation_years > 5 * 300),
  which(dat0$Fragment_age_category %in% ">500" & num_hasnona & num_Generation_years > 5 * 500),
  which(dat0$Fragment_age_category %in% ">1000" & num_hasnona & num_Generation_years > 5 * 1000)
)
temp[ids] <- 3 # > 5 generations


#   * Case 3: numeric fragment age && ordinal generation time
num_hasnona <- !is.na(dat0$Fragment_age_yr)
ids <- c(
  which(dat0$Animals_Generation_years %in% "<<1" & num_hasnona & dat0$Fragment_age_yr <= 1)
)
temp[ids] <- 1 # <= 1 generation

ids <- c(
  which(dat0$Animals_Generation_years %in% "<<1" &
    num_hasnona & dat0$Fragment_age_yr > 1 & dat0$Fragment_age_yr <= 5 * 1)
)
temp[ids] <- 2 # >1 to <= 5 generations

ids <- c(
  which(dat0$Animals_Generation_years %in% "<<1" & num_hasnona & dat0$Fragment_age_yr > 5 * 1)
)
temp[ids] <- 3 # > 5 generations


#   * Case 4: ordinal fragment age && ordinal generation time
ids <- c(
  which(!is.na(dat0$Fragment_age_category) & dat0$Animals_Generation_years %in% "<<1")
)
temp[ids] <- 3 # > 5 generations



dat0[, "NGenSinceFrag_ord"] <- factor(temp,
  levels = seq_along(relgen_categories), labels = relgen_categories, ordered = TRUE)

stopifnot(identical(
  sum(table(dat0[, "FragmentAgeCat_ord_yrs"], dat0[, "NGenSinceFrag_ord"])),
  sum(!is.na(dat0[, "NGenSinceFrag_ord"]))))

# replace missing data with 'NlifeSinceFrag_ord'
print(table(dat0[, "NlifeSinceFrag_ord"], dat0[, "NGenSinceFrag_ord"]))

ids <- is.na(dat0[, "NGenSinceFrag_ord"])
dat0[ids, "NGenSinceFrag_ord"] <- factor(as.integer(dat0[ids, "NlifeSinceFrag_ord"]),
  levels = seq_along(relgen_categories), labels = relgen_categories, ordered = TRUE)

utemp <- unique(dat0[, c("ID_effect", "NGenSinceFrag_ord")])
print(table(utemp[, "NGenSinceFrag_ord"]))

utemp <- unique(dat0[, c("ID_effect", "NGenSinceFrag_ord", "NlifeSinceFrag_ord")])
print(table(utemp[, "NlifeSinceFrag_ord"], utemp[, "NGenSinceFrag_ord"]))

print(table(dat0[, "SpeciesAgeCat_ord_yrs"], dat0[, "NGenSinceFrag_ord"]))
print(table(dat0[, "FragmentAgeCat_ord_yrs"], dat0[, "NGenSinceFrag_ord"]))


#--- Effect level value
dat0_sorted <- dat0[order(dat0$ID_effect), ]
ntimes_effect <- table(dat0_sorted$ID_effect)

temp <- with(dat0_sorted, tapply(NGenSinceFrag_ord, ID_effect, function(f) {
  # nearest even order statistic
  quantile(as.integer(f), probs = 0.5, type = 3, na.rm = TRUE)
}))
temp <- factor(temp, levels = seq_along(relgen_categories),
  labels = relgen_categories, ordered = TRUE)
stopifnot(names(temp) == names(ntimes_effect))
dat0_sorted[, "Effect_NGenSinceFrag_ord"] <- rep(temp, times = ntimes_effect)

# Restore original order
temp <- dat0$ID_record_orig
dat0 <- dat0_sorted[order(dat0_sorted$ID_record_orig), ]
stopifnot(temp == dat0$ID_record_orig)


#------------------------------------
#--- Organism group
temp0 <- dat0$Organism_group1
print(table(temp0, useNA = "always"))
#      amphibia          ant          bat          bee       beetle         bird    bryophyte    butterfly
#            63           29           64           13           18           68           29           42
#    herbaceous  hymenoptera       lichen       mammal    marsupial   orthoptera      primate      reptile
#           563           12           16            8           59            6           32           51
#        rodent        shrub small mammal       spider         tree         <NA>
#            45          174            4           10          329            0

# Suggestions Dec 15, 2016:
#   - arthropods:
#       - one group
#       - flying vs. non-flying arthropods
#   - plants:
#       - woody vs. non-woody
#       - shrub vs. tree vs. herbaceous vs. (lichen & bryophytes)
#   - mammals:
#       - one group vs. bats
#       - top predators vs. herbivores vs. scavengers s.l. vs. bats
#   - herps
#       - one group
#       - amphibia vs. reptilia
#   - birds

# - arthropods
temp <- rep("non-arthropod", length(temp0))
temp[dat0$Class %in% c("Insecta", "Arachnida")] <- "arthropods"
dat0$Group_arthropods1 <- temp

temp <- rep("non-arthropod", length(temp0))
#TODO: group_arthropods_flying, group_arthropods_nonflying
group_arthropods <- c("ant", "bee", "beetle", "butterfly", "hymenoptera", "orthoptera", "spider")
group_arthropods_flying <- c()
group_arthropods_nonflying <- c()
temp[temp0 %in% group_arthropods_flying] <- "flying-arthropod"
temp[temp0 %in% group_arthropods_nonflying] <- "nonflying-arthropod"
dat0$Group_arthropods2 <- temp

# - plants
temp <- rep("non-plant", length(temp0))
group_plant <- dat0$Class %in% c("Jungermanniopsida", "Lecanoromycetes", "Liliopsida",
  "Magnoliopsida")
group_woody <- c("tree", "shrub")
group_nonwoody <- c("herbaceous", "bryophyte", "lichen")
temp[group_plant & temp0 %in% group_woody] <- "woody-plants"
temp[group_plant & temp0 %in% group_nonwoody] <- "nonwoody-plants"
dat0$Group_plants1 <- temp

temp <- rep("non-plant", length(temp0))
ids <- temp0 %in% c("tree", "shrub", "herbaceous")
group_lb <- c("bryophyte", "lichen")
temp[ids] <- temp0[ids]
temp[dat0$Class %in% c("Jungermanniopsida", "Lecanoromycetes") & temp0 %in% group_lb] <-
  "lichen-bryophytes"
dat0$Group_plants2 <- temp

# - mammals
temp <- rep("non-mammal", length(temp0))
group_mammals <- c("Chiroptera", "Rodentia", "Peramelemorphia", "Diprotodontia",
  "Didelphimorphia","Primates", "Carnivora", "Artiodactyla")
temp[dat0$Class %in% "Mammalia"] <- "mammals"
dat0$Group_mammals1 <- temp

temp <- rep("non-mammal", length(temp0))
#TODO: top predators vs. herbivores vs. scavengers s.l. vs. bats
group_bat <- c("bats")
group_toppredator <- c()
group_herbivore <- c()
group_scavenger <- c()
temp[temp0 %in% group_toppredator] <- "toppredator-mammal"
temp[temp0 %in% group_herbivore] <- "herbivore-mammal"
temp[temp0 %in% group_scavenger] <- "scavenger-mammal"
temp[temp0 %in% group_bat] <- "bat-mammal"
dat0$Group_mammals2 <- temp

# - herps
temp <- rep("non-herps", length(temp0))
temp[dat0$Class %in% c("Reptilia", "Amphibia")] <- "herps"
dat0$Group_herps1 <- temp

# - birds
temp <- rep("non-birds", length(temp0))
temp[dat0$Class %in% "Aves"] <- "birds"
dat0$Group_birds1 <- temp


# - combine
temp <- temp0
i_nonarthropod <- dat0$Group_arthropods1 == "non-arthropod"
i_nonplant <- dat0$Group_plants1 == "non-plant"
i_nonmammal <- dat0$Group_mammals1 == "non-mammal"
i_nonherps <- dat0$Group_herps1 == "non-herps"
i_nonbirds <- dat0$Group_birds1 == "non-birds"
ids <- !i_nonarthropod & i_nonplant & i_nonmammal & i_nonherps & i_nonbirds
temp[ids] <- dat0[ids, "Group_arthropods1"]
ids <- i_nonarthropod & !i_nonplant & i_nonmammal & i_nonherps & i_nonbirds
temp[ids] <- dat0[ids, "Group_plants1"]
ids <- i_nonarthropod & i_nonplant & !i_nonmammal & i_nonherps & i_nonbirds
temp[ids] <- dat0[ids, "Group_mammals1"]
ids <- i_nonarthropod & i_nonplant & i_nonmammal & !i_nonherps & i_nonbirds
temp[ids] <- dat0[ids, "Group_herps1"]
ids <- i_nonarthropod & i_nonplant & i_nonmammal & i_nonherps & !i_nonbirds
temp[ids] <- dat0[ids, "Group_birds1"]

dat0$Organism_group2 <- temp
stopifnot(!anyNA(dat0$Organism_group2))

utemp <- unique(dat0[, c("ID_effect", "Organism_group2")])
print(table(utemp$Organism_group2))
#     arthropods           birds           herps         mammals nonwoody-plants    woody-plants
#             14               9               9              26              44              48

print(table(dat0$Organism_group1, dat0$Organism_group2))



# Organism group 3: vertebrates vs. arthropods
temp <- dat0$Organism_group2
temp[!group_plant & i_nonarthropod] <- "vertebrates"
dat0$Organism_group3 <- temp
stopifnot(!anyNA(dat0$Organism_group3))

utemp <- unique(dat0[, c("ID_effect", "Organism_group3")])
print(table(utemp$Organism_group3))
#     arthropods nonwoody-plants     vertebrates    woody-plants
#             14              44              44              48

print(table(dat0$Organism_group1, dat0$Organism_group3))


# Organism group 4: woody vs. non-woody vs animals
temp <- dat0$Organism_group2
temp[!group_plant] <- "animals"
dat0$Organism_group4 <- temp
stopifnot(!anyNA(dat0$Organism_group4))

utemp <- unique(dat0[, c("ID_effect", "Organism_group4")])
print(table(utemp$Organism_group4))
#     animals nonwoody-plants     woody-plants
#          58              44               48

print(table(dat0$Organism_group1, dat0$Organism_group4))



# Organism group 5: woody vs. non-woody vs mammals/birds vs. other animals
temp <- dat0$Organism_group2
temp[temp %in% c("birds", "mammals")] <- "mammals&birds"
temp[temp %in% c("arthropods", "herps")] <- "arthropods&herps"
dat0$Organism_group5 <- temp
stopifnot(!anyNA(dat0$Organism_group5))

rtemp <- list(
  mA = !is.na(dat0[, "mA_standardized"]),
  HeSh = !is.na(dat0[, "Response_He"]) | !is.na(dat0[, "Response_Sh"]),
  PLp = !is.na(dat0[, "Response_PLp"]),
  Fis = !is.na(dat0[, "Response_Fis"]))
for (k in seq_along(rtemp)) {
  utemp <- unique(dat0[rtemp[[k]], c("ID_effect", "Organism_group5")])
  print(names(rtemp)[k])
  print(table(utemp$Organism_group5))
}


#------------------------------------
#--- Genetic marker types

print(sort(unique(dat0[, "Marker"])))
# [1] "AFLP"     "cDNA"     "Isozymes" "ISSR" "Microsat" "mtDNA"    "RAPD"

# cDNA should be cpDNA (and not complementary DNA)

# Distinguish categories:
#   1 = codominant DNA-based: microsatellites
#   2 = dominant DNA-based: aflp, RAPD, issr
#   3 = protein-based: isozymes
#   4 = uniparentally inherited DNA: cpDNA, mtDNA
marker_classes <- c("codominant", "dominant", "proteins", "uniparental")

temp <- sapply(dat0[, "Marker"], function(x)
  switch(x, AFLP = 2, cDNA = 4, Isozymes = 3, ISSR = 2, Microsat = 1, mtDNA = 4,
    RAPD = 2, NA))
dat0[, "Marker_cat"] <- factor(temp, levels = seq_along(marker_classes),
  labels = marker_classes, ordered = FALSE)



#------------------------------------
#--- STUDY SYSTEM
temp <- dat0$Habitat
print("Habitat types:")
print(sort(table(temp), decreasing = TRUE))

# Suggestion Dec 15, 2016: Ecoregions by Olson
# --> doesn't work because geographic coordinates are not study/plot locations

# Suggestion: Use forest types by Jaboury Ghazoul
# --> too many types
forests_jg <- c("Tropical rainforest", "Tropical montane cloud forest",
  "Tropical season forest", "Tropical dry forest and savannah woodlands",
  "Mediterranean woodlands", "Temperate deciduous forests", "Temperate rainforest",
  "Boreal forest")

# Sugestion Feb 20, 2017: pre-defined categories
habitat_types <- c(
  "Grasslands",
  "Scrublands and savannas",
  "Tropical moist and wet forests",
  "Temperate moist and wet forests",
  "Dry forests",
  "Boreal and alpine forests",
  "Other habitats")

if (FALSE) {
  # Prepare empty table to enter values for 'Habitat2' as extracted from articles
  write.csv(unique(dat0[, c("Paper_Nr", "ID_article", "FirstAuthor", "ID_effect",
    "Species_accepted", "Habitat", "Habitat2")]), file = fname_habitatreclass)
}

temp <- read.csv(fname_habitatreclass)

# copy 'Habitat2' to datafile 'dat0'
itemp <- match(dat0$ID_effect, temp$ID_effect, nomatch = 0)
stopifnot(sum(itemp > 0) == nrow(dat0))
dat0[, "Habitat2"] <- as.character(temp[itemp, "Habitat2"])

utemp <- unique(dat0[, c("Habitat", "Habitat2")])
utemp <- utemp[order(utemp[, "Habitat2"]), ]

utemp <- unique(dat0[, c("ID_effect", "Habitat", "Habitat2")])
print(as.matrix(table(utemp[, "Habitat2"])))
#  Dry forests                       17
#  Grasslands                        29
#  Other habitats                    10
#  Scrublands and savannas           21
#  Temperate moist and wet forests   38
#  Tropical moist and wet forests    40

temp <- with(utemp, table(Habitat, Habitat2))
# write.csv(temp, file.path(dir_in, "20170227_HabitatClassificationCrosstable.csv"))


# Sugestion March 27, 2017: pooled categories
temp <- rep(NA, nrow(dat0))
temp[dat0$Habitat2 %in% c("Dry forests", "Temperate moist and wet forests",
  "Tropical moist and wet forests")] <- "Forests"
temp[dat0$Habitat2 %in% c("Grasslands", "Other habitats", "Scrublands and savannas")] <-
  "Non-Forests"
dat0[, "Habitat3"] <- temp

utemp <- unique(dat0[, c("ID_effect", "Habitat", "Habitat3")])

print(as.matrix(table(utemp[, "Habitat3"])))
#  Forests                            95
#  Non-Forests                        60



#--- MATRIX
print("Matrix types:")
print(sort(table(dat0$Matrix)))

# Suggestions Dec 15, 2016:
#   - estimate distance between habitat and matrix
#   - sub-optimal former habitat vs. non-lethal non-habitat vs. lethal non-habitat (where lethality = f(taxon)
#   - some aggregate of ESA GlobCover neighbors

# --> doesn't work because geographic coordinates are not study/plot locations

# Suggestion Feb 6, 2017:pre-defined categories
matrix_types <- c(
  "ag-dominated",
  "tree-dominated",
  "urban/water-dominated")

temp2 <- rep(NA, dim(dat0)[1])

temp2[dat0$Matrix %in% c("agriculture", "agriculture/succession", "agriculture/tourism",
  "agriculture/logging/settlement", "agriculture/settlement/afforestation",
  "agriculture/logging", "agriculture/afforestation", "agriculture/settlement",
  "grassland", "grassland/scrubland/forest", "land-use intesification", "drainage")] <-
  "ag-dominated"

temp2[dat0$Matrix %in% c("secondary forest", "oil palm trees", "forest recovery",
  "logging/afforestation", "pine forest", "logging/human disturbance", "logging",
  "forest")] <- "tree-dominated"

temp2[dat0$Matrix %in% c("urbanization", "mining/tourism", "logging/settlement",
  "settlement", "water")] <- "urban/water-dominated"


# copy 'Matrix2' to datafile 'dat0'
dat0[, "Matrix2"] <- temp2

utemp <- unique(dat0[, c("ID_effect", "Matrix", "Matrix2")])

print(as.matrix(table(utemp[, "Matrix2"])))
#  ag-dominated           115
#  tree-dominated          19
#  urban/water-dominated   21

temp <- with(utemp, table(Matrix, Matrix2))
# write.csv(temp, file.path(dir_in, "20170224_MatrixClassificationCrosstable.csv"))


#--- Latitudinal bins
temp <- dat0$Effect_latitude
hist(temp)

temp <- rep(NA, nrow(dat0))
temp[dat0$Effect_latitude > -50 & dat0$Effect_latitude <= -25] <- 1
temp[dat0$Effect_latitude > -25 & dat0$Effect_latitude <= 0] <- 2
temp[dat0$Effect_latitude > 0 & dat0$Effect_latitude <= 25] <- 3
temp[dat0$Effect_latitude > 25 & dat0$Effect_latitude <= 50] <- 4
temp[dat0$Effect_latitude > 50 & dat0$Effect_latitude <= 75] <- 5

dat0[, "Effect_latitude_bin"] <- factor(temp, levels = 1:5,
  labels = c("25-50 S", "0-25 S", "0-25 N", "25-50 N", ">50 N"), ordered = TRUE)



#------------------------------------
#--- RESHAPE TO LONG FORMAT
icol_gen0 <- c("Ar", "mA", "Np", "H0", "He", "FST", "Fis", "PLp", "Sh")
coln_gen0 <- paste0("Response_", icol_gen0)
coln_gen0s <- paste0(icol_gen0, "_standardized")


dat0_final <- list()
temp <- names(dat0)[!(names(dat0) %in% c(coln_gen0, coln_gen0s))]
dat_nonmeasure <- dat0[, temp]

for (k in seq_along(icol_gen0)) {
  dat0_final[[k]] <- data.frame(dat_nonmeasure,
    Response = icol_gen0[k],
    Value = dat0[, coln_gen0[k]],
    is_standardized = dat0[, coln_gen0s[k]],
    stringsAsFactors = FALSE)
}

dat0_final <- do.call(rbind, dat0_final)


#------------------------------------
#--- SAVE CLEANED DATA
saveRDS(dat0, file = foutw)
saveRDS(dat0_final, file = foutl)
