
#--- PACKAGES
library("gdata")
library("stringr")
library("curl")
library("taxize")
library("maps")



#--- SETTINGS
range_N_is_standardized <- c(0.8, 1.2) # decision from meeting Dec 21, 2016

dir_prj <- "Prj04_GenetFragment_MetaAnalysis"
dir_in1 <- file.path(dir_prj, "2_Data", "1_Literature")
dir_in2 <- file.path(dir_prj, "2_Data", "4_DataExtracted")
dir_out <- file.path(dir_prj, "2_Data", "5_DataCleaned")
dir_pdfs <- file.path(dir_prj, "2_Data", "3_PDFs_IncludedArticles")

# File names
fin1 <- file.path(dir_in2, "170131-Data-fixed.xls")
fin2 <- file.path(dir_in1, "20170111_CompleteDB_AllRecords.csv")

if (!exists("tag_newdate")) tag_newdate <- format(Sys.Date(), "%Y%m%d")
fspeciesnames <- file.path(dir_in2, "20170130_speciesnames.rda")
fspeciesnames_new <- file.path(dir_in2, paste0(tag_newdate, "_speciesnames.rda"))
foutw <- file.path(dir_out, paste0(tag_newdate, "_ExtractedDataCleaned1_wide.rds"))

fname_fragsize_ranks <- file.path(dir_in2, "20171105_FragmentSize_Ranked.csv")
fname_fragsize_ranks_new <- file.path(dir_in2, paste0(tag_newdate, "_FragmentSize_Ranked.csv"))


# Categories
age_categories <- c("(0, 2] yrs", "(2-10] yrs", "(10-15] yrs", "(15-50] yrs",
  "(50-100] yrs", "(100-200] yrs", ">200 yrs")
age_categories2 <- c("(0, 10] yrs", "(10-50] yrs", ">50 yrs")


#--- READ DATA
print(paste(Sys.time(), "reading input", shQuote(basename(fin1))))
if (FALSE) {
  dat <- readxl::read_excel(fin1,
    col_types = c("numeric", "text", "text", "text", "text", "text", "text", "text",
      "numeric", "numeric", "numeric", "numeric", "text", "numeric", "text", "text",
      "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
      "numeric", "numeric", "text", "numeric", "numeric", "numeric", "numeric", "text",
      "text"))

  dat <- openxlsx::read.xlsx(fin1)
}

dat_colnames <- names(gdata::read.xls(fin1, as.is = TRUE, nrows = 1))
dat <- gdata::read.xls(fin1, as.is = TRUE, skip = 3, header = FALSE)
names(dat) <- dat_colnames
dat[dat == ""] <- NA

datr <- read.csv(fin2)


#------------------------------------
#--- REMOVE ADDITIONAL STUDIES
# Dec 19, 2017: we exclude mosses and lichens:
temp <- dat[, "Organism.group"] %in% c("bryophyte", "lichen")
to_remove_ID_article <- unique(dat[temp, "Study"])
#   -> bryophyte: remove study '14-15-7321-Pohjamoetal2008.pdf'
#   -> lichen: remove study '26-7369-Otaloraetal2011.pdf'
dat <- dat[!(dat$Study %in% to_remove_ID_article), ]

#------------------------------------
#--- CLEAN DATA
dat0 <- list()
dat0$ID_record_orig <- seq_len(dim(dat)[1])
dat0$ID_article <- as.integer(dat$Study)
dat0$FirstAuthor <- as.character(dat$Study.1)
dat0 <- as.data.frame(dat0, stringsAsFactors = FALSE)


#------------------------------------
#--- PUBLICATION DATA
ids <- match(dat0$ID_article, datr$RecordID, nomatch = 0)
dat0[ids > 0, "Pub_Year"] <- datr[ids, "Year"]


#------------------------------------
#--- Extract "Paper_Nr" by Bruno and Hans-Peter and link to 'ID_article'
temp <- list.files(dir_pdfs)
temp <- sapply(strsplit(temp, split = "-", fixed = TRUE), function(x) {

    x2 <- as.integer(x[2])
    if (is.na(x2)) {
      # Filename format: PaperNR-ArticleID_FirstAuthor_Year
      x <- c(x[1], as.integer(unlist(strsplit(x[2], split = "_", fixed = TRUE)[[1]])), x[3])
      x2 <- as.integer(x[2])
    }

    x3 <- as.integer(x[3])
    res <- if (!is.na(x3)) {
      # Filename format: PaperNR1-PaperNR2-ArticleID
      matrix(c(x[1], x2, x3, x3), nrow = 2, ncol = 2)
    } else {
      matrix(c(x[1], x2), nrow = 1, ncol = 2)
    }

    res
  })
temp <- do.call(rbind, temp)

ids <- data.frame(pnr = temp[, 1], ida = temp[, 2], stringsAsFactors = FALSE)

repeat {
  itemp <- which(is.na(ids$ida))
  if (length(itemp) == 0)
    break

  ids[itemp, "ida"] <- ids[itemp - 1, "ida"]
}

ids2 <- aggregate(ids$pnr, by = list(ids$ida), paste0, collapse = "/")
colnames(ids2) <- c("ida", "pnr")

id_x <- match(dat0[, "ID_article"], ids2[, "ida"], nomatch = 0)
use_r <- id_x > 0
# Missing PDFs, article_ID: 471  875 1836 6530 6632 7321
print(sort(unique(dat0[!use_r, "ID_article"])))

dat0[, "Paper_Nr"] <- NA
dat0[use_r, "Paper_Nr"] <- ids2[id_x, "pnr"]


#--- LOCATION
dat0$Fragment_longitude <- as.numeric(dat$Coordinates.1)
dat0$Fragment_latitude <- as.numeric(dat$Coordinates)

#- Fix country names
id <- which(dat[, "Country"] %in% "USA, Oregon")
dat[id, "Country"] <- "USA"

#- Plot
zooms <- list(
  Global = list(xlim = NULL, ylim = NULL),
  Europe = list(xlim = c(-7, 25), ylim = c(35, 65)),
  Australia = list(xlim = c(110, 160), ylim = c(-40, -10)),
  MesoAmerica = list(xlim = c(-105, -80), ylim = c(5, 20)),
  NorthAmerica = list(xlim = c(-130, -60), ylim = c(30, 50)),
  SouthAmerica = list(xlim = c(-90, -25), ylim = c(-40, 5))
)

for (k in seq_along(zooms)) {
  ftemp <- file.path(dir_out, "Figures", paste0("Fig_Locations_", names(zooms)[k], ".pdf"))

  pdf(height = 2 * 14, width = 20, file = ftemp)
    temp <- if (all(sapply(zooms[[k]], is.null))) {
        !duplicated(dat0$ID_article)
      } else {
        rep(TRUE, dim(dat0)[1])
      }
    lab_country <- factor(dat[, "Country"])
    cols <- rainbow(n = nlevels(lab_country))
    old_op <- par(mfrow = c(2, 1))

    for (p in 1:2) {
      with(dat0, plot(Fragment_longitude, Fragment_latitude,
        xlim = zooms[[k]][["xlim"]], ylim = zooms[[k]][["ylim"]],
        pch = 16, col = cols[lab_country]))
      map("world", col = "gray40", add = TRUE)

      if (p == 1) {
        with(dat0[temp, ], text(Fragment_longitude, Fragment_latitude, pos = 3,
          labels = lab_country[temp]))
      } else if (p == 2) {
        with(dat0[temp, ], text(Fragment_longitude, Fragment_latitude, pos = 3,
          labels = Paper_Nr))
      }
    }

    par(old_op)
  dev.off()

}



#--- SPECIES NAME
do_speciesnames <- TRUE
test_speciesnames <- NULL

temp <- stringr::str_trim(dat$Species) # clean spaces at beginning or end
temp <- gsub("[[:space:]]+", " ", temp) # clean multiple spaces between words
dat0$Species <- temp
specieslist <- unique(dat0$Species)

if (file.exists(fspeciesnames)) {
  test_speciesnames <- new.env(parent = emptyenv())
  load(fspeciesnames, envir = test_speciesnames)
  do_speciesnames <- !all(specieslist %in% test_speciesnames$specieslist)
}

if (!do_speciesnames) {
  ids <- match(specieslist, test_speciesnames$specieslist, nomatch = 0)

  if (identical(specieslist, test_speciesnames$specieslist[ids])) {
    species_resolved <- get("species_resolved", envir = test_speciesnames)[ids, ]
    unresolved_species <- get("unresolved_species", envir = test_speciesnames)[ids]
    species_accepted <- get("species_accepted", envir = test_speciesnames)[ids]

  } else {
    do_speciesnames <- TRUE
  }
}

if (do_speciesnames &&
  (requireNamespace("curl") && curl::has_internet()) || !requireNamespace("curl")) {

  # Resolve species names based on http://gni.globalnames.org/
  species_resolved <- taxize::gnr_resolve(names = specieslist, canonical = TRUE,
    best_match_only = TRUE, cap_first = TRUE, fields = "minimal")
  unresolved_species <- attr(species_resolved, "not_known")

  # Find synonyms and extract accepted names from ITIS
  species_resolved[, "species_accepted"] <- species_resolved[, "db_accepted"] <- NA
  species_accepted <- list()

  needs_accepted <- species_resolved[is.na(species_resolved[, "species_accepted"]),
    "matched_name2"]
  species_accepted[["ITIS"]] <- taxize::synonyms(x = needs_accepted, db = "itis",
    accepted = FALSE)

  for (k in seq_along(species_accepted[["ITIS"]])) {
    res <- species_accepted[["ITIS"]][[k]]
    irow <- which(species_resolved[, "matched_name2"] == needs_accepted[k])

    if (length(irow) == 1 && length(res) > 0 && !is.na(res)) {
      id_res <- names(res)

      x <- if ("acc_name" %in% id_res) {
        # if 'acc_name' exists, then this is the accepted name by ITIS
        res[1, "acc_name"]

      } else if ("message" %in% id_res && identical(res[1, "message"], "no syns found")) {
        # if 'message' exists and contains "no syns found", then queried name is accepted
        species_resolved[irow, "matched_name2"]

      } else NA

      species_resolved[irow, "species_accepted"] <- x
      species_resolved[irow, "db_accepted"] <- if (is.na(x)) NA else "ITIS"
    }
  }

  # Find synonyms and extract accepted names from Catalogue of Life
  needs_accepted <- species_resolved[is.na(species_resolved[, "species_accepted"]),
    "matched_name2"]
  species_accepted[["CoL"]] <- taxize::synonyms(x = needs_accepted, db = "col", rows = 1)

  for (k in seq_along(species_accepted[["CoL"]])) {
    res <- species_accepted[["CoL"]][[k]]
    irow <- which(species_resolved[, "matched_name2"] == needs_accepted[k])

    if (length(irow) == 1 && length(res) > 0 && !is.na(res) &&
      !is.null(res[["name_status"]])) {
      # if 'name_status' contains "synonyms", then queried name is accepted by CoL
      has_acc <- NROW(res) > 0 && res$name_status == "synonym"

      if (has_acc) {
        species_resolved[irow, "species_accepted"] <- species_resolved[irow, "matched_name2"]
        species_resolved[irow, "db_accepted"] <- "Catalogue of Life"
      }
    }
  }

  # Extract accepted names from Tropicos
  needs_accepted <- species_resolved[is.na(species_resolved[, "species_accepted"]),
    "matched_name2"]
  tropicos <- list()

  for (k in seq_along(needs_accepted)) {
    res <- taxize::tp_search(name = needs_accepted[k])
    irow <- which(species_resolved[, "matched_name2"] == needs_accepted[k])
    tpan <- NA

    if (length(irow) == 1 && length(res) > 0 && !is.na(res)) {
      if (all(c("nomenclaturestatusname", "rankabbreviation", "scientificname") %in%
        names(res))) {
        temp <- res[, "nomenclaturestatusname"] %in% "Legitimate"
        has_legitimate <- NROW(res) > 0 && any(temp)

        if (has_legitimate) {
          temp <- res[temp & res[, "rankabbreviation"] == "sp.", "scientificname"]

          if (!is.na(temp) && NROW(temp) == 1 && is.character(temp)) {
            species_resolved[irow, "species_accepted"] <- temp
            species_resolved[irow, "db_accepted"] <- "Tropicos"
          }
        }
      }

      if (is.na(species_resolved[irow, "species_accepted"])) {
        if (all(c("nameid") %in% names(res))) {
          tpan <- taxize::tp_accnames(id = res[2, "nameid"])

          if (!("Error" %in% names(tpan))) {
            print(paste(k, needs_accepted[k], "'http://services.tropicos.org/Name/<id>/AcceptedNames'",
              "works again: implement code"))
            print(tpan)

          } else {
            print(tpan)
          }

        } else {
          print(paste(k, needs_accepted[k], ": Tropicos' 'nameid' not found"))
        }
      }
    }

    tropicos[[needs_accepted[k]]] <- list(tps = res, tpan = tpan)
  }

  species_accepted[["Tropicos"]] <- tropicos


  # Find accepted names: missing
  needs_accepted <- species_resolved[is.na(species_resolved[, "species_accepted"]),
    "matched_name2"]
  if (length(needs_accepted) > 0)
    print(paste("Following species names are resolved, but no accepted name was found:",
      paste(needs_accepted, collapse = ", ")))

  # Save to disk
  save(specieslist, species_resolved, unresolved_species, species_accepted,
    file = fspeciesnames_new)
  if (!identical(fspeciesnames, fspeciesnames_new))
    unlink(fspeciesnames)
}


# Resolved names
irow_corrected <- apply(species_resolved, 1,
  function(x) x["user_supplied_name"] != x["matched_name2"])
irow_corrected <- which(irow_corrected)

print(species_resolved[irow_corrected, c("user_supplied_name", "data_source_title",
  "matched_name2")])
#                    user_supplied_name      data_source_title                  matched_name2
#104     Medicago sativa subsp. falcata                    EOL        Medicago sativa falcata
#105 Leucochrysum albicans var tricolor GBIF Backbone Taxonomy Leucochrysum albicans tricolor
#113       Ardisia crenata var. bicolor                   NCBI        Ardisia crenata bicolor
#128                   Trilium reliquum      Catalogue of Life              Trillium reliquum
#135                     Stippa pennata      Catalogue of Life                  Stipa pennata

print(unresolved_species)

dat0$Species_resolved <- dat0$Species
dat0$Species_source <- rep("publication", dim(dat0)[1])

for (k in seq_along(irow_corrected)) {
  i <- irow_corrected[k]
  irows <- which(species_resolved[i, "user_supplied_name"] == dat0$Species)

  dat0$Species_source[irows] <- species_resolved[i, "data_source_title"]
  dat0$Species_resolved[irows] <- species_resolved[i, "matched_name2"]
}

# Accepted names
dat0$Species_accepted <- dat0$Accepted_source <- rep(NA, dim(dat0)[1])

for (k in seq_len(dim(dat0)[1])) {
  irows <- which(species_resolved[, "user_supplied_name"] == dat0$Species[k])

  dat0$Species_accepted[k] <- species_resolved[irows, "species_accepted"]
  dat0$Accepted_source[k] <- species_resolved[irows, "db_accepted"]
}




#--- Organism group
temp0 <- tolower(dat$Organism.group)
# Dec 5, 2016: Bruno/Hans-Peter: three groups: 'shrub', 'tree', and 'vascular plant' (other)
temp0[temp0 == "plant"] <- "vascular plant"
temp0[temp0 == "vascular plant"] <- "herbaceous"
dat0$Organism_group1 <- temp0


#--- UNITS OF ANALYSIS
# Calculate ID for each unit ([article x country] x species)
# data.frame 'dat0' contains rows for each unit x fragment
temp <- ceiling(log10(1 + dim(dat0)[1]))
flag_article <- formatC(dat0$ID_article, width = temp, flag = "0")
temp1 <- gsub(" ", "_", dat$Country)
temp2 <- gsub(" ", "_", dat0$Species_resolved)
dat0$ID_unit <- paste(flag_article, temp1, temp2, sep = "_")

# Fragments
# Articles which share fragment names: identical or different fragments?
temp <- table(dat0$ID_article, dat$Fragment)
temp <- apply(temp, 2, function(x) sum(x > 0))
temp <- names(temp)[temp > 1]

for (k in temp) {
  icols <- which(k == dat$Fragment)
  x <- cbind(dat0[icols, c("ID_article", "FirstAuthor", "Species_resolved",
    "Fragment_longitude", "Fragment_latitude")],
    dat[icols, c("Fragment", "Fragment_No")])

  # Fragements in similar location (e.g. < 100 km apart)
  x_dists_km <- sp::spDists(as.matrix(x[, c("Fragment_longitude", "Fragment_latitude")]),
    longlat = TRUE)

  lt <- lower.tri(x_dists_km)
  ids_near <- which(lt, arr.ind = TRUE)[which(x_dists_km[lt] < 100),]

  if (NROW(ids_near) > 0) {
    cat(paste("Fragment descriptor", shQuote(k), "\n"))
    xx <- x[sort(unique(as.vector(ids_near))), ]
    print(xx[do.call(order, xx[, c("Fragment_longitude", "Fragment_latitude")]), ])
    cat("\n")
  }

  # Fragment descriptor 'Roche'
  #    ID_article FirstAuthor     Species_resolved Fragment Fragment_No
  # 667       7336      Honnay Anthyllis vulneraria    Roche           3
  # 859       7133   Jacquemyn       Cirsium acaule    Roche          11
  #    Coordinates Coordinates.1
  # 667    50.07253      4.539539
  # 859    50.09092      4.595247

  # Fragment descriptor 'GL'
  #      ID_article FirstAuthor   Species_resolved
  # 696        7128          Li   Isoodon obesulus
  # 1315       7017    Malekian Petaurus breviceps
  #      Fragment_longitude Fragment_latitude Fragment Fragment_No
  # 696            140.5333         -37.63333       GL          14
  # 1315           140.7083         -37.75000       GL           9

  # Fragment descriptor 'C'
  #      ID_article FirstAuthor    Species_resolved
  # 1000       7326        Noel  Plethodon cinereus
  # 249        7545       Rogic Peromyscus leucopus
  #      Fragment_longitude Fragment_latitude Fragment Fragment_No
  # 1000          -73.61620          45.50540        C           2
  # 249           -73.23100          45.48183        C           3

  # Fragment descriptor 'A'
  #     ID_article FirstAuthor    Species_resolved
  # 999       7326        Noel  Plethodon cinereus
  # 247       7545       Rogic Peromyscus leucopus
  #     Fragment_longitude Fragment_latitude Fragment Fragment_No
  # 999           -73.6162          45.50540        A           1
  # 247           -73.2310          45.48183        A           1

  # Fragment descriptor '9'
  #      ID_article FirstAuthor       Species_resolved
  # 354        7521   Heelemann     Hemimeris racemosa
  # 372        7521   Heelemann Eriocephalus africanus
  # 1388       7391   Heelemann        Nemesia barbata
  #      Fragment_longitude Fragment_latitude Fragment Fragment_No
  # 354            18.66222         -33.76389        9           9
  # 372            18.66222         -33.76389        9           9
  # 1388           18.76611         -33.91472        9           8

  }

# Impossible to tell based on extracted data ==> assume each fragment x article is a
# unique combination
dat0$ID_fragment <- paste0(dat$Fragment, "_art", flag_article)


#--- STUDY SYSTEM
temp <- stringr::str_trim(dat$Habitat)
temp[temp == "tropical rain forest"] <- "tropical rainforest"
temp <- tolower(temp)
sort(unique(temp))
dat0$Habitat <- temp

temp <- stringr::str_trim(dat$Matrix)
temp[temp == "secondary forestes"] <- "secondary forest"
temp <- tolower(temp)
sort(unique(temp))
dat0$Matrix <- temp


#--- FRAGMENTS
# Fragment type
dat0$is_fragment <- dat$Control == 0

# Fragmentation age
temp <- dat$Age_history
dat0$Fragment_age_yr <- as.numeric(temp)
dat0$Fragment_age_source <- rep("publication", dim(dat0)[1])
temp <- ifelse(is.na(dat0$Fragment_age_yr) & !is.na(temp), temp, NA)
dat0$Fragment_age_category <- gsub(" ", "", temp)

ina <- is.na(dat0$Fragment_age_yr)
if (any(ina)) {
  temp <- dat$Assumed.age
  dat0$Fragment_age_yr[ina] <- as.numeric(temp[ina])
  iest <- ina & !is.na(temp)
  dat0$Fragment_age_source[iest] <- "estimated"

  ina <- is.na(dat0$Fragment_age_category)
  temp <- ifelse(is.na(dat0$Fragment_age_yr) & !is.na(temp), temp, NA)
  dat0$Fragment_age_category[ina] <- gsub(" ", "", temp[ina])
}

# Fragment area
temp <- dat$Area
dat0$Fragment_area_ha <- as.numeric(temp) # unit = hectare
temp <- ifelse(is.na(dat0$Fragment_area_ha) & !is.na(temp), temp, NA)
temp <- gsub(" ", "", temp)
dat0$Fragment_area_category <- ifelse(!is.na(temp), paste(temp, "hectare"), NA)

# Population size
temp <- dat$Population.size
dat0$Pop_size_n <- as.integer(temp)
temp <- ifelse(is.na(dat0$Pop_size_n) & !is.na(temp), temp, NA)
temp <- gsub(" ", "", temp)
dat0$Pop_size_category <- with(dat, ifelse(is.na(Population.size.1), temp,
  ifelse(is.na(temp), Population.size.1, paste(Population.size.1, "=", temp, "n"))))

# Population density
dat0$Pop_density_nPERha <- 1e4 * dat$Population.density # unit of 'Population.density' = # m-2


#--- GENETIC DATA
dat0$Marker <- dat$Marker.type

# N of examined individuals
temp <- dat$No.Indivduals.examined
dat0$Ind_sampled_nPERfragment <- as.numeric(temp)
temp <- ifelse(is.na(dat0$Ind_sampled_nPERfragment) & !is.na(temp), temp, NA)
dat0$Ind_sampled_category <- gsub(" ", "", temp)

# Genetic indices
icol_gen <- c("NumberAlleles", "NumberAllelesPerLocus", "Np", "H0", "He", "FST", "Fis",
  "X.PolymLoci", "Genetic.diversity")
icol_gen0 <- c("Ar", "mA", "Np", "H0", "He", "FST", "Fis", "PLp", "Sh")
coln_gen0 <- paste0("Response_", icol_gen0)
stopifnot(length(icol_gen) == length(icol_gen0))
temp <- dat[, icol_gen]

gen_notnumeric <- which(!sapply(temp, is.numeric))

for (k in gen_notnumeric) {
  temp[, k] <- as.numeric(temp[, k])
}

for (k in seq_along(icol_gen)) {
  dat0[[coln_gen0[k]]] <- temp[, icol_gen[k]]
}


#--- UNITS OF ANALYSIS (part 2)
# One effect size for each article x species x marker
dat0$ID_effect <- apply(dat0[, c("ID_unit", "Marker")], 1, paste, collapse = "-")
dat0_sorted <- dat0[order(dat0$ID_effect), ]

ntimes_effect <- table(dat0_sorted$ID_effect)
#ntimes_orig <- tapply(as.integer(dat$Fragment_No), dat0_sorted$ID_effect, length)
#ntimes_orig <- tapply(as.integer(dat$Fragment_No), dat0_sorted$ID_effect, max, na.rm = TRUE)
#cbind(as.integer(ntimes_effect), as.integer(ntimes_orig))
#stopifnot(as.integer(ntimes_effect), as.integer(ntimes_orig))


#--- AGGREGATE PER EFFECT
# Location (part 2)
# Because Bruno and Hans-Peter entered fragment coordinates as study averages,
# effect-aggregated coordinates with variances > 0 suggest data entry mistakes
temp <- with(dat0_sorted, tapply(Fragment_longitude, ID_effect, var, na.rm = TRUE))
if (any(temp > 0)) {
  stop(paste("Following effects have more than one value of longitude:\n",
    paste("*", names(temp)[temp > 0], collapse = "\n")))
}
temp <- with(dat0_sorted, tapply(Fragment_latitude, ID_effect, var, na.rm = TRUE))
if (any(temp > 0)) {
  stop(paste("Following effects have more than one value of latitude:\n",
    paste("*", names(temp)[temp > 0], collapse = "\n")))
}

temp <- with(dat0_sorted, tapply(Fragment_longitude, ID_effect, mean, na.rm = TRUE))
stopifnot(names(temp) == names(ntimes_effect))
dat0_sorted$Effect_longitude <- rep(temp, times = ntimes_effect)

temp <- with(dat0_sorted, tapply(Fragment_latitude, ID_effect, mean, na.rm = TRUE))
stopifnot(names(temp) == names(ntimes_effect))
dat0_sorted$Effect_latitude <- rep(temp, times = ntimes_effect)


# Fragmentation age (part 2)
temp <- with(dat0_sorted, tapply(Fragment_age_yr, ID_effect, mean, na.rm = TRUE))
stopifnot(names(temp) == names(ntimes_effect))
dat0_sorted$Effect_mean_age_yr <- rep(temp, times = ntimes_effect)


temp <- with(dat0_sorted, tapply(Fragment_age_yr, ID_effect, function(x) diff(range(x, na.rm = TRUE))))
temp[!is.finite(temp)] <- NA
stopifnot(names(temp) == names(ntimes_effect))
dat0_sorted$Effect_range_age_yr <- rep(temp, times = ntimes_effect)


#--- Fragment age (ordinal)
table(dat0_sorted$Fragment_age_yr)
table(dat0_sorted$Fragment_age_category)
#  <49   <50  >100 >1000   >22   >30  >300   >40   >49  >500 20-90
#    1     1    31    15     1    12    44     4     3     4     7

temp <- rep(NA, nrow(dat0_sorted))
temp[dat0_sorted$Fragment_age_yr <= 2] <- 1 # 1-2 yr
temp[dat0_sorted$Fragment_age_yr %in% 3:10] <- 2 # 3-10 yrs
temp[dat0_sorted$Fragment_age_yr %in% 11:15] <- 3 # 11-15 yrs
temp[dat0_sorted$Fragment_age_yr %in% 16:50 |
  dat0_sorted$Fragment_age_category %in% c("<49", "<50", ">22", ">30", ">40")] <- 4 # 16-50 yrs
temp[dat0_sorted$Fragment_age_yr %in% 51:100 |
  dat0_sorted$Fragment_age_category %in% c(">49", "20-90")] <- 5 # 51-100 yrs
temp[dat0_sorted$Fragment_age_yr %in% 101:200 |
  dat0_sorted$Fragment_age_category %in% c(">100")] <- 6 # 101-200 yrs
temp[dat0_sorted$Fragment_age_yr > 200 |
  dat0_sorted$Fragment_age_category %in% c(">1000", ">300", ">500")] <- 7 # >200 yrs
dat0_sorted[, "FragmentAgeCat_ord_yrs"] <- factor(temp,
    levels = seq_along(age_categories), labels = age_categories, ordered = TRUE)

table(dat0_sorted$Fragment_age_yr, dat0_sorted[, "FragmentAgeCat_ord_yrs"], useNA = "always")
table(dat0_sorted$Fragment_age_category, dat0_sorted[, "FragmentAgeCat_ord_yrs"], useNA = "always")

temp <- with(dat0_sorted, tapply(FragmentAgeCat_ord_yrs, ID_effect, function(f) {
  # nearest even order statistic
  quantile(as.integer(f), probs = 0.5, type = 3, na.rm = TRUE)
}))
temp <- factor(temp, levels = seq_along(age_categories), labels = age_categories, ordered = TRUE)
stopifnot(names(temp) == names(ntimes_effect))
dat0_sorted[, "EffectAgeCat_ord_yrs"] <- rep(temp, times = ntimes_effect)


# Aggregated age categories
temp <- rep(NA, nrow(dat0_sorted))
temp[dat0_sorted[, "EffectAgeCat_ord_yrs"] %in% c("(0, 2] yrs", "(2-10] yrs")] <- 1
temp[dat0_sorted[, "EffectAgeCat_ord_yrs"] %in% c("(10-15] yrs", "(15-50] yrs")] <- 2
temp[dat0_sorted[, "EffectAgeCat_ord_yrs"] %in% c("(50-100] yrs", "(100-200] yrs", ">200 yrs")] <- 3
dat0_sorted[, "EffectAgeCat2_ord_yrs"] <- factor(temp, levels = seq_along(age_categories2),
  labels = age_categories2, ordered = TRUE)




#--- RANK PER EFFECT
#--- Combine numeric and ordinal input data for: 'fragment area' and 'population size'
# Note: these ranks are only used for correlation calculations per effect, i.e., ranking
#   is ok to occur only within effect because it does not need to be comparable among effects

# First attempt at ranking (ignoring ordinal data)
ID_effect <- dat0_sorted$ID_effect
for (k in seq_along(ID_effect)) {
  itemp <- dat0_sorted$ID_effect == ID_effect[k]
  dtemp <- dat0_sorted[itemp, ]

  dat0_sorted[itemp, "Fragment_area_rank"] <- rank(dtemp[, "Fragment_area_ha"],
    na.last = "keep")
  dat0_sorted[itemp, "Pop_size_rank"] <- rank(dtemp[, "Pop_size_n"], na.last = "keep")
}


# Export for manual integration of ordinal data
if (!file.exists(fname_fragsize_ranks)) {
  rtemp <- dat0_sorted[, c("ID_record_orig", "ID_article", "FirstAuthor", "ID_effect",
    "ID_fragment", "Ind_sampled_nPERfragment", "Ind_sampled_category",
    "Fragment_area_ha", "Fragment_area_category", "Fragment_area_rank",
    "Pop_size_n", "Pop_size_category", "Pop_size_rank")]
  rtemp$rowID <- seq_len(nrow(rtemp))
  write.csv(rtemp, file = fname_fragsize_ranks_new, row.names = FALSE)

} else {
  rtemp <- read.csv(fname_fragsize_ranks)
}

if (all(dat0_sorted$ID_record_orig %in% rtemp$ID_record_orig)) {
  ids <- match(dat0_sorted$ID_record_orig, rtemp$ID_record_orig, nomatch = 0)

  stopifnot(identical(dat0_sorted$ID_record_orig[ids > 0], rtemp$ID_record_orig[ids]))

  dat0_sorted[ids > 0, "Fragment_area_rank"] <- rtemp[ids, "Fragment_area_rank"]
  dat0_sorted[ids > 0, "Pop_size_rank"] <- rtemp[ids, "Pop_size_rank"]

  temp1 <- rtemp[ids, "Ind_sampled_nPERfragment"]
  temp2 <- rtemp[ids, "Ind_sampled_nPERfragment_drs"]
  dat0_sorted[ids > 0, "Ind_sampled_nPERfragment"] <- ifelse(is.na(temp2), temp1, temp2)

} else {
    stop("'fname_fragsize_ranks$ID_record_orig' and 'dat0_sorted$ID_record_orig' disagree")
}





#----------------------------
#--- Restore original order
temp <- dat0$ID_record_orig
dat0 <- dat0_sorted[order(dat0_sorted$ID_record_orig), ]
stopifnot(temp == dat0$ID_record_orig)


#--- COMMENTS
# 20161205: Bruno & Hans-Peter: "each comment per ArticleID applies to every row of that ArticleID"
irow_order <- order(dat0$ID_article)
dat0_sorted <- dat0[irow_order, ]
coms1 <- dat$Remarks[irow_order]

temp <- tapply(coms1, dat0_sorted$ID_article, function(x) paste(unique(x), collapse = "; "))
ntimes_article <- table(dat0_sorted$ID_article)
stopifnot(names(temp) == names(ntimes_article))
coms_unip <- rep(temp, times = ntimes_article)

stopifnot(lengths(tapply(coms_unip, dat0_sorted$ID_effect, unique)) == 1)


#--- Are genetic measures standardized by sampling size/effort?
# 'all data standardized' is true if all samples per effect had (almost-)identical numbers
# of individuals
irow_orig <- seq_len(dim(dat0)[1])
irow_order2 <- order(dat0_sorted$ID_effect)
dat0_sorted2 <- dat0_sorted[irow_order2, ]
has_nearconstN <- tapply(dat0_sorted2$Ind_sampled_nPERfragment, dat0_sorted2$ID_effect,
  function(x) {
    x <- na.exclude(x)
    mtemp <- mean(x)
    all(is.finite(mtemp), x > range_N_is_standardized[1] * mtemp,
      x < range_N_is_standardized[2] * mtemp)
  })
temp <- tapply(coms_unip[irow_order2], dat0_sorted2$ID_effect, function(x)
  paste(unique(x), collapse = "; "))
temp[has_nearconstN] <- paste(temp[has_nearconstN], "all data standardized", sep = "; ")

ntimes_effect <- table(dat0_sorted2$ID_effect)
temp <- rep(temp, times = ntimes_effect)
coms_unip <- temp[order(irow_orig[irow_order2])]

# init
template <- rep(NA, dim(dat0_sorted)[1])
coln_gen0s <- paste0(icol_gen0, "_standardized")

tag_comment <- function(comments, dat, tag, val) {
  temp <- lapply(tag, grep, x = comments)
  temp <- unique(unlist(temp))
  dat[temp] <- val
  dat
}

# Allelic richness; mean allelic richness (A per locus)
temp1 <- tag_comment(coms_unip, template, c("data not corrected"), FALSE)
temp1 <- tag_comment(coms_unip, temp1, c("all data standardized", "allelic richness standardized",
  "rarified allelic richness", "standardized allelic richness"), TRUE)
dat0_sorted$Ar_standardized <- ifelse(is.na(dat0_sorted$Response_Ar), NA, temp1)
dat0_sorted$mA_standardized <- ifelse(is.na(dat0_sorted$Response_mA), NA, temp1)

# Polymorphic loci
temp1 <- tag_comment(coms_unip, template, c("data not corrected"), FALSE)
temp1 <- tag_comment(coms_unip, temp1, c("all data standardized", "rarified polymorphic loci"), TRUE)
dat0_sorted$P_standardized <- ifelse(is.na(dat0_sorted$Response_P), NA, temp1)

# He
temp1 <- tag_comment(coms_unip, template, c("data not corrected"), FALSE)
temp1 <- tag_comment(coms_unip, temp1, c("all data standardized", "rarified He"), TRUE)
dat0_sorted$He_standardized <- ifelse(is.na(dat0_sorted$Response_He), NA, temp1)

# Other icol_gen0
for (k in coln_gen0s) {
  icol <- k == names(dat0_sorted)
  if (!any(icol)) {
    dat0_sorted[[k]] <- template
  }
}


#--- Study type: survey, (manipulative) experiment, simulation
temp <- rep("Survey", dim(dat0_sorted)[1])
dat0_sorted$Study_type <- tag_comment(coms_unip, temp, c("fragmentation experiment"), "Experiment")
dat0_sorted$Study_type <- tag_comment(coms_unip, temp, c("simulation experiment"), "Simulation")


#------------------------------------
#--- SAVE CLEANED DATA
temp <- dat0$ID_record_orig
dat0 <- dat0_sorted[order(dat0_sorted$ID_record_orig), ]
stopifnot(temp == dat0$ID_record_orig)

saveRDS(dat0, file = foutw)
