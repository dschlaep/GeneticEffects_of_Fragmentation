
#--- SETTINGS
redo <- FALSE
do_ms <- TRUE
do_targets <- TRUE

if (!exists("tag_usedate")) tag_usedate <- "20171107"


# Paths
dir_prj1 <- "Prj04_GenetFragment_MetaAnalysis"
dir_ana <- file.path(dir_prj, "3_Analysis")

#--- Setup
source(file.path(dir_ana, "Step5-0_Manuscript-methods.R"))

foutw <- file.path(dir_prj, "2_Data", "5_DataCleaned",
  paste0(tag_usedate, "_ExtractedDataCleaned3_wide.rds"))
dat0 <- readRDS(foutw)


#--- Number of species and effects
# Number of effects
dim(deffects)[1]
# 155

# number of animal and plant species
temp <- unique(dmoderators[, c("Organism_group4", "Species_resolved")])
ttemp <- table(temp[, "Organism_group4"])
# animals nonwoody-plants    woody-plants
# 52              42              43


#--- Table Habitat classification
utemp <- unique(dat0[, c("Habitat", "Habitat2", "Habitat3")])
utemp <- utemp[do.call(order, list(utemp[, "Habitat3"], utemp[, "Habitat2"], utemp[, "Habitat"])), ]

write.csv(utemp, file = file.path(dir_res0, "Tables", "Table_Habitat_Input_vs_Classification.csv"))
