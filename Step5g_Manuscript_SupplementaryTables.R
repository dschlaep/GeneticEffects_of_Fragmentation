
#--- SETTINGS
redo <- FALSE
do_ms <- TRUE
do_targets <- TRUE

if (!exists("tag_usedate")) tag_usedate <- "20180112"


# Paths
dir_prj1 <- "Prj04_GenetFragment_MetaAnalysis"
dir_ana <- file.path(dir_prj, "3_Analysis")

#--- Setup
source(file.path(dir_ana, "Step5-0_Manuscript-methods.R"))

foutw <- file.path(dir_prj, "2_Data", "5_DataCleaned",
  paste0(tag_usedate, "_ExtractedDataCleaned3_wide.rds"))
dat0 <- readRDS(foutw)


#--------------------------------------
#--- Number of species and effects
# Number of effects
dim(deffects)[1]
# 155

# number of animal and plant species
temp <- unique(dmoderators[, c("Organism_group4", "Species_resolved")])
ttemp <- table(temp[, "Organism_group4"])
print(ttemp)
# animals nonwoody-plants    woody-plants
# 52              42              43


#--------------------------------------
#--- Phylogeny figures

source(file.path(dir_ana, "Step3_2_Analysis-Phylogeny.R"))

dir_temp <- file.path(dir_res0, "Figs_Phylogeny")
dir.create(dir_temp, showWarnings = FALSE)

pdf(height = 5, width = 6, file = file.path(dir_temp, "Fig_Phylogeny_animals_woBL.pdf"))
ape::plot.phylo(dphyloanimals, use.edge.length = FALSE, node.depth = 2,
  no.margin = TRUE, cex = 0.5, label.offset = 0.25)
dev.off()

pdf(height = 8, width = 6, file = file.path(dir_temp, "Fig_Phylogeny_plants_woBL.pdf"))
ape::plot.phylo(dphyloplants, use.edge.length = FALSE, node.depth = 2,
  no.margin = TRUE, cex = 0.5, label.offset = 0.25)
dev.off()

pdf(height = 8, width = 6, file = file.path(dir_temp, "Fig_Phylogeny_plants_wBL.pdf"))
ape::plot.phylo(dphyloplants,
  no.margin = TRUE, cex = 0.5, label.offset = 5)
ape::add.scale.bar()
dev.off()

#--------------------------------------
#--- Table Habitat classification
utemp <- unique(dat0[, c("Habitat", "Habitat2", "Habitat3")])
utemp <- utemp[do.call(order, list(utemp[, "Habitat3"], utemp[, "Habitat2"], utemp[, "Habitat"])), ]

write.csv(utemp, file = file.path(dir_res0, "Tables", "Table_Habitat_Input_vs_Classification.csv"))
