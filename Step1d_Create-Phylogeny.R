
# Paths
dir_prj <- "Prj04_GenetFragment_MetaAnalysis"
dir_in <- file.path(dir_prj, "2_Data", "4_DataExtracted")
dir_phylo <- file.path(dir_in, "Phylogeny")
dir_SPhyloMaker <- file.path(dir_phylo, "S.PhyloMaker-master-f3948f2")
dir_ana <- file.path(dir_prj, "3_Analysis")

#--- Setup
library("ape")
library("phytools")

source(file.path(dir_ana, "Step3_0-Analysis-settings.R"))
source(file.path(dir_ana, "Step3_2_Analysis-Phylogeny.R"))


# Load 'S.PhyloMaker' https://github.com/jinyizju/S.PhyloMaker
source(file.path(dir_SPhyloMaker, "R_codes for S.PhyloMaker.R"))
PhytoPhylo_MegaTree <- read.tree(file.path(dir_SPhyloMaker, "PhytoPhylo.tre"))
PhytoPhylo_Nodes <- read.delim(file.path(dir_SPhyloMaker, "nodes.tab"), header = TRUE)

# example <- read.delim(file.path(dir_SPhyloMaker, "example.splist.tab"), header = TRUE)       # read in the example species list.
dspecies <- dmoderators[, c("Species_resolved", "Genus", "Family")]
colnames(dspecies) <- c("species", "genus", "family")
dspecies_plants <- dspecies[dmoderators[, "Kingdom"] %in% "Plantae", ]
dspecies_plants <- unique(dspecies_plants)

# run the function S.PhyloMaker
res_PhytoPhylo <- S.PhyloMaker(spList = dspecies_plants, tree = PhytoPhylo_MegaTree,
  nodes = PhytoPhylo_Nodes, scenarios = "S3")

par_prev <- par(mfrow = c(1, 2), mar = c(0, 0, 1, 0), cex = 1.1)
  plot(dphyloplants, main = "dphyloplants")
  plot(res_PhytoPhylo$Scenario.3, main = "Scenarion Three")
par(par_prev)


