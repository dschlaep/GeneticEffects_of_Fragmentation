
# Paths
dir_prj <- "Prj04_GenetFragment_MetaAnalysis"
dir_in <- file.path(dir_prj, "2_Data", "4_DataExtracted")
dir_phylo <- file.path(dir_in, "Phylogeny")

# File names
f_phylo_plants <- file.path(dir_phylo, "20171004_MetaTotNrPlantsNew.nwk")
f_phylo_animals <- file.path(dir_phylo, "20171010_NewickTreeAnimalsMeta.nwk")


#--- READ DATA
print(paste(Sys.time(), "reading input", shQuote(basename(f_phylo_plants))))
dphyloplants <- ape::read.tree(f_phylo_plants)

print(paste(Sys.time(), "reading input", shQuote(basename(f_phylo_animals))))
dphyloanimals <- ape::read.tree(f_phylo_animals)

#--- Validate phylogenetic trees
validate_tree <- function(tree) {
  check1 <- inherits(try(ape::checkValidPhylo(tree)), "try-error")
  check2 <- !ape::is.ultrametric(tree)
  check3 <- anyDuplicated(tree[["tip.label"]]) > 0

  if (check1 || check2 || check3) {
    stop(shQuote(match.call()[[2]]), " phylogenetic tree is not ultrametric,",
      " it may not fit a BM model of evolution.")
  }

  TRUE
}

validate_tree(dphyloplants)
validate_tree(dphyloanimals)

