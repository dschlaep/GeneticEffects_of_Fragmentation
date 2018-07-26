#--- SETTINGS
yin <- "E"
vin <- "var"

fragment_sizes <- c("Fragment_area_ha-published", "Pop_size_n-published",
    "Fragment_area_rank-published", "Pop_size_rank-published")
cor_methods <- c("pearson", "kendall", "spearman")
cor_transforms <- c("none", "unbiased", "ztransform")
weight_methods <- c("winvar", "whierarch")
withNonStandardized <- c("TRUE", "FALSE")
withNonStandardizedL <- c(TRUE, FALSE)
withControls <- c("TRUE", "FALSE")
withControlsL <- c(TRUE, FALSE)
effect_groups <- c("all", "Forests", "Non-Forests")

marker_classes <- c("codominant", "dominant", "proteins", "uniparental")

if (!exists("tag_outdate")) tag_outdate <- "20180709" # "20180112"
if (!exists("tag_usedate")) tag_usedate <- "20180705" #"20180112" #"20171107"


# Paths
dir_prj <- "Prj04_GenetFragment_MetaAnalysis"
dir_in <- file.path(dir_prj, "2_Data", "4_DataExtracted")
dir_out1 <- file.path(dir_prj, "2_Data", "5_DataCleaned")
dir_out2 <- file.path(dir_prj, paste0("4_Results_v", tag_outdate))
dir.create(dir_out2, recursive = TRUE, showWarnings = FALSE)
dir_big <- file.path("~", "Downloads", "BigData", paste0("4_Results_v", tag_usedate))
dir.create(dir_big, recursive = TRUE, showWarnings = FALSE)

dir_res_ <- file.path(dir_big, "__AnalysisOutput")
dir_res0 <- file.path(dir_out2, "0_Miscellaneous")
dir_res12 <- file.path(dir_out2, "1-2_GeneticDiversity&Inbreeding")
dir_res3 <- file.path(dir_out2, "3_TimeEffects")

temp <- sapply(c(dir_res_, dir_res0, dir_res12, dir_res3),
  function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE))

# File names
fout_esr <- file.path(dir_out1, paste0(tag_usedate, "_EffectSize_Correlation.rda"))
fout_mod <- file.path(dir_out1, paste0(tag_usedate, "_ModeratorPredictors.rds"))


#--- READ DATA
print(paste(Sys.time(), "reading input", shQuote(basename(fout_esr))))
load(fout_esr) # load: es_cor, var_options
deffects <- es_cor

print(paste(Sys.time(), "reading input", shQuote(basename(fout_mod))))
dmoderators <- readRDS(fout_mod)


#--- ANALYSIS DESIGN
responses_all <- c("Ar", "mA", "Np", "Fis", "FST", "H0", "He", "Sh", "He+Sh", "PLp")
# March 13, 2017: meeting decided on this subset
responses_H1 <- c("mA", "He+Sh", "PLp")
responses_H2 <- "Fis"
responses_Hall <- c(responses_H1, responses_H2)

cov_taxonomy <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species_resolved")

moderators_all <- names(dmoderators)[-(1:6)]
moderators_Hs <- c("Habitat3", "Matrix2", "Organism_group4")
#moderators_Hs_highres <- c("Habitat2", "Matrix2", "Organism_group5")
moderators_Hs_highres <- c("Habitat2")
#moderators_Htime <- c("SpeciesAgeCat2_ord_yrs", "EffectAgeCat2_ord_yrs",
#  "Effect_NGenSinceFrag_ord", "Effect_NlifeSinceFrag_ord")
moderators_Htime <- c("EffectAgeCat2_ord_yrs", "Effect_NGenSinceFrag_ord")
moderators_Hmisc <- c("Study_type", "Marker_cat",
  "Effect_latitude_bin", "Effect_range_age_yr",
  "Fragment_area_ha_range", "Pop_size_n_range")
moderators_plant_traits1 <- c("Plants_Lifeform", "Plants_Sex", "Plants_Self_compatibility")
moderators_plant_traits2 <- c("Plants_Mating_system", "Plants_Pollination_vector",
  "Plants_Seed_dispersal_vector")
moderators_plant_traits <- c(moderators_plant_traits1, moderators_plant_traits2)
moderators_animal_traits <- c("Animals_Length_cat", "Animals_Mass_cat",
  "Animals_Dispersal_ability", "Animals_Trophic_group")
moderators_traits <- c(moderators_plant_traits, moderators_animal_traits)
moderators_H12 <- unique(c(moderators_Hs, moderators_Hs_highres, moderators_traits))

std_design <- list(
  s1_wcontr = FALSE,
  s2_wnonnorm = TRUE,
  s3_responses = responses_Hall,
  s4_fragsize = "Fragment_area_ha-published",
  s5_cormethod = "pearson",
  s6_cortransform = "ztransform",
  d7_weightmethod = "whierarch")


full_design <- list(
  s1_wcontr = withControlsL,
  s2_wnonnorm = withNonStandardizedL,
  s3_responses = responses_Hall,
  s4_fragsize = fragment_sizes,
  s5_cormethod = cor_methods,
  s6_cortransform = cor_transforms,
  d7_weightmethod = weight_methods)


design_arguments <- list(
  args_target = list(
    withControlsL. = std_design[["s1_wcontr"]],
    only_useadj_standardized = FALSE, withNonStandardizedL. = std_design[["s2_wnonnorm"]],
    fragment_sizes. = std_design[["s4_fragsize"]],
    cor_methods. = std_design[["s5_cormethod"]],
    cor_transforms. = std_design[["s6_cortransform"]],
    weight_methods. = std_design[["d7_weightmethod"]]
  ),

  args_full = list(
    withControlsL. = full_design[["s1_wcontr"]],
    only_useadj_standardized = FALSE, withNonStandardizedL. = std_design[["s2_wnonnorm"]],
    fragment_sizes. = full_design[["s4_fragsize"]],
    cor_methods. = full_design[["s5_cormethod"]],
    cor_transforms. = full_design[["s6_cortransform"]],
    weight_methods. = full_design[["d7_weightmethod"]]
  ))
