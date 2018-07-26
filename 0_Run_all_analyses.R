#!/usr/bin/env Rscript

do_prepare_data <- FALSE # prepare data for meta-analysis
do_calculate_effects <- FALSE # calculate effect sizes
do_meta_models <- FALSE # fit meta-analysis models and check model assumptions
do_output_misc <- FALSE # create miscellaneous output figures
do_output_prepare <- FALSE # prepare data for output
do_outputs <- FALSE # create output figures and tables

tag_newdate <- format(Sys.Date(), "%Y%m%d")
tag_usedate <- "20180705" #"20180112" # "20171107"
tag_outdate <- "20180709" #"20180112"


print(paste(Sys.time(), "0_Run_all_analysis.R: 'tag_newdate' =", tag_newdate,
  "- 'tag_usedate' =", tag_usedate, "- 'tag_outdate' =", tag_outdate))

if (do_prepare_data) {
  print(paste(Sys.time(), "Step1a_Clean-data.R"))
  source("Step1a_Clean-data.R")

  print(paste(Sys.time(), "Step1b_Convertions.R"))
  source("Step1b_Convertions.R")

  print(paste(Sys.time(), "Step1c_Moderators.R"))
  source("Step1c_Moderators.R")
}

if (do_calculate_effects) {
  print(paste(Sys.time(), "Step2_Calculate-effect-sizes.R"))
  source("Step2_Calculate-effect-sizes.R")
}


if (do_meta_models) {
  print(paste(Sys.time(), "Step3a_Analyse_effects.R"))
  source("Step3a_Analyse_effects.R")

  print(paste(Sys.time(), "Step3b_Analysis-checks.R"))
  source("Step3b_Analysis-checks.R")
}


if (do_output_misc) {
  print(paste(Sys.time(), "Step5a1_Manuscript_PRISM-diagram.R"))
  # Appendix S1: FIG. S1. PRISM diagram
  #   4_Results/0_Miscellaneous/Figs_PRISM/Fig_PRISM-plot.pptx and
  #   4_Results_v20180709/0_Miscellaneous/Figs_PRISM/Fig_PRISM-plot_v20180713.pdf
  source("Step5a1_Manuscript_PRISM-diagram.R")

  print(paste(Sys.time(), "Step5a2_Manuscript_Figure_Map.R"))
  # FIG. 1. Geographic distribution of study locations included in our meta-analysis
  #   4_Results_v20180709/0_Miscellaneous/Figs_Maps/Fig_MainMap_Fragment_area_ha-published_withoutControls_withNonStandardized_CORpearson-ztransform_whierarch_map.pdf
  # Appendix S5: FIG. S1. Geographic distribution of study locations separated by response variables and habitat types
  #   4_Results_v20180709/0_Miscellaneous/Figs_Maps/
  # Appendix S5: TABLE S1. Number of unique study locations by continent, response variables, and habitat types
  #   4_Results_v20180709/0_Miscellaneous/Figs_Maps/Table_MainMap_Fragment_area_ha-published_withoutControls_withNonStandardized_CORpearson-ztransform_whierarch_byHabitat3_map_perResponse.csv
  source("Step5a2_Manuscript_Figure_Map.R")

  print(paste(Sys.time(), "Step5a3_Manuscript_Figure_Checks.R"))
  # Appendix S4: TABLE S5. Tests for publication bias
  #   4_Results_v20180709/0_Miscellaneous/Tables/Table_PublicationBias_Fragment_area_ha-published_withoutControls_withNonStandardized_CORpearson-ztransform_whierarch_by_all.csv
  # Appendix S4: FIG. S1-S4. Diagnostic plots
  #   4_Results_v20180709/0_Miscellaneous/Figs_Checks/Fig_Diagnostics_Response_*_Fragment_area_ha-published_withoutControls_withNonStandardized_CORpearson-ztransform_whierarch_by_all.pdf
  source("Step5a3_Manuscript_Figure_Checks.R")
}

if (do_output_prepare) {
  print(paste(Sys.time(), "Step5b1_Manuscript-Data.R"))
  source("Step5b1_Manuscript-Data.R")

  print(paste(Sys.time(), "Step5b2_Manuscript_Sensitivity.R"))
  source("Step5b2_Manuscript_Sensitivity.R")
}

if (do_outputs) {
  print(paste(Sys.time(), "Step5c_Manuscript_Figure_MainForest.R"))
  # FIG. 2. Effect sizes of A, He, and PLP
  #   4_Results_v20180709/1-2_GeneticDiversity&Inbreeding/Figs_MainForest/Fig_MainForest_1-2_Diversity&Inbreeding-target-resp-gendiv_Fragment_area_ha-published_withoutControls_withNonStandardized_modXhabitat_CORpearson-ztransform_whierarch.pdf
  # Appendix S5: FIG. S2. Effect sizes of A, He, and PLP as a function of main effects
  #   4_Results_v20180709/1-2_GeneticDiversity&Inbreeding/Figs_MainForest/Fig_MainForest_1-2_Diversity&Inbreeding-target-resp-gendiv_Fragment_area_ha-published_withoutControls_withNonStandardized_CORpearson-ztransform_whierarch.pdf
  # Appendix S5: FIG. S3. Effect sizes of A, He, and PLP as a function of highly resolved habitat types.
  #   4_Results_v20180709/1-2_GeneticDiversity&Inbreeding/Figs_MainForest/Fig_MainForest_1-2_Diversity&Inbreeding-highres-target-resp-gendiv_Fragment_area_ha-published_withoutControls_withNonStandardized_CORpearson-ztransform_whierarch.pdf
  # Appendix S5: FIG. S4. Effect sizes of A, He, and PLP as a function of plant traits (part 1).
  #   4_Results_v20180709/1-2_GeneticDiversity&Inbreeding/Figs_MainForest/Fig_MainForest_1-2_Diversity&Inbreeding-PlantTraits1-target-resp-gendiv_Fragment_area_ha-published_withoutControls_withNonStandardized_modXhabitat_CORpearson-ztransform_whierarch.pdf
  # Appendix S5: FIG. S5. Effect sizes of A, He, and PLP as a function of plant traits (part 2).
  #   4_Results_v20180709/1-2_GeneticDiversity&Inbreeding/Figs_MainForest/Fig_MainForest_1-2_Diversity&Inbreeding-PlantTraits2-target-resp-gendiv_Fragment_area_ha-published_withoutControls_withNonStandardized_modXhabitat_CORpearson-ztransform_whierarch.pdf
  # Appendix S5: FIG. S6. Effect sizes of A, He, and PLP as a function of animal traits.
  #   4_Results_v20180709/1-2_GeneticDiversity&Inbreeding/Figs_MainForest/Fig_MainForest_1-2_Diversity&Inbreeding-AnimalTraits-target-resp-gendiv_Fragment_area_ha-published_withoutControls_withNonStandardized_modXhabitat_CORpearson-ztransform_whierarch.pdf
  # Appendix S5: FIG. S7. Effect sizes of Fis as a function of habitat and interactions
  #   4_Results_v20180709/1-2_GeneticDiversity&Inbreeding/Figs_MainForest/Fig_MainForest_1-2_Diversity&Inbreeding-target-resp-inbred_Fragment_area_ha-published_withoutControls_withNonStandardized_modXhabitat_CORpearson-ztransform_whierarch.pdf
  # Appendix S5: FIG. S8. Effect sizes of Fis as a function of highly resolved habitat types
  #   4_Results_v20180709/1-2_GeneticDiversity&Inbreeding/Figs_MainForest/Fig_MainForest_1-2_Diversity&Inbreeding-highres-target-resp-inbred_Fragment_area_ha-published_withoutControls_withNonStandardized_modXhabitat_CORpearson-ztransform_whierarch.pdf
  # Appendix S5: FIG. S9. Effect sizes of Fis as a function of plant traits (part 1).
  #   4_Results_v20180709/1-2_GeneticDiversity&Inbreeding/Figs_MainForest/Fig_MainForest_1-2_Diversity&Inbreeding-PlantTraits1-target-resp-inbred_Fragment_area_ha-published_withoutControls_withNonStandardized_modXhabitat_CORpearson-ztransform_whierarch.pdf
  # Appendix S5:
  #   4_Results_v20180709/1-2_GeneticDiversity&Inbreeding/Figs_MainForest/
  # Appendix S5:
  #   4_Results_v20180709/1-2_GeneticDiversity&Inbreeding/Figs_MainForest/
  # Appendix S5:
  #   4_Results_v20180709/1-2_GeneticDiversity&Inbreeding/Figs_MainForest/
  source("Step5c_Manuscript_Figure_MainForest.R")

  print(paste(Sys.time(), "Step5e_Manuscript_Table_Heterogeneity.R"))
  # Appendix S4: TABLES S1-S4. Tests of residual heterogeneity (QE) and omnibus F-tests of moderators (QM)
  #   4_Results_v20180709/1-2_GeneticDiversity&Inbreeding/Tables/Table_HeterogeneityTests_1-2_Diversity&Inbreeding-*.csv
  #   4_Results_v20180709/3_TimeEffects/Tables/Table_HeterogeneityTests_3_TimeEffects-*.csv
  source("Step5e_Manuscript_Table_Heterogeneity.R")

  print(paste(Sys.time(), "Step5g_Manuscript_SupplementaryTables.R"))
  # Appendix S2: TABLE S1. Assignments of habitat categories to habitat information
  #   4_Results_v20180709/0_Miscellaneous/Tables/Table_Habitat_Input_vs_Classification.csv
  # Appendix S3: FIG. S1. Animal phylogeny
  #   4_Results_v20180709/0_Miscellaneous/Figs_Phylogeny/Fig_Phylogeny_animals_woBL.pdf
  # Appendix S3: FIG. S2. Plant phylogeny
  #   4_Results_v20180709/0_Miscellaneous/Figs_Phylogeny/Fig_Phylogeny_plants_woBL.pdf
  source("Step5g_Manuscript_SupplementaryTables.R")

  print(paste(Sys.time(), "Step5h_Manuscript_Table_KBloomberg.R"))
  # Appendix S3: TABLE S1. Strength of the phylogenetic signal (Blomberg's K)
  #   4_Results_v20180709/1-2_GeneticDiversity&Inbreeding/Tables/Table_BlombergsK_1-2_*.csv
  #   4_Results_v20180709/3_TimeEffects/Tables/Table_BlombergsK_3_*.csv
  source("Step5h_Manuscript_Table_KBloomberg.R")
}

print(paste(Sys.time(), "0_Run_all_analysis.R: completed."))
