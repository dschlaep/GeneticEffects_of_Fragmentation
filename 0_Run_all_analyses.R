#!/usr/bin/env Rscript

do_prepare_data <- FALSE # prepare data for meta-analysis
do_calculate_effects <- FALSE # calculate effect sizes
do_meta_models <- FALSE # fit meta-analysis models and check model assumptions
do_output_misc <- FALSE # create miscellaneous output figures
do_output_prepare <- FALSE # prepare data for output
do_outputs <- FALSE # create output figures and tables

tag_newdate <- format(Sys.Date(), "%Y%m%d")
tag_usedate <- "20180112" # "20171107"


if (do_prepare_data) {
  source("Step1a_Clean-data.R")
  source("Step1b_Convertions.R")
  source("Step1c_Moderators.R")
}

if (do_calculate_effects) {
  source("Step2_Calculate-effect-sizes.R")
}


if (do_meta_models) {
  source("Step3a_Analyse_effects.R")
  source("Step3b_Analysis-checks.R")
}


if (do_output_misc) {
  source("Step5a1_Manuscript_PRISM-diagram.R")
  source("Step5a2_Manuscript_Figure_Map.R")
  source("Step5a3_Manuscript_Figure_Checks.R")
}

if (do_output_prepare) {
  source("Step5b1_Manuscript-Data.R")
  source("Step5b2_Manuscript_Sensitivity.R")
}

if (do_outputs) {
  source("Step5c_Manuscript_Figure_MainForest.R")
  source("Step5e_Manuscript_Table_Heterogeneity.R")
  source("Step5g_Manuscript_SupplementaryTables.R")
  source("Step5h_Manuscript_Table_KBloomberg.R")
}
