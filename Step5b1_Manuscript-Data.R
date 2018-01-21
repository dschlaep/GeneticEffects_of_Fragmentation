#--- SETTINGS
redo <- FALSE
do_ms_data <- TRUE
interactions <- c(TRUE, FALSE)


# Paths
dir_prj1 <- "Prj04_GenetFragment_MetaAnalysis"
dir_ana <- file.path(dir_prj1, "3_Analysis")


#--- Setup
source(file.path(dir_ana, "Step5-0_Manuscript-methods.R"))


#--- Put together dataset for analysis

if (do_ms_data) {
  #--- targets
  template_args <-
    list(
      only_wo_controls = TRUE, withControlsL. = FALSE,
      only_useadj_standardized = FALSE, withNonStandardizedL. = TRUE,
      fragment_sizes. = std_design[["s4_fragsize"]],
      cor_methods. = std_design[["s5_cormethod"]],
      cor_transforms. = std_design[["s6_cortransform"]],
      weight_methods. = std_design[["d7_weightmethod"]],
      dir_res = dir_res_, dir_out = dir_ms_out
    )

  template_args <- expand.grid(template_args, stringsAsFactors = FALSE)

  for (k in seq_len(NROW(template_args))) {
    do_ms_loops(fun = "msres_getdatatogether", args = template_args[k, ],
      do_targets = TRUE, do_interactions = interactions)
  }


  #--- all
  template_args <-
    # list(
    #   only_wo_controls = FALSE, withControlsL. = full_design[["s1_wcontr"]],
    #   only_useadj_standardized = FALSE, withNonStandardizedL. = full_design[["s2_wnonnorm"]],
    #   fragment_sizes. = full_design[["s4_fragsize"]],
    #   cor_methods. = full_design[["s5_cormethod"]],
    #   cor_transforms. = full_design[["s6_cortransform"]],
    #   weight_methods. = full_design[["d7_weightmethod"]],
    #   dir_res = dir_res_, dir_out = dir_ms_out
    # )
    list(
      only_wo_controls = FALSE, withControlsL. = FALSE,
      only_useadj_standardized = FALSE, withNonStandardizedL. = TRUE,
      fragment_sizes. = full_design[["s4_fragsize"]],
      cor_methods. = full_design[["s5_cormethod"]],
      cor_transforms. = std_design[["s6_cortransform"]],
      weight_methods. = std_design[["d7_weightmethod"]],
      dir_res = dir_res_, dir_out = dir_ms_out
    )

  template_args <- expand.grid(template_args, stringsAsFactors = FALSE)

  for (k in seq_len(NROW(template_args))) {
    do_ms_loops(fun = "msres_getdatatogether", args = template_args[k, ],
      do_targets = FALSE, do_interactions = TRUE)
  }

}
