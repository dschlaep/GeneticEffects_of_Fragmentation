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
  template_args <- c(
      list(only_wo_controls = TRUE),
      design_arguments[["args_target"]],
      list(dir_res = dir_res_, dir_out = dir_ms_out)
    )

  template_args <- expand.grid(template_args, stringsAsFactors = FALSE)

  for (k in seq_len(NROW(template_args))) {
    do_ms_loops(fun = "msres_getdatatogether", args = template_args[k, ],
      do_targets = TRUE, do_interactions = interactions)
  }


  #--- all
  template_args <- c(
      list(only_wo_controls = FALSE),
      design_arguments[["args_full"]],
      list(dir_res = dir_res_, dir_out = dir_ms_out)
    )

  template_args <- expand.grid(template_args, stringsAsFactors = FALSE)

  for (k in seq_len(NROW(template_args))) {
    do_ms_loops(fun = "msres_getdatatogether", args = template_args[k, ],
      do_targets = FALSE, do_interactions = TRUE)
  }

}
