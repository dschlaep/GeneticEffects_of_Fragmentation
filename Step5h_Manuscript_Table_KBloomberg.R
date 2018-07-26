#--- SETTINGS
redo <- FALSE
do_ms <- TRUE
do_targets <- TRUE


# Paths
dir_prj1 <- "Prj04_GenetFragment_MetaAnalysis"
dir_ana <- file.path(dir_prj1, "3_Analysis")


#--- Setup
source(file.path(dir_ana, "Step5-0_Manuscript-methods.R"))


#--- Main figures
msres_table_BloombergsK <- function(responses, moderators,
  fragment_sizes. = fragment_sizes, withControlsL. = withControlsL, only_wo_controls = TRUE,
  withNonStandardizedL. = withNonStandardizedL, only_useadj_standardized = TRUE,
  cor_methods. = cor_methods, cor_transforms. = cor_transforms,
  weight_methods. = weight_methods, dir_res, dir_out, ftag = "", ...) {

  #--- Load analysis outputs
  msdat <- msres_getdatatogether(responses, moderators, interaction_wHabitat3 = FALSE,
    fragment_sizes., withControlsL., only_wo_controls, withNonStandardizedL.,
    only_useadj_standardized, cor_methods., cor_transforms., weight_methods., dir_res,
    dir_out = dir_ms_out, ftag)

  stopifnot(NROW(msdat[["res"]]) == length(moderators),
    NCOL(msdat[["res"]]) == length(responses))

  #--- Gather data for table
  estimates <- c("K", "P")
  fit_types <- c("aphylo_wBL", "aphylo_woBL", "pphylo_wBL", "pphylo_woBL")
  fit_types_legend <- c("Animal-Phylo (BL)", "Animal-Phylo (w/o BL)",
    "Plant-Phylo (BL)", "Plant-Phylo (w/o BL)")

  dat_table <- list(
    est = array(NA, dim = c(length(responses), length(moderators), length(fit_types),
      length(estimates)), dimnames = list(responses, moderators, fit_types_legend, estimates))
  )

  for (ir in seq_along(responses)) {
    for (ig in seq_along(moderators)) {
      x <- msdat[["res"]][moderators[ig], responses[ir]][[1]]

      for (ift in seq_along(fit_types)) {
        fit <- switch(fit_types[ift],
          aphylo_wBL = NULL,
          aphylo_woBL = unlist(x[["aphylo"]][["add"]][["BlombergsK_woBL"]]),
          pphylo_wBL = unlist(x[["pphylo"]][["add"]][["BlombergsK_wBL"]]),
          pphylo_woBL = unlist(x[["pphylo"]][["add"]][["BlombergsK_woBL"]]),
          NA)

        if (all(estimates %in% names(fit))) {
          dat_table[["est"]][ir, ig, ift, ] <- fit[estimates]
        }
      }
    }
  }



  #--- Create manuscript table
  temp <- reshape2::melt(dat_table[["est"]])
  dat_out <- reshape2::dcast(temp, Var2 + Var3 ~ Var1 + Var4)
  names(dat_out)[1:2] <- c("Moderator", "Phylogeny")
  dat_out[, "Moderator"] <- clean_labels(as.character(dat_out[, "Moderator"]),
    remove_digits = TRUE, break_lines = FALSE)

  #--- Save table to file
  dir <- file.path(dir_out, "Tables")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(dat_out, file = file.path(dir, paste0("Table_BlombergsK_", ftag, "_",
    msdat[["tag_fix"]], ".csv")))

  invisible(TRUE)
}



if (do_ms) {
  template_args <- if (do_targets) {
      c(
        list(only_wo_controls = TRUE),
        design_arguments[["args_target"]],
        list(dir_res = dir_res_, panels_vertical = FALSE)
      )

    } else {
      c(
        list(only_wo_controls = FALSE),
        design_arguments[["args_full"]],
        list(dir_res = dir_res_, panels_vertical = FALSE)
      )
    }

  template_args <- expand.grid(template_args, stringsAsFactors = FALSE)

  args_along <- list(
      dir_out = list(dir_res0, dir_res12, dir_res12, dir_res12, dir_res3)
    )

  for (k in seq_len(NROW(template_args))) {
    do_ms_loops(fun = "msres_table_BloombergsK", args = template_args[k, ],
      args_along = args_along, do_targets = do_targets)
  }
}

