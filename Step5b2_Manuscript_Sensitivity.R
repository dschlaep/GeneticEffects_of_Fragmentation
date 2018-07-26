#--- SETTINGS
redo <- FALSE
do_ms <- TRUE
interactions <- c(TRUE, FALSE)


# Paths
dir_prj1 <- "Prj04_GenetFragment_MetaAnalysis"
dir_ana <- file.path(dir_prj1, "3_Analysis")


#--- Setup
source(file.path(dir_ana, "Step5-0_Manuscript-methods.R"))


#--- Sensitivity of results
msres_sensitivity <- function(responses, moderators, interaction_wHabitat3 = FALSE,
  args_target, args_full, dir_res, dir_out, ftag = "", panels_vertical = TRUE,
  ...) {

  temp_xHabitat3 <- if (interaction_wHabitat3) "modXhabitat" else ""
  tag_fix <- if (interaction_wHabitat3) {
      paste0(ftag, "_", temp_xHabitat3)
    } else {
      ftag
    }

  design <- expand.grid(args_full, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  is_target <- rep(TRUE, NROW(design))
  for (k in seq_along(args_target)) {
    is_target <- is_target & design[, names(args_target)[[k]]] == args_target[[k]]
  }
  stopifnot(sum(is_target) == 1)
  design[, "is_target"] <- is_target

  #--- Load analysis outputs: only 'responses' and 'moderators' allowed to vary
  msdat <- vector(mode = "list", length = NROW(design))
  for (k in seq_len(NROW(design))) {
    msdat[[k]] <- msres_getdatatogether(responses, moderators,
      interaction_wHabitat3,
      fragment_sizes. = design[k, "fragment_sizes."],
      withControlsL. = design[k, "withControlsL."],
      only_wo_controls = FALSE, # i.e., use value(s) of `withControlsL.`
      withNonStandardizedL. = design[k, "withNonStandardizedL."],
      only_useadj_standardized = design[k, "only_useadj_standardized"],
      cor_methods. = design[k, "cor_methods."],
      cor_transforms. = design[k, "cor_transforms."],
      weight_methods. = design[k, "weight_methods."],
      dir_res, dir_out = dir_ms_out, ftag)
  }

  modcats <- lapply(msdat, function(x) x[["modcats"]])
  stopifnot(any(lengths(modcats) > 0))
  modcats <- sort(unique(unlist(modcats)))

  #--- data
  dir_sens <- file.path(dir_res, "Datasets_Sensitivity")
  dir.create(dir_sens, showWarnings = FALSE)
  ftemp <- file.path(dir_sens, paste0("Data_Sensitivity_", tag_fix, ".rds"))

  if (redo || !file.exists(ftemp)) {
    id_target <- sapply(names(args_target), function(name)
      which(args_full[[name]] == args_target[[name]]))
    stopifnot(lengths(args_target) == 1, length(interaction_wHabitat3) == 1,
      !any(is.null(id_target)))

    outcomes <- c("CI_lt0", "CI_with0", "CI_gt0")
    fit_types <- c("all_un", "all_mv", "aphylo", "pphylo")
    fit_types_legend <- c("Random-Effect", "Multi-Level", "Animal-Phylo", "Plant-Phylo")

    #--- Get model results
    out_res <- array(NA,
      dim = c(length(responses), length(modcats), NROW(design), length(fit_types), length(outcomes)),
      dimnames = list(responses, modcats, NULL, fit_types_legend, outcomes))

    for (id in seq_len(NROW(design))) {
      for (ir in seq_along(responses)) {
        for (ig in seq_along(moderators)) {
          x <- msdat[[id]][["res"]][moderators[ig], responses[ir]][[1]]

          fits <- list(
            all_un = x[["all"]][["fit"]][["uni"]][["mod"]],
            all_mv = x[["all"]][["fit"]][["mv"]][["mod"]],
            aphylo = x[["aphylo"]][["fit"]][["woBL"]][["mod"]],
            pphylo = x[["pphylo"]][["fit"]][["woBL"]][["mod"]])

          has_fits <- sapply(fits, function(f) f[["has_fit"]])

          if (any(has_fits)) {
            coefnames <- list()
            for (ft in fit_types) if (has_fits[ft]) {
              temp <- names(fits[[ft]][["fit"]][["coef.na"]])[!fits[[ft]][["fit"]][["coef.na"]]]
              coefnames[[ft]] <- if (interaction_wHabitat3) {
                  paste0(substr(temp, 2, 4), ig, substr(temp, 5, nchar(temp)))
                } else {
                  temp <- substr(temp, 2, nchar(temp))
                  for (i in seq_along(msdat[[id]][["moderators"]])) {
                    temp <- sub(msdat[[id]][["moderators"]][i],
                      paste0(msdat[[id]][["moderators"]][i], ig), temp)
                  }
                  temp
                }
            }

            for (ift in seq_along(fit_types)) if (has_fits[ift]) {
              ft <- fit_types[ift]
              ids <- match(modcats, coefnames[[ft]], nomatch = 0)
              ids0 <- ids > 0

              # Estimates per level
              temp <- as.matrix(data.frame(fits[[ft]][["fit"]][c("ci.lb", "ci.ub")]))
              stopifnot(sum(ids0) == NROW(temp))

              out_res[ir, ids0, id, ift, "CI_lt0"] <- apply(temp, 1,
                function(x) all(x < 0))
              out_res[ir, ids0, id, ift, "CI_with0"] <- apply(temp, 1,
                function(x) x[1] <= 0 && x[2] >= 0)
              out_res[ir, ids0, id, ift, "CI_gt0"] <- apply(temp, 1,
                function(x) all(x > 0))
            }
          }
        }
      }
    }

    #--- Target outcome:
    # wrap subset in `array` to make sure no dimension is beeing dropped
    # (e.g., if responses = "Fis")
    target <- array(out_res[, , design[, "is_target"], , ], dim = dim(out_res)[-3],
      dimnames = dimnames(out_res)[-3])

    #--- Calculate sensivity:
    res <- list(
      design = design,
      N = nrow(design),
      # Sensitivity: sum across 'design' dimension of 'out_res
      sensitivity = apply(out_res, c(1:2, 4:5), sum, na.rm = TRUE),
      target = target,
      out_res = out_res
    )

    saveRDS(res, file = ftemp)

  } else {
    res <- readRDS(ftemp)
  }


  #-- Convert 'sensitivity' array into table and write to disk file
  dir_temp <- file.path(dir_out, "Tables")
  dir.create(dir_temp, recursive = TRUE, showWarnings = FALSE)

  ftemp <- file.path(dir_temp, paste0("Table_Sensitivity_", tag_fix, ".csv"))

  if (redo || !file.exists(ftemp)) {
    temp <- reshape2::melt(res[["sensitivity"]])
    table_sens <- reshape2::dcast(temp, Var1 + Var2 + Var3 ~ Var4)
    colnames(table_sens)[1:3] <- c("Response", "Moderator", "Model")

    temp <- reshape2::melt(res[["target"]])
    table_target <- reshape2::dcast(temp, Var1 + Var2 + Var3 ~ Var4)
    colnames(table_target) <- c(colnames(table_sens)[1:3],
      paste0("target_", colnames(table_sens)[4:6]))

    stopifnot(table_sens[, 1:3] == table_target[, 1:3])

    out <- cbind(table_sens, N = res[["N"]], table_target[, -(1:3)])

    write.csv(out, file = ftemp, row.names = FALSE)
  }

  invisible(TRUE)
}



if (do_ms) {

  args <- list(
    args_target = design_arguments[["args_target"]],
    args_full = design_arguments[["args_full"]],

    dir_res = dir_res_, panels_vertical = FALSE
  )

  args_along <- list(
    dir_out = list(dir_res0, dir_res12, dir_res12, dir_res12, dir_res3)
  )

  # 144 (full_design[["s1_wcontr"]]), 72 (std_design[["s1_wcontr"]])
  print(paste("N =", prod(sapply(args[["args_full"]], length))))

  do_ms_loops(fun = "msres_sensitivity", args = args,
    args_along = args_along, do_targets = TRUE, do_interactions = interactions)
}

