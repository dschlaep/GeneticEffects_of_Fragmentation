#--- SETTINGS
redo <- FALSE
do_ms <- TRUE
do_targets <- TRUE
interactions <- TRUE


# Paths
dir_prj1 <- "Prj04_GenetFragment_MetaAnalysis"
dir_ana <- file.path(dir_prj1, "3_Analysis")


#--- Setup
source(file.path(dir_ana, "Step5-0_Manuscript-methods.R"))


#--- Main figures
msres_table_Heterogeneity <- function(responses, moderators, interaction_wHabitat3 = FALSE,
  fragment_sizes. = fragment_sizes, withControlsL. = withControlsL, only_wo_controls = TRUE,
  withNonStandardizedL. = withNonStandardizedL, only_useadj_standardized = TRUE,
  cor_methods. = cor_methods, cor_transforms. = cor_transforms,
  weight_methods. = weight_methods, dir_res, dir_out, ftag = "", ...) {

  #--- Load analysis outputs
  msdat <- msres_getdatatogether(responses, moderators, interaction_wHabitat3,
    fragment_sizes., withControlsL., only_wo_controls, withNonStandardizedL.,
    only_useadj_standardized, cor_methods., cor_transforms., weight_methods., dir_res,
    dir_out = dir_ms_out, ftag)

  stopifnot(NROW(msdat[["res"]]) == length(moderators),
    NCOL(msdat[["res"]]) == length(responses))

  #--- Gather data for table
  estimates1 <- c("QE", "QEp", "QM", "QMp")
  estimates2 <- c("Species_N", "QEdf", "QMdf1", "QMdf2", "I2_overall", "heterogeneity_total",
    "Htotal_explained")
  estimatesC1 <- c("beta", "ci.lb", "ci.ub", "pval")
  estimatesC2 <- c("Species_N", "tval")

  fit_types <- c("all_un", "all_mv", "aphylo", "pphylo")
  fit_types_legend <- c("Random-Effect", "Multi-Level", "Animal-Phylo", "Plant-Phylo")

  if (all(nchar(msdat[["moderatorsX"]]) == 3)) {
    mods_legend <- paste(msdat[["moderators"]], "x Habitat")
    modcats_legend <- msdat[["modcats"]]
    for (k in seq_along(msdat[["moderatorsX"]])) {
      modcats_legend <- gsub(paste0(msdat[["moderatorsX"]][k], k),
        paste0(msdat[["moderators"]][k], k), modcats_legend)
    }
  } else {
    mods_legend <- msdat[["moderators"]]
    modcats_legend <- msdat[["modcats"]]
  }

  dat_het <- list(
    omnibus = array(NA, dim = c(length(responses), length(mods_legend), length(fit_types),
      length(estimates1) + length(estimates2)), dimnames = list(responses, mods_legend,
      fit_types_legend, c(estimates1, estimates2))),

    modlevs = array(NA, dim = c(length(responses), length(modcats_legend),
      length(fit_types), length(estimatesC1) + length(estimatesC2)), dimnames =
      list(responses, modcats_legend, fit_types_legend, c(estimatesC1, estimatesC2)))
  )

  for (ir in seq_along(responses)) {
    for (ig in seq_along(moderators)) {
      x <- msdat[["res"]][moderators[ig], responses[ir]][[1]]
      fits <- list(
        all_un = x[["all"]][["fit"]][["uni"]],
        all_mv = x[["all"]][["fit"]][["mv"]],
        aphylo = x[["aphylo"]][["fit"]][["woBL"]],
        pphylo = x[["pphylo"]][["fit"]][["woBL"]])

      x2 <- msdat[["data"]][moderators[ig], responses[ir]][[1]]
      dats <- list(
        all_un = x2[["all"]],
        all_mv = x2[["all"]],
        aphylo = x2[["aphylo"]],
        pphylo = x2[["pphylo"]])

      has_fits <- sapply(fits, function(f) f[["mod"]][["has_fit"]])

      if (any(has_fits)) {
        for (ft in seq_along(fit_types)) {
          fti <- fit_types[ft]
          ftl <- fit_types_legend[ft]

          if (has_fits[fti]) {
            xm <- fits[[fti]][["mod"]][["fit"]]
            xd <- dats[[fti]]

            #- Omnibus-type output
            dat_het[["omnibus"]][ir, ig, ftl, estimates1] <- as.numeric(xm[estimates1])

            # Should be identical: identical(length(unique(xd[["d"]][, "Species"])), as.integer(xm[["dfs"]]))
            dat_het[["omnibus"]][ir, ig, ftl, "Species_N"] <-
              length(unique(xd[["d"]][, "Species"]))

            dat_het[["omnibus"]][ir, ig, ftl, c("I2_overall", "heterogeneity_total")] <-
              if (inherits(xm, "rma.uni")) {
                as.numeric(xm[c("I2", "tau2")])
              } else {
                c(sum(as.numeric(xm[["I2_overall"]])), as.numeric(xm[["heterogeneity_total"]]))
              }

            dat_het[["omnibus"]][ir, ig, ftl, "Htotal_explained"] <-
              as.numeric(fits[[fti]][["Htotal_explained"]])

            dat_het[["omnibus"]][ir, ig, ftl, "QEdf"] <- as.numeric(xm[["k"]] - xm[["p"]])

            dat_het[["omnibus"]][ir, ig, ftl, "QMdf1"] <- as.numeric(xm[["m"]])
            if (is.element(xm[["test"]], c("knha", "adhoc", "t"))) {
              dat_het[["omnibus"]][ir, ig, ftl, "QMdf2"] <- as.numeric(xm[["dfs"]])
            }


            #- Output per moderator level
            temp <- names(xm[["coef.na"]])[!xm[["coef.na"]]]
            coefnames <- if (interaction_wHabitat3) {
                paste0(substr(temp, 2, 4), ig, substr(temp, 5, nchar(temp)))
              } else {
                temp <- substr(temp, 2, nchar(temp))
                for (i in seq_along(msdat[["moderators"]])) {
                  temp <- sub(msdat[["moderators"]][i],
                    paste0(msdat[["moderators"]][i], ig), temp)
                }
                temp
              }

            ids <- match(msdat[["modcats"]], coefnames, nomatch = 0)

            temp <- table(xd[["d"]][, "Species"], xd[["d"]][, xd[["mods"]][["names"]]])
            temp <- apply(temp, 2, function(x) sum(x > 0))
            names(temp) <- paste0(xd[["mods"]][["names"]], names(temp))
            dat_het[["modlevs"]][ir, ids > 0, ftl, "Species_N"] <-
              temp[colnames(xm[["X"]])][ids]

            if (is.element(xm[["test"]], c("knha", "adhoc", "t"))) {
              dat_het[["modlevs"]][ir, ids > 0, ftl, "tval"] <- as.numeric(xm[["zval"]])[ids]
            } else {
              print(paste("Code to lookup z-value for test =", xm[["test"]], "not implemented."))
            }

            for (k in estimatesC1) {
              dat_het[["modlevs"]][ir, ids > 0, ftl, k] <- as.numeric(xm[[k]])[ids]
            }
          }
        }
      }
    }
  }



  #--- Create manuscript table
  temp <- reshape2::melt(dat_het[["omnibus"]])
  dat_out <- reshape2::dcast(temp, Var1 + Var2 + Var3 ~ Var4)
  names(dat_out)[1:3] <- c("Response", "Moderator", "Model")
  dat_out[, "Moderator"] <- clean_labels(dat_out[, "Moderator"], remove_digits = TRUE,
    break_lines = FALSE)

  #--- Save table to file
  dir <- file.path(dir_out, "Tables")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(dat_out, file = file.path(dir, paste0("Table_HeterogeneityTests_", ftag, "_",
    msdat[["tag_fix"]], ".csv")))


  temp <- reshape2::melt(dat_het[["modlevs"]])
  dat_out <- reshape2::dcast(temp, Var1 + Var2 + Var3 ~ Var4)
  names(dat_out)[1:3] <- c("Response", "Moderator", "Model")

  temp <- paste0(names(msdat[["cats"]]), seq_along(msdat[["cats"]]))
  temp <- data.frame(Moderator = rep(temp, msdat[["n_cats_max"]]),
    Level = unlist(msdat[["cats"]]), stringsAsFactors = FALSE)
  stopifnot(as.character(dat_out[, "Moderator"]) ==
    rep(unname(apply(temp, 1, paste, collapse = "")), each = length(fit_types)))
  dat_out <- data.frame(Response = dat_out[, "Response"],
    Moderator = rep(temp[, "Moderator"], each = length(fit_types)),
    Level = rep(temp[, "Level"], each = length(fit_types)),
    dat_out[, -which(names(dat_out) %in% c("Response", "Moderator"))],
    stringsAsFactors = FALSE)

  dat_out[, "Moderator"] <- clean_labels(dat_out[, "Moderator"], remove_digits = TRUE,
    break_lines = FALSE)
  dat_out[, "Level"] <- clean_labels(dat_out[, "Level"], remove_digits = FALSE,
    break_lines = FALSE)


  #--- Save table to file
  dir <- file.path(dir_out, "Tables")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  write.csv(dat_out, file = file.path(dir, paste0("Table_ModeratorLevelTests_", ftag, "_",
    msdat[["tag_fix"]], ".csv")))

  invisible(TRUE)
}



if (do_ms) {
  template_args <- if (do_targets) {
      list(
        only_wo_controls = TRUE, withControlsL. = std_design[["s1_wcontr"]],
        only_useadj_standardized = FALSE, withNonStandardizedL. = std_design[["s2_wnonnorm"]],
        fragment_sizes. = std_design[["s4_fragsize"]],
        cor_methods. = std_design[["s5_cormethod"]],
        cor_transforms. = std_design[["s6_cortransform"]],
        weight_methods. = std_design[["d7_weightmethod"]],
        dir_res = dir_res_, panels_vertical = FALSE
      )

    } else {
      # list(
      #   only_wo_controls = FALSE, withControlsL. = full_design[["s1_wcontr"]],
      #   only_useadj_standardized = FALSE, withNonStandardizedL. = full_design[["s2_wnonnorm"]],
      #   fragment_sizes. = full_design[["s4_fragsize"]],
      #   cor_methods. = full_design[["s5_cormethod"]],
      #   cor_transforms. = full_design[["s6_cortransform"]],
      #   weight_methods. = full_design[["d7_weightmethod"]],
      #   dir_res = dir_res_
      # )
      list(
        only_wo_controls = FALSE, withControlsL. = FALSE,
        only_useadj_standardized = FALSE, withNonStandardizedL. = TRUE,
        fragment_sizes. = full_design[["s4_fragsize"]],
        cor_methods. = full_design[["s5_cormethod"]],
        cor_transforms. = std_design[["s6_cortransform"]],
        weight_methods. = std_design[["d7_weightmethod"]],
        dir_res = dir_res_, panels_vertical = FALSE
      )
    }

  template_args <- expand.grid(template_args, stringsAsFactors = FALSE)

  args_along <- list(
      dir_out = list(dir_res0, dir_res12, dir_res12, dir_res12, dir_res3)
    )

  for (k in seq_len(NROW(template_args))) {
    do_ms_loops(fun = "msres_table_Heterogeneity", args = template_args[k, ],
      args_along = args_along, do_targets = do_targets, do_interactions = interactions)
  }
}

