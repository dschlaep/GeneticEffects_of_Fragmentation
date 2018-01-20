# Overall checks for simple meta-analysis models without including moderators or accounting for dependencies

###TODO:
  # - phylogenetic relatedness
  #   - Egger's regression for phylo: https://stackoverflow.com/questions/34270651/trimfill-for-rma-mv-models-in-metafor-packlage

  #--- Multicollinearity among moderators


# Paths
dir_prj <- "Prj04_GenetFragment_MetaAnalysis"
dir_ana <- file.path(dir_prj, "3_Analysis")

#--- Setup
source(file.path(dir_ana, "Step3_0-Analysis-settings.R"))
source(file.path(dir_ana, "Step3_1_Analysis-methods.R"))


#--- SETTINGS
do_targets <- TRUE
do_checks <- TRUE



#------ Testing for publication bias
attempt_trimfill <- function(x_nomod, estimator = "R0") {
  if (x_nomod$k > 2) {

    tf <- try(metafor::trimfill(x_nomod, estimator = estimator), silent = TRUE)

    if (inherits(tf, "try-error")) {
      tf <- try(metafor::trimfill(x_nomod, estimator = estimator,
        control = list(maxiter = 1000, stepadj = 0.5)), silent = TRUE)
    }

    if (inherits(tf, "try-error")) {
      print(tf)
    }
  } else {
    tf <- try(stop("nrow(data) <= 2"), silent = TRUE)
  }

  tf
}


add_publicationbias <- function(res) {

  .local <- function() {
    if (res[["has_fit"]]) {
      #--- Three-parameter selection model test
      res[["3psm"]] <- threePSM.est(d = res[["fit"]][["yi"]], v = res[["fit"]][["vi"]],
        min.pvalues = 3)

      #--- Funnel plot tests
      res[["trim-fill"]] <- attempt_trimfill(res[["fit"]])

      # Mixed-effects 'Egger test'
      res[["eggers-test-mixed"]] <- metafor::regtest(res[["fit"]], model = "rma",
        predictor = "sei")

      # Begg & Mazumbar's 1994 rank correlation test for funnel plot asymmetry
      res[["begg-mazumbar"]] <- metafor::ranktest(res[["fit"]])
    }
  }

  suppressWarnings(temp <- try(.local()))

  res[["has"]][["publicationbias"]] <- !inherits(temp, "try-error")

  res
}


#------ Testing for temporal pattern in effect size
add_temporalpatterns <- function(res) {
  # Cumulative meta-analysis (in the order of publication year)
  # overall:
  if (res[["has_fit"]]) {
    yrs <- res[["data"]][, "Pub_Year"]
    oyrs <- order(yrs)
    res[["cumtime"]] <- list(
      fit = metafor::cumul(res[["fit"]], order = oyrs),
      yrs_sorted = yrs[oyrs])
  }

  res[["has"]][["temporalpatterns"]] <- TRUE

  res
}

#------ Model sensitivity 2
add_modelsensitivity2 <- function(res, yin = "E", vin = "var") {
  if (res[["has_fit"]]) {
    # Fail-safe numbers
    res[["Nfs"]] <- calc_failsafeN(res[["data"]], yin, vin)

    # Influence measures
    res[["infl"]] <- calc_influence(res)
  }

  res[["has"]][["modelsensitivity"]] <- TRUE

  res
}


run_metachecks <- function(responses, fragment_sizes. = fragment_sizes,
  withControlsL. = withControlsL, only_wo_controls = TRUE,
  withNonStandardizedL. = withNonStandardizedL, only_useadj_standardized = TRUE,
  cor_methods. = cor_methods, cor_transforms. = cor_transforms,
  weight_methods. = weight_methods, effect_groups. = effect_groups, deffects, dmoderators,
  dir_res, yin = "E", vin = "var") {

  temp_withControlsL <- if (only_wo_controls) FALSE else withControlsL.

  for (ir in responses) {
    dtag <- paste0("Response_", ir)
    dir_resout <- file.path(dir_res, "Results", dtag)
    dir.create(dir_resout, recursive = TRUE, showWarnings = FALSE)

    temp_withNonStandardizedL <- if (only_useadj_standardized || ir %in% c("Ar", "mA")) {
        # always use only standardized values for Ar and mA
        FALSE
      } else {
        withNonStandardizedL.
      }
    for (indep in fragment_sizes.) {
      for (wcontr in temp_withControlsL) {
        temp1 <- if (wcontr) "withControls" else "withoutControls"
        for (wnonnorm in temp_withNonStandardizedL) {
          temp2 <- if (wnonnorm) "withNonStandardized" else "onlyStandardized"
          for (method in cor_methods.) {
            temp_cor_transforms <- if (method == "pearson") {
                cor_transforms.
              } else {
                c("none", "ztransform")
              }
            for (tf in temp_cor_transforms) {
              for (win in weight_methods.) {
                for (ieg in effect_groups.) {

                  ftag <- paste(ir, "by", indep, temp1, temp2, paste0("COR", method, "-", tf),
                    win, "by", ieg, sep = "_")

                  print(paste(Sys.time(), "'run checks':", ftag))

                  ftemp <- file.path(dir_resout, paste0(ftag, "_checks.rds"))

                  res <- new.env()

                  # Data
                  temp <- data.frame(ID_effect = dmoderators[, "ID_effect"],
                    deffects[, 2 - as.integer(wcontr), 2 - as.integer(wnonnorm), ir, indep, method, tf, ],
                    Pub_Year = dmoderators[, "Pub_Year"], stringsAsFactors = FALSE)

                  if (ieg == "all") {
                    temp <- temp

                  } else if (ieg == "Forests") {
                    temp <- temp[dmoderators[, "Habitat3"] %in% ieg, , drop = FALSE]

                  } else if (ieg == "Non-Forests") {
                     temp <- temp[dmoderators[, "Habitat3"] %in% ieg, , drop = FALSE]

                  } else stop("Not implemented: ", ieg)

                  res[["data"]] <- temp[complete.cases(temp[, c(yin, vin, win)]), ]

                  if (nrow(res[["data"]]) > 1) {
                    wi <- if (!identical(win, "winvar")) {
                        res[["data"]][, win]
                      } else NULL

                    # Fit overall random-effects model without moderators
                    res[["fit"]] <- attempt_uni_rma(data = res[["data"]],
                      yi = res[["data"]][, yin], vi = res[["data"]][, vin], fmods = NULL,
                      wi = wi)
                    res[["has_fit"]] <- !inherits(res[["fit"]], "try-error")


                    # Add checks
                    res <- add_publicationbias(res)
                    res <- add_temporalpatterns(res)
                    res <- add_modelsensitivity2(res)
                  }

                  saveRDS(res, ftemp)
                }
              }
            }
          }
        }
      }
    }
  }
}

if (do_checks) {
  template_args <- if (do_targets) {
      list(
        only_wo_controls = TRUE, withControlsL. = std_design[["s1_wcontr"]],
        only_useadj_standardized = FALSE, withNonStandardizedL. = std_design[["s2_wnonnorm"]],
        fragment_sizes. = std_design[["s4_fragsize"]],
        cor_methods. = std_design[["s5_cormethod"]],
        cor_transforms. = std_design[["s6_cortransform"]],
        weight_methods. = std_design[["d7_weightmethod"]],
        effect_groups. = effect_groups,
        deffects = deffects, dmoderators = dmoderators, dir_res = dir_res_
      )

    } else {
      list(
        only_wo_controls = FALSE, withControlsL. = full_design[["s1_wcontr"]],
        only_useadj_standardized = FALSE, withNonStandardizedL. = full_design[["s2_wnonnorm"]],
        fragment_sizes. = full_design[["s4_fragsize"]],
        cor_methods. = full_design[["s5_cormethod"]],
        cor_transforms. = full_design[["s6_cortransform"]],
        weight_methods. = full_design[["d7_weightmethod"]],
        effect_groups. = effect_groups,
        deffects = deffects, dmoderators = dmoderators, dir_res = dir_res_
      )
    }

  do.call("run_metachecks", args = c(template_args,
    responses = list(responses_Hall)))
}

