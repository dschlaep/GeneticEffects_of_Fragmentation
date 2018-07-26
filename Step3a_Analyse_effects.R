Sys.setenv(LANG = "en_US.UTF-8")

###DONE:
  # Ar: use only standardized data
  # PLp, He&Sh, Fis: use with non-standardized data



library("metafor")
library("multcomp")
library("ape")

#--- SETTINGS
do_targets <- FALSE
do_rma <- TRUE

do_interactions <- c(TRUE, FALSE)
#do_interactions <- TRUE

do_adds <- list(
  add_confint = FALSE,
  add_BlombergsK = TRUE,
  add_modelsensitivity = TRUE,
  add_permutations = FALSE,
  add_gosh = FALSE
)
redos_rma <- list(
  read_data = FALSE,
  fit_metaanalysis = FALSE,
  add_confint = FALSE,
  add_BlombergsK = FALSE,
  add_modelsensitivity = FALSE,
  add_permutations = FALSE,
  add_gosh = FALSE
)

#--- Versions
do_check_outdatedversion <- FALSE
  fit_version <- numeric_version("2.2.3")
  # Version 2.0.0: first version with version control
  # Version 2.0.1: wrap content of 'add_publicationbias' in try()
  # Version 2.1.0: if weights == "winvar", then weights = NULL; prior: weights = dat[, "winvar"]
  # Version 2.2.0:
  #   - Added confidence intervals for tau^2, I2, and H2
  #   - Percentage of total amount of heterogeneity can be accounted for by including moderators
  #   - Use DerSimonian-Laird estimator if second attempt with REML fails ("because it is
  #     non-iterative and therefore faster and guaranteed to provide an estimate for
  #     tau^2 (the REML estimator requires iterative estimation and can occasionally
  #     fail to converge")
  #   - Added 3PSM (three-parameter selection model; Iyengar & Greenhouse, 1988; McShane,
  #     Böckenholt, & Hansen, 2016) to 'add_publicationbias': 3PSM "generally performs
  #     better than trim-and-fill, p-curve, p-uniform, PET, PEESE, or PET-PEESE, and that
  #     some of these other methods should typically not be used at all"
  #     (Carter et al. 2017 preprint)
  # Version 2.2.1:
  #   - Multi-level analyses corrected for dependencies: phylogeny, marker classes, articles
  # Version 2.2.2:
  #   - Multi-variate analysis for marker classes, i.e., allow sigma (residual heterogeneity)
  #     to vary between marker classes
  # Version 2.2.3:
  #   - `fit_multilevel`: calls `metafor:rma.mv` with a larger maximum of iteration steps
  #     and adds the quasi-Newton method "BFGS" from 'optim' as a fourth optimization option


# Paths
dir_prj <- "Prj04_GenetFragment_MetaAnalysis"
dir_ana <- file.path(dir_prj, "3_Analysis")

#--- Setup
source(file.path(dir_ana, "Step3_0-Analysis-settings.R"))
source(file.path(dir_ana, "Step3_1_Analysis-methods.R"))

source(file.path(dir_ana, "Step3_2_Analysis-Phylogeny.R"))




# Fit meta-analysis
fit_metaanalysis <- function(mdata, res, yin, vin, win) {

  #--- Random/Mixed-effect Models with Knapp & Hartung Adjustment: not accounting for
  #    non-independence
  res[["all"]][["fit"]][["uni"]] <- fit_uni(mdata = mdata, yin = yin, vin = vin, win = win)

  #--- Random/Mixed-effect Models with multi-level covariance structures accounting for
  #    marker categories and species and articles with multiple effects
  res[["all"]][["fit"]][["mv"]] <- fit_multilevel(sdata = mdata[["all"]], yin = yin, vin = vin,
    win = win, account_repeats = TRUE, account_phylo = FALSE)

  #--- Random/Mixed-effect Models with multi-level covariance structures accounting for
  #    marker categories, articles with multiple effects, and phylogenetic independence
  res[["pphylo"]][["fit"]][["wBL"]] <- fit_multilevel(sdata = mdata[["pphylo"]], yin = yin,
    vin = vin, win = win, account_repeats = TRUE, account_phylo = TRUE,
    tree_branchLengths = TRUE)
  res[["pphylo"]][["fit"]][["woBL"]] <- fit_multilevel(sdata = mdata[["pphylo"]], yin = yin,
    vin = vin, win = win, account_repeats = TRUE, account_phylo = TRUE,
    tree_branchLengths = FALSE)

  res[["aphylo"]][["fit"]][["wBL"]] <- NULL
  res[["aphylo"]][["fit"]][["woBL"]] <- fit_multilevel(sdata = mdata[["aphylo"]], yin = yin,
    vin = vin, win = win, account_repeats = TRUE, account_phylo = TRUE)

  res
}



#' Estimate confidence intervals
add_confint <- function(res) {
  # fit_uni
  res[["all"]][["fit"]][["uni"]][["nomod"]] <- calc_confint(res[["all"]][["fit"]][["uni"]][["nomod"]])
  res[["all"]][["fit"]][["uni"]][["mod"]] <- calc_confint(res[["all"]][["fit"]][["uni"]][["mod"]])

  # fit_multilevel
  res[["all"]][["fit"]][["mv"]][["nomod"]] <- calc_confint(res[["all"]][["fit"]][["mv"]][["nomod"]])
  res[["all"]][["fit"]][["mv"]][["mod"]] <- calc_confint(res[["all"]][["fit"]][["mv"]][["mod"]])
  res[["all"]][["fit"]][["mv"]][["mod_uni"]] <- calc_confint(res[["all"]][["fit"]][["mv"]][["mod_uni"]])

  res[["pphylo"]][["fit"]][["wBL"]][["nomod"]] <- calc_confint(res[["pphylo"]][["fit"]][["wBL"]][["nomod"]])
  res[["pphylo"]][["fit"]][["wBL"]][["mod"]] <- calc_confint(res[["pphylo"]][["fit"]][["wBL"]][["mod"]])
  res[["pphylo"]][["fit"]][["wBL"]][["mod_uni"]] <- calc_confint(res[["pphylo"]][["fit"]][["wBL"]][["mod_uni"]])

  res[["pphylo"]][["fit"]][["woBL"]][["nomod"]] <- calc_confint(res[["pphylo"]][["fit"]][["woBL"]][["nomod"]])
  res[["pphylo"]][["fit"]][["woBL"]][["mod"]] <- calc_confint(res[["pphylo"]][["fit"]][["woBL"]][["mod"]])
  res[["pphylo"]][["fit"]][["woBL"]][["mod_uni"]] <- calc_confint(res[["pphylo"]][["fit"]][["woBL"]][["mod_uni"]])

  res[["aphylo"]][["fit"]][["wBL"]][["nomod"]] <- calc_confint(res[["aphylo"]][["fit"]][["wBL"]][["nomod"]])
  res[["aphylo"]][["fit"]][["wBL"]][["mod"]] <- calc_confint(res[["aphylo"]][["fit"]][["wBL"]][["mod"]])
  res[["aphylo"]][["fit"]][["wBL"]][["mod_uni"]] <- calc_confint(res[["aphylo"]][["fit"]][["wBL"]][["mod_uni"]])

  res[["aphylo"]][["fit"]][["woBL"]][["nomod"]] <- calc_confint(res[["aphylo"]][["fit"]][["woBL"]][["nomod"]])
  res[["aphylo"]][["fit"]][["woBL"]][["mod"]] <- calc_confint(res[["aphylo"]][["fit"]][["woBL"]][["mod"]])
  res[["aphylo"]][["fit"]][["woBL"]][["mod_uni"]] <- calc_confint(res[["aphylo"]][["fit"]][["woBL"]][["mod_uni"]])

  res[["has"]][["confint"]] <- TRUE

  res
}


#' Estimate phylogenetic signal with Blomberg's K
add_BlombergsK <- function(data, res, yin = "E", yinse = "E_se") {
  datasets <- names(data)

  for (k in datasets) {
    if (!is.null(data[[k]][["dagg"]])) {
      if (!is.null(data[[k]][["tree"]])) {
        temp <- try(calc_BlombergsK(data = data[[k]][["dagg"]], yin = yin, yinse = yinse,
          tree = data[[k]][["tree"]]))
        if (!inherits(temp, "try-error")) {
          res[[k]][["add"]][["BlombergsK_wBL"]] <- temp
        }
      }

      if (!is.null(data[[k]][["tree_woBL"]])) {
        temp <- try(calc_BlombergsK(data = data[[k]][["dagg"]], yin = yin, yinse = yinse,
          tree = data[[k]][["tree_woBL"]]))
        if (!inherits(temp, "try-error")) {
          res[[k]][["add"]][["BlombergsK_woBL"]] <- temp
        }
      }
    }
  }

  res[["has"]][["BlombergsK"]] <- TRUE

  res
}



#------ Model sensitivity
add_modelsensitivity <- function(data, res, yin = "E", vin = "var") {
  datasets <- names(data)

  for (k in datasets) {
    # Fail-safe numbers
    res[[k]][["add"]][["Nfs"]] <- calc_failsafeN(data[[k]][["d"]], yin, vin)

    # Influence measures
    fits <- names(res[[k]][["fit"]])
    for (m in fits) {
      res[[k]][["add"]][[m]][["infl_nomod"]] <- calc_influence(res[[k]][["fit"]][[m]][["nomod"]])
      res[[k]][["add"]][[m]][["infl_mod"]] <- calc_influence(res[[k]][["fit"]][[m]][["mod"]])
    }
  }

  res[["has"]][["modelsensitivity"]] <- TRUE

  res
}



#--- Requiring heavy computation
add_permutations <- function(res, exact = FALSE, iter = 2000, permci = FALSE) {
  if (res[["uni"]][["mod"]][["has_fit"]] && nrow(res[["data"]][["d"]]) >= 6) {
    res[["uni"]][["mod"]][["permutest"]] <- metafor::permutest(
      res[["uni"]][["mod"]][["fit"]], exact = exact, iter = iter, permci = permci)
  }

  res[["has"]][["permutations"]] <- TRUE

  res
}

#' Calculate GOSH plots
#'
#' "Olkin et al. (2012) proposed the GOSH (graphical display of study heterogeneity) plot,
#' which is based on examining the results of a fixed-effects model in all possible
#' subsets of size 1, …, k of the k studies included in a meta-analysis. In a homogeneous
#' set of studies, the model estimates obtained this way should form a roughly symmetric,
#' contiguous, and unimodal distribution. On the other hand, when the distribution is
#' multimodal, then this suggests the presence of heterogeneity, possibly due to outliers
#' and/or distinct subgroupings of studies." (Viechtbauer, metafor v2.0-0)
add_gosh <- function(res, parallel = "multicore", ncpus = 20) {
  if (res[["uni"]][["mod"]][["has_fit"]] && nrow(res[["data"]][["d"]]) >= 6) {
    res[["uni"]][["mod"]][["gosh"]] <- metafor::gosh(res[["uni"]][["mod"]][["fit"]],
      parallel = parallel, ncpus = ncpus)
  }

  res[["has"]][["gosh"]] <- TRUE

  res
}

#' Prepare dataset structure for meta-analysis without and with plant and animal
#' phylogenetic corrections
#'
#' @return A list with three elements representing 3 datasets \itemize{
#'   \item all = all[ data]
#'   \item pphylo = p[lant]-phylo[genetic subset]
#'   \item aphylo = a[nimal]-phylo[genetic subset]
#' }
#' Each dataset element is a list with five elements \itemize{
#'   \item d = 'full' data
#'   \item dred = a reduced/aggregated dataset that matches the tree's species
#'   \item tree = a phylogenetic tree
#'   \item tree_woBL = the same phylogenetic tree but with uniformly fixed branch lengthts,
#'   \item mods = moderators
#'}
prepare_datasets <- function(deffects, dmoderators, moderators, response,
  interacting_moderator = NULL, wcontr, wnonnorm, indep, method, tf, win, tree_plants,
  tree_animals) {

  temp <- list(d = list(), dagg = list(), tree = list(), tree_woBL = list(), mods = list())
  data <- list(all = temp, pphylo = temp, aphylo = temp)

  #--- Data
  stopifnot(identical(dimnames(deffects)[[1]], as.character(dmoderators[, "ID_effect"])))
  temp_dat <- data.frame(
    deffects[, 2 - as.integer(wcontr), 2 - as.integer(wnonnorm), response, indep, method, tf, ],
    dmoderators, stringsAsFactors = FALSE)
  col_multilevels <- c("Marker_cat", "ID_article", "Species_resolved")
  col_complete <- unique(c(dimnames(deffects)[[length(dim(deffects))]], moderators,
    col_multilevels, interacting_moderator))
  id_use <- complete.cases(temp_dat[, col_complete, drop = FALSE])
  temp_dat <- temp_dat[id_use, , drop = FALSE]
  temp_dat[, "fSpecies"] <- factor(temp_dat[, "Species_resolved"])
  temp_dat[, "Species"] <- gsub(" ", "_", as.character(temp_dat[, "Species_resolved"]))

  # Columns to initialize random-effects
  temp_dat[, "betweenStudyVariance"] <- as.factor(temp_dat[, "ID_article"])
  temp_dat[, "phylogenyVariance"] <- as.factor(temp_dat[, "Species"])

  # Add interaction levels
  moderatorsX <- moderators
  if (!is.null(interacting_moderator)) {
    for (k in seq_along(moderators)) {
      mods <- c(moderators[k], interacting_moderator)
      moderatorsX[k] <- name_interaction(mods)
      temp_dat[, moderatorsX[k]] <- factor(apply(temp_dat[, mods], 1, paste, collapse = "_x_"))
    }
  }

  #--- Subset trees to available species
  data[["all"]][["tree"]] <- data[["all"]][["tree_woBL"]] <- NULL

  ids_plants <- !(tree_plants[["tip.label"]] %in% temp_dat[, "Species"])
  if (any(ids_plants)) {
    data[["pphylo"]][["tree"]] <- ape::drop.tip(tree_plants,
      tip = tree_plants[["tip.label"]][ids_plants])
    # Ignore branch lengths of tree
    temp <- data[["pphylo"]][["tree"]]
    temp[["edge.length"]][] <- 1
    data[["pphylo"]][["tree_woBL"]] <- temp

  } else {
    data[["pphylo"]][["tree"]] <- data[["pphylo"]][["tree_woBL"]] <- NULL
  }

  ids_animals <- !(tree_animals[["tip.label"]] %in% temp_dat[, "Species"])
  if (any(ids_animals)) {
    data[["aphylo"]][["tree"]] <- ape::drop.tip(tree_animals,
      tip = tree_animals[["tip.label"]][ids_animals])
    # Ignore branch lengths of tree
    temp <- data[["aphylo"]][["tree"]]
    temp[["edge.length"]][] <- 1
    data[["aphylo"]][["tree_woBL"]] <- temp

  } else {
    data[["aphylo"]][["tree"]] <- NULL
  }

  # Sort data according to species in subset trees and allowing for multiple occurrences
  # of species in dataset (species occur in trees at most once)
  temp_dat1p <- merge(data.frame(Species = data[["pphylo"]][["tree"]][["tip.label"]]),
    temp_dat, sort = FALSE)
  if (NROW(temp_dat1p) > 0) {
    temp_dat1p <- data.frame(temp_dat1p[, -1, drop = FALSE], Species = as.character(temp_dat1p[, 1]),
      stringsAsFactors = FALSE)
    temp_dat1p[, "with_phylo"] <- TRUE
  }
  temp_dat1a <- merge(data.frame(Species = data[["aphylo"]][["tree"]][["tip.label"]]),
    temp_dat, sort = FALSE)
  if (NROW(temp_dat1a) > 0) {
    temp_dat1a <- data.frame(temp_dat1a[, -1, drop = FALSE], Species = as.character(temp_dat1a[, 1]),
      stringsAsFactors = FALSE)
    temp_dat1a[, "with_phylo"] <- TRUE
  }

  ids2_sp <- !(temp_dat[, "Species"] %in%
    c(data[["pphylo"]][["tree"]][["tip.label"]], data[["aphylo"]][["tree"]][["tip.label"]]))
  temp_dat2 <- temp_dat[ids2_sp, , drop = FALSE]
  if (NROW(temp_dat2) > 0) {
    temp_dat2[, "with_phylo"] <- FALSE
  }

  #--- Complete dataset for analysis
  data[["all"]][["d"]] <- rbind(temp_dat1p, temp_dat1a, temp_dat2)
  data[["pphylo"]][["d"]] <- temp_dat1p
  data[["aphylo"]][["d"]] <- temp_dat1a

  #--- Reduced dataset to match phylogenetic tree
  data[["all"]][["dagg"]] <- agg_by_sp(data[["all"]][["d"]], data[["all"]][["tree"]])
  data[["pphylo"]][["dagg"]] <- agg_by_sp(data[["pphylo"]][["d"]], data[["pphylo"]][["tree"]])
  data[["aphylo"]][["dagg"]] <- agg_by_sp(data[["aphylo"]][["d"]], data[["aphylo"]][["tree"]])

  print(paste("overall dataset has nrow = ", nrow(data[["all"]][["d"]])))
  print(paste("plant-phylogenetic dataset has nrow = ", nrow(data[["pphylo"]][["d"]])))
  print(paste("animal-phylogenetic dataset has nrow = ", nrow(data[["aphylo"]][["d"]])))

  #--- Categories for multiple moderators
  data[["all"]][["mods"]] <- get_categories(data[["all"]][["d"]], moderatorsX)
  data[["pphylo"]][["mods"]] <- get_categories(data[["pphylo"]][["d"]], moderatorsX)
  data[["aphylo"]][["mods"]] <- get_categories(data[["aphylo"]][["d"]], moderatorsX)

  data
}


run_metaanalysis <- function(responses, moderators, interaction_wHabitat3 = FALSE,
  fragment_sizes. = fragment_sizes, withControlsL. = withControlsL, only_wo_controls = TRUE,
  withNonStandardizedL. = withNonStandardizedL, only_useadj_standardized = TRUE,
  cor_methods. = cor_methods, cor_transforms. = cor_transforms,
  weight_methods. = weight_methods, deffects, dmoderators, dir_res, dos, redos) {

  temp_withControlsL <- if (only_wo_controls) FALSE else withControlsL.

  for (ir in responses) {
    dtag <- paste0("Response_", ir)
    dir_data <- file.path(dir_res, "Data", dtag)
    dir.create(dir_data, recursive = TRUE, showWarnings = FALSE)
    dir_resout <- file.path(dir_res, "Results", dtag)
    dir.create(dir_resout, recursive = TRUE, showWarnings = FALSE)

    temp_withNonStandardizedL <- if (only_useadj_standardized || ir %in% c("Ar", "mA")) {
        # always use only standardized values for Ar and mA
        FALSE
      } else {
        withNonStandardizedL.
      }
    for (ig in moderators) {
      for (ii in fragment_sizes.) {
        for (ic in temp_withControlsL) {
          for (inn in temp_withNonStandardizedL) {
            for (im in cor_methods.) {
              temp_cor_transforms <- if (im == "pearson") {
                  cor_transforms.
                } else {
                  c("none", "ztransform")
                }
              for (it in temp_cor_transforms) {
                for (iw in weight_methods.) {
                  # File identification
                  ftag <- filetag_ID(interaction_wHabitat3, ii,
                    ic, only_wo_controls = FALSE, inn, im, it, iw)
                  ftag <- paste(ir, ig, "by", ftag, sep = "_")

                  print(paste(Sys.time(), "'run_metaanalysis':", ftag))

                  fdata <- file.path(dir_data, paste0(ftag, "_data.rds"))
                  ftemp <- file.path(dir_resout, paste0(ftag, "_output.rds"))
                  save_to_disk <- FALSE

                  do_fit_rma <- !file.exists(ftemp) || !file.exists(ftemp) ||
                    redos[["read_data"]] || redos[["fit_metaanalysis"]]

                  if (file.exists(fdata)) {
                    mdata <- readRDS(fdata)
                  }

                  if (!do_fit_rma && file.exists(ftemp)) {
                    temp <- readRDS(ftemp)
                    do_fit_rma <- do_fit_rma || !identical(temp[["version"]], fit_version)
                  }

                  if (do_fit_rma) {
                    mdata <- prepare_datasets(deffects, dmoderators,
                      moderators = ig, response = ir,
                      interacting_moderator = if (interaction_wHabitat3) "Habitat3" else NULL,
                      wcontr = ic, wnonnorm = inn, indep = ii, method = im, tf = it, win = iw,
                      tree_plants = dphyloplants, tree_animals = dphyloanimals)

                    dir.create(dirname(fdata), recursive = TRUE, showWarnings = FALSE)
                    saveRDS(mdata, fdata)

                    res <- new.env(parent = baseenv())
                    res <- fit_metaanalysis(mdata, res, yin = "E", vin = "var", win = iw)
                    save_to_disk <- TRUE

                  } else {
                    res <- readRDS(ftemp)
                  }

                  if (!("has" %in% names(res))) {
                    res[["has"]] <- list()
                  }

                  if (redos[["add_BlombergsK"]] && isTRUE(res[["has"]][["BlombergsK"]]) ||
                    dos[["add_BlombergsK"]] && !isTRUE(res[["has"]][["BlombergsK"]])) {

                    res <- add_BlombergsK(mdata, res, yin = "E", yinse = "E_se")
                    save_to_disk <- TRUE
                  }

                  if (redos[["add_modelsensitivity"]] && isTRUE(res[["has"]][["modelsensitivity"]]) ||
                    dos[["add_modelsensitivity"]] && !isTRUE(res[["has"]][["modelsensitivity"]])) {

                    res <- add_modelsensitivity(mdata, res, yin = "E")
                    save_to_disk <- TRUE
                  }

                  if (redos[["add_confint"]] && isTRUE(res[["has"]][["confint"]]) ||
                    dos[["add_confint"]] && !isTRUE(res[["has"]][["confint"]])) {

                    res <- add_confint(res)
                    save_to_disk <- TRUE
                  }


                  if (redos[["add_permutations"]] && isTRUE(res[["has"]][["permutations"]]) ||
                    dos[["add_permutations"]] && !isTRUE(res[["has"]][["permutations"]])) {

                    res <- add_permutations(res)
                    save_to_disk <- TRUE
                  }

                  if (redos[["add_gosh"]] && isTRUE(res[["has"]][["gosh"]]) ||
                    dos[["add_gosh"]] && !isTRUE(res[["has"]][["gosh"]])) {

                    res <- add_gosh(res)
                    save_to_disk <- TRUE
                  }

                  if (save_to_disk) {
                    dir.create(dirname(ftemp), recursive = TRUE, showWarnings = FALSE)

                    res[["version"]] <- fit_version
                    saveRDS(res, ftemp)

                    print(warnings())
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}



if (do_rma) {
  template_args <- if (do_targets) {
      c(
        list(only_wo_controls = TRUE),
        design_arguments[["args_target"]],
        list(deffects = deffects, dmoderators = dmoderators, dir_res = dir_res_,
        dos = do_adds, redos = redos_rma)
      )

    } else {
      c(
        list(only_wo_controls = FALSE),
        design_arguments[["args_full"]],
        list(deffects = deffects, dmoderators = dmoderators, dir_res = dir_res_,
        dos = do_adds, redos = redos_rma)
      )
    }

  for (intH3 in do_interactions) {
    # 0_Miscellaneous
    do.call("run_metaanalysis", args = c(template_args,
      responses = list(responses_Hall),
      moderators = list(moderators_Hmisc),
      interaction_wHabitat3 = list(intH3)))

    # 1_DecreasedGeneticDiversity
    do.call("run_metaanalysis", args = c(template_args,
      responses = list(responses_H1),
      moderators = list(moderators_H12),
      interaction_wHabitat3 = list(intH3)))

    # 2_IncreasedInbreeding
    do.call("run_metaanalysis", args = c(template_args,
      responses = list(responses_H2),
      moderators = list(moderators_H12),
      interaction_wHabitat3 = list(intH3)))

    # 3_TimeEffects
    do.call("run_metaanalysis", args = c(template_args,
      responses = list(responses_Hall),
      moderators = list(moderators_Htime),
      interaction_wHabitat3 = list(intH3)))
  }
}



#------ Check for outdated version of analysis objects
if (do_check_outdatedversion) {
  ftemp <- list.files(file.path(dir_res_, "Results"), pattern = "_output.rds",
    full.names = TRUE, recursive = TRUE)

  for (f in ftemp) {
    x <- readRDS(f)

    if (!(is.environment(x) && "version" %in% names(x) &&
        x[["version"]] >= fit_version)) {

      # move file to 'outdated' folder
      dir_outdated <- paste0(dirname(f), "_outdated")
      dir.create(dir_outdated, recursive = TRUE, showWarnings = FALSE)
      file.rename(from = f, to = file.path(dir_outdated,
        sub("_output.rds", "_output_outdated.rds", basename(f))))
    }
  }
}
