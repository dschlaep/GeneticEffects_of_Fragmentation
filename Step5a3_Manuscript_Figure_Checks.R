# Overall checks for simple meta-analysis models without including moderators or accounting for dependencies

#--- SETTINGS
redo <- FALSE
do_ms <- TRUE
do_targets <- TRUE

library("metafor")


# Paths
dir_prj1 <- "Prj04_GenetFragment_MetaAnalysis"
dir_ana <- file.path(dir_prj, "3_Analysis")

#--- Setup
source(file.path(dir_ana, "Step5-0_Manuscript-methods.R"))


# Plot functions
# Plot cumulative meta-analysis (in the order of publication year)
plot_cumul_ma <- function(data, responses, panels_vertical, dir_out, ftag_fix,
    cor_transforms., cor_methods.) {
  pub_years <- NULL
  for (ir in seq_along(responses)) {
    pub_years <- c(pub_years, data[[responses[ir]]][["cumtime"]][["yrs_sorted"]])
  }
  pub_years <- seq(from = min(pub_years), to = max(pub_years), by = 1)

  x_names <- c("estimate", "ci.lb", "ci.ub")
  res_cumtime <- array(list(), dim = c(length(pub_years), length(responses)),
    dimnames = list(pub_years, responses))

  for (ir in seq_along(responses)) {
    x <- data[[responses[ir]]][["cumtime"]]

    for (y in pub_years) {
      ids <- x[["yrs_sorted"]] %in% y
      if (any(ids)) {
        x_data <- matrix(NA, nrow = sum(ids), ncol = 1 + length(x_names),
          dimnames = list(NULL, c("Year", x_names)))
        x_data[, "Year"] <- y

        for (k in x_names) {
          x_data[, k] <- x[["fit"]][[k]][ids]
        }

        res_cumtime[as.character(y), responses[ir]] <- list(x_data)
      }
    }
  }

  # Prepare data for plot
  yt <- rev(seq_along(pub_years))

  res <- list()
  for (ir in seq_along(responses)) {
    x <- res_cumtime[, responses[ir]]
    N <- sapply(x, NROW)
    xdf <- data.frame(do.call(rbind, x))
    N_per_yr <- rep(N, N) # number of effects per year
    i_within_yr <- unlist(lapply(N, seq_len)) # runningn count of effect within each year
    y_per_yr <- rep(yt, N)
    xdf[, "y"] <- y_per_yr +
      (N_per_yr - i_within_yr - (N_per_yr - 1) / 2) / (2 * N_per_yr)
    res[[responses[ir]]] <- xdf
  }


  # Prepare plot
  ylims <- c(0.5, length(yt) + 0.5)

  dim_fig <- if (panels_vertical) c(length(responses), 1) else c(1, length(responses))
  cex <- 1; cex_txt <- 0.7
  h.panel <- 0.25 * length(pub_years); h.edgeL <- 0.5; h.edgeU <- 0.3
  w.panel <- 3; w.edgeL <- 0.5; w.edgeR <- 0.25

  # Plot
  ftemp0 <- file.path(dir_out, "Figs_Checks")
  dir.create(ftemp0, recursive = TRUE, showWarnings = FALSE)

  pdf(height = h.edgeL + h.panel * dim_fig[1] + h.edgeU,
    width = w.edgeL + w.panel * dim_fig[2] + w.edgeR,
    file = file.path(ftemp0, paste0("Fig_CumulTime_", ftag_fix, ".pdf")))

  temp <- c(
    rep(0, 1 + dim_fig[1] + 1),
    sapply(seq_len(dim_fig[2]), function(k)
      c(0, (k - 1) * dim_fig[1] + seq_len(dim_fig[1]), 0)),
    rep(0, 1 + dim_fig[1] + 1))
  temp <- matrix(temp, nrow = 1 + dim_fig[1] + 1, ncol = 1 + dim_fig[2] + 1,
    byrow = FALSE)

  layout(temp,
    heights = c(h.edgeU, rep(h.panel, times = dim_fig[1]), h.edgeL),
    widths = c(w.edgeL, rep(w.panel, times = dim_fig[2]), w.edgeR))
  par_prev <- par(mgp = c(1, 0.25, 0), mar = c(0.5, 0.5, 0.5, 0.2), tcl = 0.3, cex = cex)

  for (ir in seq_along(responses)) {
    xdf <- res[[responses[ir]]]

    xlims <- range(xdf[, c("ci.lb", "ci.ub")], na.rm = TRUE)
    xlims <- c(max(-2, xlims[1]), min(2, xlims[2]))

    plot(xdf[,  "estimate"], xdf[, "y"], pch = 4, lwd = 1, xlim = xlims, ylim = ylims,
      ann = FALSE, axes = FALSE)
    arrows(x0 = xdf[, "ci.lb"], x1 = xdf[, "ci.ub"], y0 = xdf[, "y"], code = 3,
      angle = 90,length = grid::unit(0.05, "lines"), lwd = 1)
    points(xdf[,  "estimate"], xdf[, "y"], pch = 4, lwd = 1, col = "darkgray")
    segments(x0 = 0, y0 = ylims[1], y1 = ylims[2], col = "gray")

    # Axes
    do_xaxis <- if (panels_vertical) ir == length(responses) else TRUE
    if (do_xaxis) {
      axis(side = 1)
      eff_lab <- paste0(
        if (identical(cor_transforms., "ztransform")) {
          "z-transformed "
        } else if (identical(cor_transforms., "unbiased")) {
          "unbiased "
        } else "",
        capwords(cor_methods.), "'s r")

      mtext(side = 1, text = eff_lab, adj = 0.5, line = 1.25)
    }

    # y-Axis: sample sizes for each fit_type
    do_yaxis <- panels_vertical || ir == 1L
    axis(side = 2, at = yt, labels = if (do_yaxis) pub_years else FALSE, las = 1)

    # Panel identification text
    temp <- lab_responses[responses[ir] == lab_responses[, "code"], "out"]
    mtext(side = 3, text = paste0("(", letters[ir], ") ", temp), font = 2,
      at = xlims[1], adj = 0)
  }

  par(par_prev)
  dev.off()
}

# Some diagnostic plots
plot_diagnostics <- function(data, responses, dir_out, ftag_fix) {

  for (ir in seq_along(responses)) {
    dtag <- paste0("Response_", responses[ir])
    x <- data[[responses[ir]]]

    ftemp <- file.path(dir_out, "Figs_Checks", paste0("Fig_Diagnostics_", dtag, "_", ftag_fix, ".pdf"))

    pdf(height = 7, width = 7, file = ftemp)
    par_prev <- par(mfrow = c(2, 2), mar = c(3, 3, 1, 1), mgp = c(1.5, 0, 0), tcl = 0.3)

    #--- Quantile-Quantile plot with pseudo confidence envelope
    metafor::qqnorm.rma.uni(x[["fit"]], bonferroni = TRUE, level = 0.95, main = "")
    mtext(side = 3, text = paste0("(", letters[1], ") normal q-q plot"), font = 2, adj = 0)

    #--- Baujat et al. (2002): diagnostic plot to detect sources of heterogeneity
    # x-axis = squared Pearson residual of a study; for fixed-effect this is the individual
    #   contribution to the overall QE-test statistic
    # y-axis = standardized squared difference between the predicted/fitted value for the
    #   study with and without the study included in the model fitting; i.e., individual
    #   influence of a study
    bj <- metafor::baujat(x[["fit"]], symbol = 16, main = "")
    mtext(side = 3, text = paste0("(", letters[2], ") Baujat plot"), font = 2, adj = 0)

    #--- Galbraith plots
    # vertical axis corresponds to standardized values, it is referred to as the z-axis
    # within this function. On the right hand side of the plot, an arc is drawn (referred to
    # as the y-axis within this function) corresponding to the individual observed effect
    # sizes or outcomes. A line projected from (0,0) through a particular point within the
    # plot onto this arc indicates the value of the individual observed effect size or
    # outcome for that point.
    metafor::galbraith(x[["fit"]], main = "")
    mtext(side = 3, text = paste0("(", letters[3], ") Galbraith plot"), font = 2, adj = 0)


    #--- Trim & fill
    metafor::funnel(x[["trim-fill"]], yaxis = "sei", addtau2 = TRUE,
      level = c(95, 99, 99.9), back = "white", shade = c("gray60", "gray75", "gray90"),
      main = "")
    mtext(side = 3, text = paste0("(", letters[4], ") trim & fill funnel plot"), font = 2, adj = 0)

    par(par_prev)
    dev.off()
  }
}


# Tabulate outliers
table_outliers <- function(data, responses, dir_out, ftag_fix, clean = TRUE) {

  ftemp0 <- file.path(dir_out, "Tables")
  dir.create(ftemp0, recursive = TRUE, showWarnings = FALSE)

  ftable <- file.path(ftemp0, paste0("Table_Outliers_", ftag_fix, ".csv"))
  if (clean) unlink(ftable)

  for (ir in seq_along(responses)) {
    dtag <- paste0("Response_", responses[ir])
    x <- data[[responses[ir]]]

    # Calculate Baujat's plot statistics
    bj <- metafor::baujat(x[["fit"]])
    dev.off()
    bj_xy <- bj[, c("x", "y")]
    bj_mahal <- mahalanobis(bj_xy, colMeans(bj_xy), cov(bj_xy))
    bj_box <- boxplot.stats(bj_mahal)[["out"]]
    ids <- as.integer(names(bj_box))

    Baujat_outlier <- seq_len(NROW(x[["data"]])) %in% ids

    #--- Collect outliers
    dat_outliers <- data.frame(x[["data"]], response = responses[ir],
      infl = x[["infl"]][["infl"]][[1]][["inf"]],
      infl_outlier = x[["infl"]][["infl"]][[1]][["is.infl"]],
      Baujat = bj_xy,
      Baujat_mahalanobis = bj_mahal, Baujat_outlier = Baujat_outlier,
      Outlier = Baujat_outlier & bj_mahal > 10 & x[["infl"]][["infl"]][[1]][["is.infl"]]
    )

    if (file.exists(ftable)) {
      write.table(dat_outliers, file = ftable, sep = ",", append = TRUE,
        col.names = FALSE, row.names = FALSE)
    } else {
      write.csv(dat_outliers, file = ftable, row.names = FALSE)
    }
  }

  invisible(ftable)
}


# Tabulate publication bias (fail-safe numbers, 3-PSM)
table_pubbias <- function(data, responses, dir_out, ftag_fix, clean = TRUE) {

  ftable <- file.path(dir_out, "Tables",
    paste0("Table_PublicationBias_", ftag_fix, ".csv"))
  if (clean) unlink(ftable)

  Nfs_names <- unique(unlist(lapply(data, function(x)
    sapply(x[["Nfs"]], function(nfs) nfs[["type"]]))))
  Nfs_lengths <- length(Nfs_names)

  cols <- c("Response",
    paste("TrimFill_MissingStudies",
      c("Estimator", "Side", "N", "se", "pvalue"), sep = "_"),

    paste("FailSafeN", Nfs_names, rep(c("N", "pvalue"), each = Nfs_lengths),
      sep = "_"),

    paste("FunnelPlotAsymmetry",
      c(paste("RankCorTest_BeggMazumbar1994", c("KendallsTau", "pvalue"),
        sep = "_"),
        paste("MixedEffectRegTest_Egger",
          c("predictor", "dfs", "zval", "pvalue"), sep = "_")
      ), sep = "_")
  )

  dat_pubbias <- data.frame(matrix(NA, nrow = length(responses),
    ncol = length(cols), dimnames = list(NULL, cols)))

  for (ir in seq_along(responses)) {
    x <- data[[responses[ir]]]

    dat_pubbias[ir, "Response"] <- responses[ir]

    temp <- x[["trim-fill"]]
    dat_pubbias[ir, "TrimFill_MissingStudies_Estimator"] <- temp[["k0.est"]]
    dat_pubbias[ir, "TrimFill_MissingStudies_Side"] <- temp[["side"]]
    dat_pubbias[ir, "TrimFill_MissingStudies_N"] <- temp[["k0"]]
    dat_pubbias[ir, "TrimFill_MissingStudies_se"] <- temp[["se.k0"]]
    dat_pubbias[ir, "TrimFill_MissingStudies_pvalue"] <- temp[["p.k0"]]

    temp <- x[["Nfs"]]
    for (k in seq_len(Nfs_lengths)) {
      temp_col <- paste0("FailSafeN_", temp[[Nfs_names[k]]][["type"]])
      dat_pubbias[ir, paste0(temp_col, "_N")] <- temp[[Nfs_names[k]]][["fsnum"]]
      dat_pubbias[ir, paste0(temp_col, "_pvalue")] <- temp[[Nfs_names[k]]][["pval"]]
    }


    temp <- x[["3psm"]]
if (any(!is.na(temp[, -(1:2)]))) {
  print("3psm is not NA yet not implemented")
  print(ftag_fix)
  print(temp)
}

    temp <- x[["begg-mazumbar"]]
    dat_pubbias[ir, "FunnelPlotAsymmetry_RankCorTest_BeggMazumbar1994_KendallsTau"] <- temp[["tau"]]
    dat_pubbias[ir, "FunnelPlotAsymmetry_RankCorTest_BeggMazumbar1994_pvalue"] <- temp[["pval"]]


    temp <- x[["eggers-test-mixed"]]
    dat_pubbias[ir, "FunnelPlotAsymmetry_MixedEffectRegTest_Egger_predictor"] <- temp[["predictor"]]
    dat_pubbias[ir, "FunnelPlotAsymmetry_MixedEffectRegTest_Egger_dfs"] <- temp[["dfs"]]
    dat_pubbias[ir, "FunnelPlotAsymmetry_MixedEffectRegTest_Egger_zval"] <- temp[["zval"]]
    dat_pubbias[ir, "FunnelPlotAsymmetry_MixedEffectRegTest_Egger_pvalue"] <- temp[["pval"]]
  }

  write.csv(dat_pubbias, file = ftable, row.names = FALSE)

  invisible(ftable)
}


# Main function
msplot_metachecks <- function(responses, fragment_sizes. = fragment_sizes,
  withControlsL. = withControlsL, only_wo_controls = TRUE,
  withNonStandardizedL. = withNonStandardizedL, only_useadj_standardized = TRUE,
  cor_methods. = cor_methods, cor_transforms. = cor_transforms,
  weight_methods. = weight_methods, effect_groups. = effect_groups, dir_res,
  dir_out, panels_vertical = FALSE) {

  temp_withControlsL <- if (only_wo_controls) FALSE else withControlsL.

  # Only 'responses' allowed to vary
  stopifnot(lengths(list(fragment_sizes., temp_withControlsL, cor_methods.,
    cor_transforms., weight_methods., effect_groups.)) == 1L)

  # File identification
  temp_withNonStandardizedL <- !only_useadj_standardized
  ftag_fix <- filetag_ID(interaction_wHabitat3 = FALSE, fragment_sizes.,
    withControlsL., only_wo_controls, temp_withNonStandardizedL, cor_methods.,
    cor_transforms., weight_methods.)
  ftag_fix <- paste(ftag_fix, "by", effect_groups., sep = "_")

  #--- Load analysis outputs
  data <- list()

  for (ir in seq_along(responses)) {
    dtag <- paste0("Response_", responses[ir])
    dir_data <- file.path(dir_res, "Data", dtag)
    dir_resout <- file.path(dir_res, "Results", dtag)

    temp_withNonStandardizedL <- if (only_useadj_standardized || responses[ir] %in% c("Ar", "mA")) {
        # always use only standardized values for Ar and mA
        FALSE
      } else {
        withNonStandardizedL.
      }

    ftag <- filetag_ID(interaction_wHabitat3 = FALSE, fragment_sizes.,
      withControlsL., only_wo_controls, temp_withNonStandardizedL, cor_methods.,
      cor_transforms., weight_methods.)

    ftag <- paste(responses[ir], "by", ftag, "by", effect_groups., sep = "_")
    fdata <- file.path(dir_resout, paste0(ftag, "_checks.rds"))

    if (file.exists(fdata)) {
      data[[responses[ir]]] <- readRDS(fdata)
    } else {
      print(paste("File not found:", shQuote(basename(fdata))))
    }
  }


  #--- Create output
  plot_cumul_ma(data, responses, panels_vertical, dir_out, ftag_fix,
    cor_transforms., cor_methods.)

  plot_diagnostics(data, responses, dir_out, ftag_fix)

  table_outliers(data, responses, dir_out, ftag_fix)

  table_pubbias(data, responses, dir_out, ftag_fix)
}



if (do_ms) {
  template_args <- if (do_targets) {
      c(
        list(only_wo_controls = TRUE),
        design_arguments[["args_target"]],
        list(effect_groups. = effect_groups,
          dir_res = dir_res_, dir_out = dir_res0)
      )

    } else {
      c(
        list(only_wo_controls = FALSE),
        design_arguments[["args_full"]],
        list(effect_groups. = effect_groups,
          dir_res = dir_res_, dir_out = dir_res0)
      )
    }

  template_args <- expand.grid(template_args, stringsAsFactors = FALSE)

  for (k in seq_len(NROW(template_args))) {
    do.call("msplot_metachecks", args = c(template_args[k, ],
      responses = list(responses_Hall)))
  }
}

