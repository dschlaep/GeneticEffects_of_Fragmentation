#--- SETTINGS
redo <- FALSE
do_ms <- TRUE
do_targets <- TRUE
do_popsize <- TRUE
interactions <- c(TRUE, FALSE)


# Paths
dir_prj1 <- "Prj04_GenetFragment_MetaAnalysis"
dir_ana <- file.path(dir_prj1, "3_Analysis")


#--- Setup
source(file.path(dir_ana, "Step5-0_Manuscript-methods.R"))


#--- Main figures
msres_metaanalysis <- function(responses, moderators, interaction_wHabitat3 = FALSE,
  fragment_sizes. = fragment_sizes, withControlsL. = withControlsL, only_wo_controls = TRUE,
  withNonStandardizedL. = withNonStandardizedL, only_useadj_standardized = TRUE,
  cor_methods. = cor_methods, cor_transforms. = cor_transforms,
  weight_methods. = weight_methods, dir_res, dir_out, ftag = "", panels_vertical = TRUE,
  ...) {

  #--- Load analysis outputs
  msdat <- msres_getdatatogether(responses, moderators, interaction_wHabitat3,
    fragment_sizes., withControlsL., only_wo_controls,
    withNonStandardizedL., only_useadj_standardized,
    cor_methods., cor_transforms., weight_methods., dir_res, dir_out = dir_ms_out, ftag)

  temp <- lab_responses[match(responses, lab_responses[, "code"]), "hypothesis"]
  msdat <- msres_addsensitivity(msdat, sided = temp, interaction_wHabitat3, dir_res, ftag)

  if (length(msdat[["modcats"]]) > 0) {
    #--- Gather data for plotting
    estimates <- c("beta", "ci.lb", "ci.ub", "pval")
    nm <- length(lab_fittypes[, "code"])

    temp <- array(0, dim = c(length(responses), length(msdat[["modcats"]]), nm),
      dimnames = list(responses, msdat[["modcats"]], lab_fittypes[, "code"]))

    dat_forest <- list(
      est = array(NA, dim = c(length(responses), length(msdat[["modcats"]]), nm,
        length(estimates)), dimnames = list(responses, msdat[["modcats"]],
        lab_fittypes[, "code"], estimates)),

      k_sample = temp,
      agree = temp
    )

    for (ir in seq_along(responses)) {
      for (ig in seq_along(moderators)) {
        x <- msdat[["res"]][moderators[ig], responses[ir]][[1]]

        fits <- list(
          all_un = x[["all"]][["fit"]][["uni"]][["mod"]],
          all_mv = x[["all"]][["fit"]][["mv"]][["mod"]],
          aphylo = x[["aphylo"]][["fit"]][["woBL"]][["mod"]],
          pphylo = x[["pphylo"]][["fit"]][["woBL"]][["mod"]])

        x <- msdat[["data"]][moderators[ig], responses[ir]][[1]]
        dats <- list(
          all_un = x[["all"]][["d"]],
          all_mv = x[["all"]][["d"]],
          aphylo = x[["aphylo"]][["d"]],
          pphylo = x[["pphylo"]][["d"]])

        agree <- msdat[["sensitivity"]]["agreement", moderators[ig], responses[ir]][[1]]
        colnames(agree) <- c("all_un", "all_mv", "aphylo", "pphylo")

        has_fits <- sapply(fits, function(f) f[["has_fit"]])

        if (any(has_fits)) {
          coefnames <- list()
          for (ft in lab_fittypes[, "code"]) if (has_fits[ft]) {
            temp <- names(fits[[ft]][["fit"]][["coef.na"]])[!fits[[ft]][["fit"]][["coef.na"]]]
            coefnames[[ft]] <- if (interaction_wHabitat3) {
                paste0(substr(temp, 2, 4), ig, substr(temp, 5, nchar(temp)))
              } else {
                temp <- substr(temp, 2, nchar(temp))
                for (i in seq_along(msdat[["moderators"]])) {
                  temp <- sub(msdat[["moderators"]][i],
                    paste0(msdat[["moderators"]][i], ig), temp)
                }
                temp
              }
          }
          coefnames_all <- sort(unique(unlist(coefnames)))

          for (ft in lab_fittypes[, "code"]) if (has_fits[ft]) {
            # Agreement of sensitivity analysis
            temp1 <- remove_modcatID(msdat[["modcats"]])
            temp2 <- remove_modcatID(dimnames(agree)[[1]])
            ids <- match(temp1, temp2, nomatch = 0)
            dat_forest[["agree"]][ir, ids > 0, ft] <- agree[ids, ft]

            # Sample sizes per level
            ids <- match(msdat[["modcats"]], coefnames[[ft]], nomatch = 0)
            dat_forest[["k_sample"]][ir, ids > 0, ft] <- colSums(fits[[ft]][["fit"]][["X"]])[ids]

            # Estimates per level
            temp <- fits[[ft]][["fit"]][estimates]
            for (k in estimates) {
              dat_forest[["est"]][ir, ids > 0, ft, k] <- as.vector(temp[[k]])[ids]
            }

          } else {
            ids <- match(msdat[["modcats"]], coefnames_all, nomatch = 0)

            dat_forest[["k_sample"]][ir, ids > 0, ft] <- if (NROW(dats[[ft]]) > 0) {
                # model fit failed
                NA
              } else {
                # no effects available
                0
              }
            dat_forest[["agree"]][ir, ids > 0, ft] <- dat_forest[["k_sample"]][ir, ids > 0, ft]
          }
        }
      }
    }

    #---
    labs_mod_levels <- clean_labels(unlist(msdat[["cats"]]), remove_digits = TRUE,
      break_lines = FALSE, remove_interaction = isTRUE(grepl("habitat",
        msdat[["moderators"]], ignore.case = TRUE)))
    mod_labs <- clean_labels(msdat[["moderators"]], remove_digits = TRUE,
      break_lines = FALSE)


    #--- Plot main figure
    dir <- file.path(dir_out, "Figs_MainForest")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)

    dim_fig <- if (panels_vertical) c(length(responses), 1) else c(1, length(responses))
    cex <- 1
    ctemp <- sum(msdat[["n_cats_max"]])
    if (ctemp < 15) {
      cex_txt <- 0.9
      h.panel <- ctemp * 0.75
    } else {
      cex_txt <- 0.75
      h.panel <- ctemp * 0.6
    }
    h.edgeL <- max(0.05, min(0.3, 0.3 - 0.0125 * ctemp))
    h.edgeU <- max(0.05, min(0.3, 0.3 - 0.01 * ctemp))
    w.panel <- 3
    temp <- max(strwidth(c(labs_mod_levels, "(pc = 00, 00, 00, 00%)"),
      units = "inches"))
    w.edgeL <- ceiling(105 * cex_txt * temp) / 100
    w.edgeR <- 0.1

    pdf(height = h.edgeL + h.panel * dim_fig[1] + h.edgeU,
      width = w.edgeL + w.panel * dim_fig[2] + w.edgeR,
      file = file.path(dir, paste0("Fig_MainForest_", ftag, "_",
      msdat[["tag_fix"]], ".pdf")))

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
    par_prev <- par(mgp = c(1, 0, 0), mar = c(0.5, 0.5, 0.5, 0.2), tcl = 0.3, cex = cex,
      xaxs = "i")

    # Panels
    temp <- dat_forest[["est"]][, , , c("ci.lb", "ci.ub")]
    xlims <- if (all(is.na(temp))) {
        c(-1, 1)
      } else {
        temp <- range(temp, na.rm = TRUE)
        c(max(-1.75, min(-1.25, temp[1])), min(1.75, max(1.25, temp[2])))
      }

    ys_delta <- rep(0.5 * rev(seq_along(msdat[["cats"]]) - 1), times = msdat[["n_cats_max"]])
    ys <- rev(seq_along(msdat[["modcats"]])) + ys_delta
    temp <- cumsum(c(0, msdat[["n_cats_max"]])) + c(0, 0.5 * seq_along(msdat[["n_cats_max"]]))
    ys_starts <- max(ys) + 0.75 - temp[-length(temp)]
    ylims <- c(min(ys) - 0.5, max(ys_starts))

    for (ir in seq_along(responses)) {
      # New panel
      plot(1, type = "n", xlim = xlims, xlab = "", ylim = ylims,
        ylab = "", axes = FALSE)

      # Axes
      do_xaxis <- if (panels_vertical) ir == length(responses) else TRUE
      if (do_xaxis) {
        axis(side = 1, pos = ylims[1])
        eff_lab <- paste0(
          if (identical(cor_transforms., "ztransform")) {
            "z-transformed "
          } else if (identical(cor_transforms., "unbiased")) {
            "unbiased "
          } else "",
          capwords(cor_methods.), "'s r")

        #mtext(side = 1, text = eff_lab, adj = 0.5, line = -1)
        text(x = 0, y = ylims[1] - 2.75 * strheight(eff_lab), labels = eff_lab,
          adj = c(0.5, 0), xpd = NA)
      }

      # y-Axis: sample sizes for each fit_type
      ks <- apply(dat_forest[["k_sample"]][responses[ir], , ], 1, function(k) {
          paste0("(n = ", paste0(k, collapse = ", "), ")")
        })
      ags <- apply(round(100 * dat_forest[["agree"]][responses[ir], , ]), 1, function(k) {
        paste0("(pc = ", paste0(k, collapse = ", "), "%)")
      })

      do_yaxis <- panels_vertical || ir == 1L
      if (do_yaxis) {
        # y-Axis: levels of moderators with sample sizes for each fit_type
        temp <- paste0(labs_mod_levels, "\n", ks, "\n", ags)
        mtext(side = 2, at = ys, text = temp, cex = cex_txt, las = 1)

        # y-Axis: moderators
        if (length(msdat[["moderators"]]) > 1) {
          atx <- xlims[1] - 2 - if (any(grepl("\n", mod_labs))) 0.5 else 0
          #temp <- (w.edgeL + par("din")[1] * par("plt")[1]) * diff(par("usr")[1:2]) / par("pin")[1]
          temp <- w.edgeL * diff(par("usr")[1:2]) / par("pin")[1]
          atx <- max(atx, par("usr")[1] - temp)
          if (FALSE) {
            tempy <- max(ys) + 0.5 - cumsum(c(0, msdat[["n_cats_max"]]))
            aty <- tempy[-length(tempy)] + diff(tempy) / 2
            text(x = atx, y = aty, labels = mod_labs, cex = cex_txt, srt = 90, adj = 0.5,
              xpd = NA)
          }
          aty <- ys_starts - 0.2
          text(x = atx, y = aty, labels = mod_labs, font = 2, cex = cex_txt, adj = 0,
            xpd = NA)
        }

      } else {
        temp <- paste0(ks, "\n", ags)
        mtext(side = 2, at = ys, text = temp, cex = cex_txt, adj = 0.5, las = 1)
      }

      # Left x-position for annotations
      datx <- xlims[1] - if (do_yaxis) 1.5 else 0
      # Determine left-most x-value of figure region in user coordinates
      finx_usr <- par("usr")[1] - par("plt")[1] * par("fin")[1] *
        diff(par("usr")[1:2]) / par("pin")[1]
      datx <- max(finx_usr, datx)

      # Panel identification text
      if (length(responses) > 1) {
        temp <- lab_responses[responses[ir] == lab_responses[, "code"], "out"]
        mtext(side = 3, text = paste0("(", letters[ir], ") ", temp), font = 2, at = datx,
          line = -1, adj = 0)
      }

      # Zero effect line
      segments(x0 = 0, y0 = ylims[1], y1 = ylims[2], lwd = 2, lty = 1, col = "gray")

      # Add for each fit_type
      for (ift in seq_along(lab_fittypes[, "code"])) {
        d <- dat_forest[["est"]][responses[ir], , lab_fittypes[, "code"][ift], ]
        Ns <- dat_forest[["k_sample"]][responses[ir], , lab_fittypes[, "code"][ift]]
        at_ys <- ys + (nm / 2 - ift) / (2 * nm)

        # Add mean and 95%-CI
        points(x = d[, "beta"], y = at_ys, pch = 4, lwd = 2,
          col = lab_fittypes[ift, "col"])
        has_n <- Ns >= 3
        arrows(x0 = d[has_n, "ci.lb"], x1 = d[has_n, "ci.ub"], y0 = at_ys[has_n],
          angle = 90, code = 3, length = 1 / (4 * nm + 3), lwd = 2,
          col = lab_fittypes[ift, "col"])
      }

      # Separating lines between moderators
      if (length(ys_starts) > 1) {
        segments(x0 = datx, x1 = xlims[2], y0 = ys_starts, lty = 2, col = "black", xpd = NA)
      }

      # Legend
      if (TRUE && ir == length(responses)) {
        legend(x = xlims[2], y = ys_starts[1] - 0.1, xjust = 1, yjust = 0.5, xpd = NA,
          bg = "white", ncol = 1, legend = lab_fittypes[, "out"], cex = cex_txt,
          col = lab_fittypes[, "col"], pch = 4, lwd = 3, pt.cex = cex)
      }

      if (FALSE && ir == 1) {
        legend(x = -0.75, y = ys_starts[1] - 0.1, xjust = 0, yjust = 0, xpd = NA,
          bg = "white", ncol = 2, legend = lab_fittypes[, "out"], cex = cex_txt,
          col = lab_fittypes[, "col"], pch = 4, lwd = 2, pt.cex = cex)
      }
    }

    par(par_prev)
    dev.off()

  } else {
    print(paste("No data to plot for:", ftag, "/", msdat[["tag_fix"]]))
  }
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

  if (do_popsize) {
    template_args[["fragment_sizes."]] <- "Pop_size_n-published"
  }

  template_args <- expand.grid(template_args, stringsAsFactors = FALSE)

  args_along <- list(
      dir_out = list(dir_res0, dir_res12, dir_res12, dir_res12, dir_res3)
    )

  for (k in seq_len(NROW(template_args))) {
    do_ms_loops(fun = "msres_metaanalysis", args = template_args[k, ],
      args_along = args_along, do_targets = do_targets, do_interactions = interactions)
  }
}

