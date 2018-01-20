
#--- SETTINGS
redo <- FALSE
do_ms <- TRUE
do_targets <- TRUE
levels_byHabitat3 <- c(FALSE, TRUE)


# Paths
dir_prj1 <- "Prj04_GenetFragment_MetaAnalysis"
dir_ana <- file.path(dir_prj1, "3_Analysis")


#--- Setup
source(file.path(dir_ana, "Step5-0_Manuscript-methods.R"))

fws <- 5 # figwidths[1]


#--- Main map
msmap_metaanalysis <- function(responses, fragment_sizes. = fragment_sizes,
  withControlsL. = withControlsL, only_wo_controls = TRUE,
  withNonStandardizedL. = withNonStandardizedL, only_useadj_standardized = TRUE,
  cor_methods. = cor_methods, cor_transforms. = cor_transforms,
  weight_methods. = weight_methods, deffects, dmoderators, byHabitat3 = FALSE,
  dir_res, dir_out, yin = "E", vin = "var") {

  temp_withControlsL <- if (only_wo_controls) FALSE else withControlsL.

  # Only 'responses' allowed to vary
  stopifnot(lengths(list(fragment_sizes., temp_withControlsL, withNonStandardizedL.,
    cor_methods., cor_transforms., weight_methods.)) == 1L)

  temp1 <- if (temp_withControlsL) "withControls" else "withoutControls"
  temp_cor_transforms <- if (cor_methods. == "pearson") {
      cor_transforms.
    } else {
      "ztransform"
    }
  temp2 <- if (!only_useadj_standardized) "withNonStandardized" else "onlyStandardized"

  tag_fix <- paste(fragment_sizes., temp1, temp2, paste0("COR", cor_methods., "-",
    temp_cor_transforms), weight_methods., sep = "_")
  if (byHabitat3) {
    tag_fix <- paste0(tag_fix, "_", "byHabitat3")
  }

  ftemp1 <- file.path(dir_out, "Figs_Maps", paste0("Fig_MainMap_", tag_fix, "_map.pdf"))
  ftemp2 <- file.path(dir_out, "Figs_Maps", paste0("Fig_MainMap_", tag_fix, "_map_perResponse.pdf"))

  # Data: locations
  locs <- dmoderators[, c("Effect_longitude", "Effect_latitude", "Habitat3")]
  locs[, "ID_locs"] <- apply(locs[, c("Effect_longitude", "Effect_latitude")], 1,
    function(x) paste(round(x, 0), collapse = "_"))
  dlocs <- data.frame(locs, matrix(0, nrow = nrow(locs), ncol = length(responses)))
  colnames(dlocs)[-seq_len(ncol(locs))] <- responses

  dlocs_uni <- locs[!duplicated(locs[, c("Effect_longitude", "Effect_latitude")]), ]
  dlocs_uni <- data.frame(dlocs_uni, matrix(0, nrow = nrow(dlocs_uni),
    ncol = length(responses)))
  colnames(dlocs_uni)[-seq_len(ncol(locs))] <- responses

  # Data: response
  resp <- deffects[, 2 - as.integer(temp_withControlsL), , responses, fragment_sizes.,
    cor_methods., temp_cor_transforms, c(yin, vin, weight_methods.)]

  for (ir in responses) {
    temp_withNonStandardizedL <- if (only_useadj_standardized || ir %in% c("Ar", "mA")) {
        # always use only standardized values for Ar and mA
        FALSE
      } else {
        withNonStandardizedL.
      }

    dlocs[, ir] <- complete.cases(resp[, 2 - as.integer(temp_withNonStandardizedL), ir, ])
    temp <- aggregate(dlocs[, ir], by = list(locs[, "ID_locs"]), sum)
    ids <- match(dlocs_uni[, "ID_locs"], temp[, 1], nomatch = 0)

    dlocs_uni[ids > 0, ir] <- temp[ids, 2]
  }

  dlocs[, "Response_any"] <- apply(dlocs[, responses], 1, sum)
  dlocs_uni[, "Response_any"] <- apply(dlocs_uni[, responses], 1, sum)
  dlocs_uni2 <- dlocs_uni[dlocs_uni[, "Response_any"] > 0, ]


  # Maps
  library("maps")
  aspr <- fws / 6

  temp_dlocs <- dlocs[dlocs[, "Response_any"] > 0, ]
  xlims <- range(temp_dlocs[, "Effect_longitude"])
  ylims <- range(temp_dlocs[, "Effect_latitude"])
  bdw <- min(c(diff(xlims), diff(ylims))) / 25
  fcols <- colorRampPalette(c("white", blues9[1:3], blues9[6], "purple3", "purple4"))
  levels_hab3 <- if (byHabitat3) levels(dlocs[, "Habitat3"]) else NA

  pdf(height = aspr * 2.75 * if (byHabitat3) 2 else 1, width = fws, file = ftemp1)
  par_prev <- if (byHabitat3) {
      par(mfrow = c(2, 1), mar = c(2.25, 2, 0.1, 0.1))
    } else {
      par(mar = c(2, 2, 0.1, 0.1))
    }
  par(mgp = c(1, 0.2, 0), tcl = 0.5, cex = 1)

  for (ip in seq_along(levels_hab3)) {
    iuse <- temp_dlocs
    iuse_uni <- dlocs_uni2
    if (byHabitat3) {
      iuse <- iuse[as.character(temp_dlocs[, "Habitat3"]) == levels_hab3[ip], ]
      iuse_uni <- iuse_uni[as.character(dlocs_uni2[, "Habitat3"]) == levels_hab3[ip], ]
    }

    with(iuse, plot(Effect_longitude, Effect_latitude, type = "n",
      asp = 1, xlim = xlims, ylim = ylims,
      xlab = if (ip == length(levels_hab3)) "Longitude" else "", ylab = "Latitude"))
    with(iuse, smoothScatter(Effect_longitude, Effect_latitude, add = TRUE,
      bandwidth = bdw, nbin = c(720, 180), xlim = c(-180, 180), ylim = c(-90, 90),
      nrpoints = 0, colramp = fcols))
    map(add = TRUE)
    with(iuse, points(Effect_longitude, Effect_latitude,
      col = "orange", lwd = 1, cex = 0.5, pch = 4))
    abline(h = 0, col = "gray")
    if (byHabitat3) {
      legend("topleft", box.col = "white", bg = "white", text.font = 2, inset = -0.04,
        legend = paste0("(", letters[ip], ") ", levels_hab3[ip]),
        "(n =", NROW(iuse_uni), ")")
    } else {
      print(paste(ftemp1, "(n =", NROW(iuse_uni), ")"))
    }
  }

  par(par_prev)
  dev.off()


  #--- Map by responses (and by habitat3)
  dim_fig <- c(if (byHabitat3) length(responses) else 2, 2)
  h.panel <- 2.75; h.edgeL <- 0.3; h.edgeU <- 0.1
  w.panel <- 6; w.edgeL <- 0.3; w.edgeR <- 0.1

  l_widths <- c(w.edgeL, rep(w.panel, times = dim_fig[2]), w.edgeR)
  w.edge <- w.edgeL + w.edgeR
  aspr <- (dim_fig[2] * fws - w.edge) / (sum(l_widths) - w.edge)
  h.panel <- aspr * h.panel
  w.panel <- aspr * w.panel

  pdf(height = h.edgeL + h.panel * dim_fig[1] + h.edgeU,
    width = w.edgeL + w.panel * dim_fig[2] + w.edgeR, file = ftemp2)

  temp <- c(
    rep(0, 1 + dim_fig[2] + 1),
    sapply(seq_len(dim_fig[1]), function(k)
      c(0, (k - 1) * dim_fig[2] + seq_len(dim_fig[2]), 0)),
    rep(0, 1 + dim_fig[2] + 1))
  temp_layout <- matrix(temp, nrow = 1 + dim_fig[1] + 1, ncol = 1 + dim_fig[2] + 1,
    byrow = TRUE)

  layout(temp_layout,
    heights = c(h.edgeU, rep(h.panel, times = dim_fig[1]), h.edgeL),
    widths = c(w.edgeL, rep(w.panel, times = dim_fig[2]), w.edgeR))
  par_prev <- par(mgp = c(1, 0, 0), mar = c(0.5, 0.5, 0, 0), tcl = 0.3, cex = 1)

  ia <- 1

  if (byHabitat3) {
    irs_xlab <- length(responses)
    irs_ylab <- seq_along(responses)
  } else {
    irs_xlab <- c(2, 4)
    irs_ylab <- c(1, 3)
  }

  for (ir in seq_along(responses)) {
    temp_dlocs2 <- dlocs_uni2[dlocs_uni2[, responses[ir]] > 0, ]

    for (ip in seq_along(levels_hab3)) {
      iuse <- temp_dlocs2
      if (byHabitat3) {
        iuse <- iuse[as.character(temp_dlocs2[, "Habitat3"]) == levels_hab3[ip], ]
      }

      with(iuse, plot(Effect_longitude, Effect_latitude, type = "n",
        asp = 1, xlim = xlims, ylim = ylims, ann = FALSE, axes = FALSE))
      if (ir %in% irs_xlab) {
        axis(side = 1, labels = TRUE)
        title(xlab = "Longitude", cex.lab = 1, xpd = NA)
      } else {
        axis(side = 1, labels = FALSE)
      }
      if (ir %in% irs_ylab && ip == 1) {
        axis(side = 2, labels = TRUE)
        title(ylab = "Latitude", cex.lab = 1, xpd = NA)
      } else {
        axis(side = 2, labels = FALSE)
      }
      with(iuse, smoothScatter(Effect_longitude, Effect_latitude, add = TRUE,
        bandwidth = bdw, nbin = c(720, 180), xlim = c(-180, 180), ylim = c(-90, 90),
        nrpoints = 0, colramp = fcols))
      map(add = TRUE)
      with(iuse, points(Effect_longitude, Effect_latitude,
        col = "orange", lwd = 1, cex = 1, pch = 4))
      abline(h = 0, col = "gray")

      #mtext(side = 3, adj = 0.025, font = 2, line = -1.5,
      #  text = paste0("(", letters[ir], ") response = ", responses[ir]))
      temp <- lab_responses[responses[ir] == lab_responses[, "code"], "out"]
      legend("topleft", box.col = "white", bg = "white", text.font = 2, inset = -0.04,
        legend = paste0("(", letters[ia], ") ",
          if (!byHabitat3 || (byHabitat3 && (ip == 1 || ir == 1))) temp,
          if (byHabitat3 && ir == 1) paste0(" / ", levels_hab3[ip], "\n"),
          " (n = ", NROW(iuse), ")"))
      ia <- ia + 1
    }
  }

  par(par_prev)
  dev.off()
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
        dir_res = dir_res_, dir_out = dir_res0
      )

    } else {
      list(
        only_wo_controls = FALSE, withControlsL. = full_design[["s1_wcontr"]],
        only_useadj_standardized = FALSE, withNonStandardizedL. = full_design[["s2_wnonnorm"]],
        fragment_sizes. = full_design[["s4_fragsize"]],
        cor_methods. = full_design[["s5_cormethod"]],
        cor_transforms. = full_design[["s6_cortransform"]],
        weight_methods. = full_design[["d7_weightmethod"]],
        dir_res = dir_res_, dir_out = dir_res0
      )
    }


  template_args <- expand.grid(template_args, stringsAsFactors = FALSE)

  for (k in seq_len(NROW(template_args))) for (h3 in levels_byHabitat3) {
    do.call("msmap_metaanalysis", args = c(template_args[k, ],
      responses = list(responses_Hall),
      deffects = list(deffects),
      dmoderators = list(dmoderators),
      byHabitat3 = list(h3)))
  }
}

