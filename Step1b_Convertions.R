
#--- SETTINGS
# dir_prj <- "~/Dropbox/Work_Share_GenetFragment_MetaAnalysis"
dir_prj <- "Prj04_GenetFragment_MetaAnalysis"

dir_in <- file.path(dir_prj, "2_Data", "4_DataExtracted")
dir_out <- file.path(dir_prj, "2_Data", "5_DataCleaned")
dir_fig <- file.path(dir_out, "Figures")

# File names
if (!exists("tag_newdate")) tag_newdate <- format(Sys.Date(), "%Y%m%d")
if (!exists("tag_usedate")) tag_usedate <- "20180705" #"20180112" # "20171107"

fin <- file.path(dir_out, paste0(tag_usedate, "_ExtractedDataCleaned1_wide.rds"))
foutw <- file.path(dir_out, paste0(tag_newdate, "_ExtractedDataCleaned2_wide.rds"))



#--- READ DATA
print(paste(Sys.time(), "reading input", shQuote(basename(fin))))
dat0 <- readRDS(fin)


#------------------------------------
#--- Compare published values against fragment-level estimates
compare_publed_vs_calc <- function(xpub, xcalc, tag_x, lims = NULL, xsource, dir.out, tol = 1e-2) {
  i_pubs <- xsource %in% "publication"
  i_notidentical <- abs((xpub - xcalc) / xpub) > tol
  i_notidentical <- ifelse(is.na(i_notidentical), FALSE, i_notidentical)
  comp_xpub <- xpub[i_pubs & i_notidentical]
  comp_xcalc <- xcalc[i_pubs & i_notidentical]
  no <- sapply(list(comp_xpub, comp_xcalc, comp_xpub & comp_xcalc, xcalc,
    xpub[i_pubs] & xcalc[i_pubs], xpub[i_pubs]), function(x) sum(is.finite(x)))
  e <- comp_xpub - comp_xcalc
  rmse <- mean(sqrt((e) ^ 2), na.rm = TRUE)
  re <- range(e, na.rm = TRUE)

  if (is.null(lims) || any(!is.finite(lims)))
    lims <- range(c(xpub, xcalc), na.rm = TRUE)

  pdf(height = 10, width = 12, file = file.path(dir.out,
    paste0("Fig_Published_vs_CalculatedAtFragment_", tag_x, ".pdf")))
  op <- par(mfcol = c(2, 2), mar = c(3, 3, 2.5, 0.5), mgp = c(1.5, 0.5, 0))

    plot(comp_xcalc, comp_xpub, xlim = lims, ylim = lims,
      xlab = paste("calculated", tag_x), ylab = paste("published", tag_x))
    title(main = paste0("Fragments with both published and calculated values:\n",
      "n(total; not shown) = ", no[5], "; n(pub != calc) = ", no[3]))
    abline(0, 1, lwd = 2, lty = 2, col = "gray")
    text(lims[1], lims[2], adj = c(0, 1), labels = paste("RMSE =", round(rmse, 3),
      "\nmin error =", round(re[1], 3), "\nmax error =", round(re[2], 3)))

    if (length(comp_xcalc) > 3) {
      xt <- seq(min(comp_xcalc, na.rm = TRUE), max(comp_xcalc, na.rm = TRUE), length = 100)
      yt <- predict(loess(e ~ comp_xcalc), data.frame(comp_xcalc = xt))
      plot(comp_xcalc, e, xlim = lims, ylim = range(c(e, yt), na.rm = TRUE),
        xlab = paste("calculated", tag_x), ylab = "Error between published and calculated")
      lines(xt, yt, col = "red", lwd = 2)
      abline(0, 0, lwd = 2, lty = 2, col = "gray")

    } else {
      plot.new()
    }

    beanplot::beanplot(xpub[i_pubs], xcalc, what = c(0, 1, 1, 0),
      col = list(c(rep("orange", 3), "black"), c(rep("dodgerblue", 3), "black")),
      overallline = "median", ylim = NULL, ylab = tag_x,
      names = paste0(c("published", "calculated"), ": n = ", no[c(6, 4)]))
    abline(0, 0, lwd = 2, lty = 2, col = "gray")

  par(op)
  dev.off()
}

identify_largest_diffs <- function(id, xcalc, xpub, xsource, method = "relative") {
  i_pubs <- xsource %in% "publication"
  e_abs <- abs(xpub - xcalc)
  e_rel <- abs(e_abs / xpub)

  dtemp <- data.frame(id = id, abs_error = e_abs, rel_error = e_rel,
    xcalc = xcalc, xpub = xpub)
  dtemp <- dtemp[complete.cases(dtemp), ]

  if (identical(method, "relative")) {
    dtemp[order(dtemp$rel_error, decreasing = TRUE), ]
  } else {
    dtemp[order(dtemp$abs_error, decreasing = TRUE), ]
  }
}



temp_ids <- unname(apply(dat0[, c("ID_effect", "ID_fragment")], 1, paste, collapse = "_"))


#--- CONVERT FRAGMENT SIZE METRICS
# fragment area (m2) <- N of individuals / density
ina <- is.na(dat0$Fragment_area_ha)
temp <- with(dat0, Pop_size_n / Pop_density_nPERha)
dat0$Fragment_area_ha[ina] <- temp[ina]

dat0$Fragment_area_source <- rep("publication", dim(dat0)[1])
dat0$Fragment_area_source[ina] <- ifelse(is.na(dat0$Fragment_area_ha[ina]), NA,
  "calc_from_density")
print(table(dat0$Fragment_area_source))
#calc_from_density       publication
#               12              1169

diffs_area <- identify_largest_diffs(id = temp_ids, xcalc = temp,
  xpub = dat0[, "Fragment_area_ha"], xsource = dat0[, "Fragment_area_source"])
diffs_area[diffs_area$rel_error > 0.8, ]
#                                                           id  abs_error rel_error      xcalc  xpub
#157        7180_Microseris_lanceolata-Isozymes_Barton_art7180  0.4500000 0.9000000 0.05000000  0.50
#641   7175_Rutidosis_leptorrhynchoides-Isozymes_Manor_art7175  0.9000000 0.9000000 0.10000000  1.00
diffs_area[diffs_area$abs_error > 10, ]

compare_publed_vs_calc(xpub = dat0[, "Fragment_area_ha"], xcalc = temp, tag_x = "Fragment area (ha)",
  lims = c(0, 50), xsource = dat0[, "Fragment_area_source"], dir.out = dir_fig)


# N of individuals <- fragment area (m2) * density
ina <- is.na(dat0$Pop_size_n)
temp <- with(dat0, Fragment_area_ha * Pop_density_nPERha)
dat0$Pop_size_n[ina] <- temp[ina]

dat0$Pop_size_source <- rep("publication", dim(dat0)[1])
dat0$Pop_size_source[ina] <- ifelse(is.na(dat0$Pop_size_n[ina]), NA,
  "calc_from_density")
print(table(dat0$Pop_size_source))
#calc_from_density       publication
#               38               599

diffs_popsize <- identify_largest_diffs(id = temp_ids, xcalc = temp,
  xpub = dat0[, "Pop_size_n"], xsource = dat0[, "Pop_size_source"])
diffs_popsize[diffs_popsize$rel_error > 0.8, ]
#                                                           id abs_error   rel_error  xcalc   xpub
#157        7180_Microseris_lanceolata-Isozymes_Barton_art7180     15300   9.0000000  17000   1700
#641   7175_Rutidosis_leptorrhynchoides-Isozymes_Manor_art7175       117   9.0000000    130     13
diffs_popsize[diffs_popsize$abs_error > 1e3, ]

compare_publed_vs_calc(xpub = dat0[, "Pop_size_n"], xcalc = temp, tag_x = "Pop size",
  xsource = dat0[, "Pop_size_source"], dir.out = dir_fig)
compare_publed_vs_calc(xpub = dat0[, "Pop_size_n"], xcalc = temp, tag_x = "Pop size",
  lims = c(0, 100000), xsource = dat0[, "Pop_size_source"], dir.out = dir_fig)



#--- CONVERT HE AND H0 TO FIS
# Fis = (He-Ho)/He
# Hans-Peter (Dec 22, 2016): Mit den He und Ho Werten auf Fragment-Ebene ergibt diese
#   Umrechnung eine sehr grobe Fis-Berechung. Korrekt sollte diese Berechnung mit jedem
#   Primer durchgeführt werden
ina <- is.na(dat0$Response_Fis)
temp <- with(dat0, (Response_He - Response_H0) / Response_He)
dat0$Response_Fis[ina] <- temp[ina]

dat0$Response_Fis_source <- rep("publication", dim(dat0)[1])
dat0$Response_Fis_source[ina] <- ifelse(is.na(dat0$Response_Fis[ina]), NA,
  "calc_from_fragment-level_HeH0")
print(table(dat0$Response_Fis_source))
#calc_from_fragment-level_HeH0                   publication
#                          426                           521

diffs_fis <- identify_largest_diffs(id = temp_ids, method = "absolute", xcalc = temp,
  xpub = dat0[, "Response_Fis"], xsource = dat0[, "Response_Fis_source"])
diffs_fis[diffs_fis$abs_error > 0.1, ]
#                                                       id abs_error  rel_error       xcalc   xpub
#570            7001_Manilkara_maxima-Microsat_LM6_art7001 0.3433333 34.3333333 -0.33333333  0.010
#569            7001_Manilkara_maxima-Microsat_LM4_art7001 0.3300000  4.1250000 -0.25000000  0.080
#1340              7125_Ocotea_porosa-Microsat_Op6_art7125 0.2029825 20.2982456  0.19298246 -0.010
#480                 6904_Quercus_robur-Isozymes_1_art6904 0.1671961  2.2903572  0.24019608  0.073
#1635 7955_Alouatta_palliata_mexicana-Microsat_GF3_art7955 0.1365217  0.7584541  0.04347826  0.180
#1314          7017_Petaurus_breviceps-Microsat_MM_art7017 0.1331475 66.5737705  0.13114754 -0.002
#568            7001_Manilkara_maxima-Microsat_LM3_art7001 0.1331250  0.3095930  0.29687500  0.430
#1307          7017_Petaurus_breviceps-Microsat_CR_art7017 0.1098571  3.3290043  0.14285714  0.033
#571            7001_Manilkara_maxima-Microsat_LM7_art7001 0.1082353  0.2081448  0.41176471  0.520

compare_publed_vs_calc(xpub = dat0[, "Response_Fis"], xcalc = temp, tag_x = "Fis",
  lims = c(-1, 1), xsource = dat0[, "Response_Fis_source"], dir.out = dir_fig)
#------------------------------------


#--- CALCULATE RANGE OF FRAGMENT SIZE PER EFFECT
range_per_effect <- function(xname, data, with_control = FALSE) {

  fdat <- if (!with_control) data[data[, "is_fragment"], ] else data
  temp <- tapply(fdat[, xname], fdat$ID_effect, function(x) diff(range(x, na.rm = TRUE)))
  temp[!is.finite(temp)] <- NA

  res2 <- rep(NA, nrow(data))
  itemp <- match(data[, "ID_effect"], names(temp), nomatch = 0)
  res2[itemp > 0] <- temp[itemp]
#  if (!with_control) {
#    res2[data[, "is_fragment"]] <- NA
#  }

  res2
}

dat0$Fragment_area_ha_range <- range_per_effect("Fragment_area_ha", dat0, FALSE)
dat0$Pop_size_n_range <- range_per_effect("Pop_size_n", dat0, FALSE)


#------------------------------------
#--- SAVE CLEANED DATA
saveRDS(dat0, file = foutw)

