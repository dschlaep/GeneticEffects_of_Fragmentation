library("gsl")

#--- OPTIONS
# Methods for variance, variance of unbiased coefficient, and variance of z-transform
var_options <- data.frame(
  #--- PEARSON
  # * r options: "Asymptotic1", "Allaire2005"
  # * unbiased-r options: "Hedges1989"
  # * ztransform-options: "Fisher1921"
  pearson = c(none = "Allaire2005", unbiased = "Hedges1989", ztransform = "Fisher1921"),

  #--- KENDALL
  # * r options: "Kendall1938", "Esscher1924"
  # * ztransfrom-options: "FiellerHartleyPearson1957"
  # Puth et al. 2015: "we can recommend the Fisher transformation approach
  #   combined with Method C [Xu et al. 2013 citing Esscher 1924] for variance
  #   estimation as an effective way to calculate confidence intervals for
  #   Kendall's tau."
  # Note: "Kendall's tb generally provides values less than 0.95 especially for large
  #   correlation values and a high percentage of ties (n = 50)"
  kendall = c(none = "Esscher1924", unbiased = "", ztransform = "FiellerHartleyPearson1957"),

  #--- SPEARMAN
  # * r options: "KendallStuart1973", "DavidKendallStuart1951", "FiellerHartleyPearson1957",
  #     "XuHouHungZou2013" (not implemented), "Borkowf2000" (not implemented)
  # * ztransfrom-options: "FiellerHartleyPearson1957", "BonettWright2000", "CarusoCliff1997"
  # Puth et al. 2015: "Method C [Caruso & Cliff 1997 for spearman's r->z] can perhaps be
  #   recommended [if there are no ties], since it offers the most consistently good
  #   performance over all the scenarios we explored."
  # Puth et al. 2015: "As with our simulations without ties, variance estimation B
  #   [Bonett & Wrigth 2000] generally provides values closest to the nominal 0.95 for
  #   different combinations of r and n."
  spearman = c(none = "FiellerHartleyPearson1957", unbiased = "", ztransform = "BonettWright2000"),

  stringsAsFactors = FALSE
)


#--- SETTINGS
dir_prj <- "Prj04_GenetFragment_MetaAnalysis"
dir_in <- file.path(dir_prj, "2_Data", "4_DataExtracted")
dir_out <- file.path(dir_prj, "2_Data", "5_DataCleaned")
dir_fig1 <- file.path(dir_out, "Figures")

# File names
if (!exists("tag_newdate")) tag_newdate <- format(Sys.Date(), "%Y%m%d")
if (!exists("tag_usedate")) tag_usedate <- "20180705" #"20180112"

fin <- file.path(dir_out, paste0(tag_usedate, "_ExtractedDataCleaned3_long.rds"))
fout_esr <- file.path(dir_out, paste0(tag_usedate, "_EffectSize_Correlation.rda"))
fout_esr_new <- file.path(dir_out, paste0(tag_newdate, "_EffectSize_Correlation.rda"))
fout_mod <- file.path(dir_out, paste0(tag_usedate, "_ModeratorPredictors.rds"))
fout_mod_new <- file.path(dir_out, paste0(tag_newdate, "_ModeratorPredictors.rds"))


#--- READ DATA
print(paste(Sys.time(), "reading input", shQuote(basename(fin))))
datl <- readRDS(fin)


#------------------------------------
#--- CALCULATE EFFECT SIZE

#--- ID of effect
# One effect size for each article x species x marker
ID_effect <- unique(datl$ID_effect)

is_fragment <- ifelse(is.na(datl$is_fragment), FALSE, datl$is_fragment)
is_standardized <- ifelse(is.na(datl$is_standardized), FALSE, datl$is_standardized)
is_published_farea <- ifelse(is.na(datl$Fragment_area_source), FALSE,
  datl$Fragment_area_source == "publication")
is_published_psize <- ifelse(is.na(datl$Pop_size_source), FALSE,
  datl$Pop_size_source == "publication")

# Loop over response (genetic measure), fragment area, correlation method,
#   include/exclude controls, include/exclude non-standardized
iuse0 <- rep(TRUE, dim(datl)[1])

# Sensitivity factors
sf <- list()
sf[["s1_wcontr"]] <- c(TRUE, FALSE)
sf[["s2_wnonnorm"]] <- c(TRUE, FALSE)
sf[["s3_responses"]] <- c(sort(unique(datl$Response)), "He+Sh")
sf[["s4_fragsize"]] <- c("Fragment_area_ha-published", "Fragment_area_ha-pub&calc",
  "Fragment_area_rank-published", "Pop_size_n-published", "Pop_size_n-pub&calc",
  "Pop_size_rank-published")
sf[["s5_cormethod"]] <- c("pearson", "kendall", "spearman")
sf[["s6_cortransform"]] <- c("none", "unbiased", "ztransform")

# Output
#   nf = number of fragments
#   ni = mean number of individuals per fragment
#   E = estimated correlation coefficient
#   var = variance of coefficient
#   winvar = weight as 1 / var
#   whierarch = weight based on hierarchical approach by Reed & Frankham 2003
sf[["out_var"]] <- c("nf", "ni", "E", "var", "winvar", "whierarch")

sf_dims <- lengths(sf)




#------ Data per effect
if (!file.exists(fout_mod)) {

  es_dat <- unique(datl[, c(
    # Identification
    "ID_effect", "ID_article", "Paper_Nr", "ID_unit", "FirstAuthor",
    # Study project
    "Pub_Year", "Study_type", "Marker", "Marker_cat",
    "Effect_longitude", "Effect_latitude", "Effect_latitude_bin",
    "EffectAgeCat2_ord_yrs", "Effect_range_age_yr",
    "Fragment_area_ha_range", "Pop_size_n_range",
    # Study subject
    "Habitat2", "Habitat3", "Matrix2", "Organism_group4", "Organism_group5",
    "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species_resolved",
    # Plant traits
    "Plants_Lifeform", "Plants_Sex", "Plants_Self_compatibility", "Plants_Mating_system",
    "Plants_Pollination_syndrome", "Plants_Pollination_vector",
    "Plants_Seed_dispersal_vector",
    # Animal traits
    "Animals_Length_cat", "Animals_Mass_cat", "Animals_Generation_cat",
    "Animals_Dispersal_ability", "Animals_Trophic_group",
    # Shared traits
    "SpeciesAgeCat2_ord_yrs",
    # Time
    "Effect_NGenSinceFrag_ord", "Effect_NlifeSinceFrag_ord"
    )])

  if (anyDuplicated(es_dat$ID_effect)) {
    ibad <- duplicated(es_dat$ID_effect) | duplicated(es_dat$ID_effect, from.last = TRUE)
    ibad_ids <- unique(es_dat$ID_effect[ibad])
    print("Effects are not uniquely described:")

    for (k in ibad_ids) {
      temp <- es_dat[es_dat$ID_effect %in% k, ]
      print(temp[do.call(order, temp), ])
    }

    stop("Describe effects uniquely before calculating effect sizes")
  }

  # Convert any character vectors to factors
  for (k in seq_len(dim(es_dat)[2])) {
    if (is.character(es_dat[, k]) && !is.factor(es_dat[, k])) {
      es_dat[, k] <- as.factor(es_dat[, k])
    }
  }

  # Save data
  saveRDS(es_dat, file = fout_mod_new)

} else {
  es_dat <- readRDS(fout_mod)
}


# Number of effects
N_effects <- length(es_dat$ID_effect)
print(N_effects) # 155

# Number of fragments for continuous versus ranked fragment size measures
temp <- datl[, c("Fragment_area_ha", "Fragment_area_rank", "Fragment_area_source")]
temp[temp$Fragment_area_source %in% "calc_from_density", "Fragment_area_ha"] <- NA
temp <- temp[, c("Fragment_area_ha", "Fragment_area_rank")]
temp <- aggregate(temp, list(datl$ID_effect), function(x) sum(!is.na(x), na.rm = TRUE))

print(temp[which(temp[, "Fragment_area_ha"] != temp[, "Fragment_area_rank"]), ])
print(temp[which(temp[, "Fragment_area_ha"] == 0 & temp[, "Fragment_area_rank"] > 0), ])
# <0 rows> (or 0-length row.names) ==> no new effect

print(stemp <- apply(temp[, -1], 2, sum))
  # Fragment_area_ha Fragment_area_rank
  #            10521              10989
diff(stemp) # 468


#------ Calculate effect sizes ES
do_calc_effects <- TRUE
if (file.exists(fout_esr)) {
  opt_methods_cur <- var_options
  # Load data: es_cor, var_options
  load(fout_esr)

  do_calc_effects <- !identical(var_options, opt_methods_cur)
}

if (do_calc_effects) {
  #--- CORRELATION COEFFICIENT
  #' @references Olkin, I. & Pratt, J.W. (1958). Unbiased estimation of certain
  #'    correlation coefficients. Ann Math Stat, 29, 201-211.
  #' @references Fisher, R.A. (1921). On the probable error of a coefficient of
  #'    correlation deduced from a small sample. Metron, 1, 3-32.
  #' @references Allaire, D. & Laurencelle, L. (2005). The Sampling Variance of r(X, Y)
  #'    for Unrelated X and Y: A Simple Algebraic Demonstration. Communications in
  #'    Statistics - Theory and Methods, 33, 2389-2392.
  #' @references Hunter, J.E. & Schmidt, F.L. (2015). Methods of meta-analysis:
  #'    correcting error and bias in research findings. Third edition. edn. SAGE,
  #'    Thousand Oaks, California.
  #' @references Reed, D.H. & Frankham, R. (2003). Correlation between fitness and
  #'    genetic diversity. Conservation Biology, 17, 230-237.
  #' @references Reed, D.H. & Frankham, R. (2001). How closely correlated are molecular
  #'    and quantitative measures of genetic variation? A meta-analysis. Evolution, 55,
  #'    1095-1103.
  #' @references Fieller, E. C., H. O. Hartley, and E. S. Pearson. (1957). Tests for Rank
  #'    Correlation Coefficients. I. Biometrika 44 (3–4): 470–81.
  #' @references Kendall, M. G. (1938). A New Measure of Rank Correlation. Biometrika
  #'    30 (1–2): 81–93.
  #' @references Bonett, Douglas G., and Thomas A. Wright. (2000). Sample Size
  #'    Requirements for Estimating Pearson, Kendall and Spearman Correlations.
  #'    Psychometrika 65 (1): 23–28.
  #' @references Puth, Marie-Therese, Markus Neuhäuser, and Graeme D. Ruxton. (2015).
  #'    Effective Use of Spearman’s and Kendall’s Correlation Coefficients for Association
  #'    between Two Measured Traits. Animal Behaviour 102 (April): 77–84.
  #' @references Xu, Weichao, Yunhe Hou, Y. S. Hung, and Yuexian Zou. (2013). A
  #'    Comparative Analysis of Spearman’s Rho and Kendall’s Tau in Normal and
  #'    Contaminated Normal Models. Signal Processing 93 (1): 261–76.
  #' @references Esscher, Fredrick. (1924). On a Method of Determining Correlation from
  #'    the Ranks of the Variates. Scandinavian Actuarial Journal 1924 (1): 201–19.
  #' @references Caruso, John C., and Norman Cliff. 1997. Empirical Size, Coverage, and
  #'    Power of Confidence Intervals for Spearman’s Rho. Educational and Psychological
  #'    Measurement 57 (4): 637–54.

  calc_effectsize_cor <- compiler::cmpfun(function(x, method, tf, var_method,
    verbose = FALSE) {
    x <- x[complete.cases(x), ]

    #--- SAMPLE SIZES
    # number of fragments
    n <- dim(x)[1]
    n <- if (is.finite(n) && n > 0) n else NA

    # number of individuals per fragment used to determine value of genetic measure
    ni <- if (is.finite(n)) mean(x[, 3]) else NA


    #--- EFFECT SIZE ESTIMATES
    # correlation coefficient (metafor::escalc(method = "COR"))
    # we use -cor instead of +cor: because we consider fragmentation effect,
    # i.e., rev(fragment size), and -cor(x, y) = cor(-x, y)
    r <- if (is.finite(n)) {
        -cor(x[, 1], x[, 2], method = method, use = "na.or.complete")
      } else NA

    rp <- if (method == "pearson") {
        r
      } else {
        -cor(x[, 1], x[, 2], method = "pearson", use = "na.or.complete")
      }

    # Fisher's r-to-z transform
    #   - Fisher, 1921; as used by metafor::escalc(method = "ZCOR") for pearson
    #   - Fieller et al. 1957: eq. 9 for kendall and spearman
    z <- if (is.finite(n) && is.finite(r)) {
        1 / 2 * log((1 + r) / (1 - r))
      } else NA

    # Olkin & Pratt's unbiased correlation coefficient (Olkin & Pratt, 1958)
    #  (metafor::escalc(method = "UCOR"))
    ru <- if (is.finite(n) && is.finite(r) && method == "pearson") {
        # bias corrected (based on equation 2.3 in Olkin & Pratt, 1958):
        #     r(ucor) = r * F_hypergeometric(1/2, 1/2, ([N - 1] - 1) / 2, 1 - r^2)
        g <- (n - 2) / 2
        if (g > 1) {
          r * gsl::hyperg_2F1(1 / 2, 1 / 2, g, 1 - r^2)
        } else NA
      } else NA


    # Estimate of correlation, unbiased correlation, or z-transform
    est <- switch(tf, none = r, unbiased = ru, ztransform = z, NA)


    #--- VARIANCE OF CORRELATION ESTIMATE
    evar <- NA

    if (is.finite(n)) {
      if (tf == "none") {
        #--- within-study sampling variance of the correlation coefficient
        # Variance of r (assuming iid and bivariate normal)
        if (method == "pearson") {
          if (var_method == "Asymptotic1" && n > 1) {
            # Large-sample approximation (metafor::escalc(vtype = "LS"))
            #   - Allaire et al. 2005: "usual asymptotic estimate, which is
            #       (1 - rho^2)^2 / n, is almost useless for small n and abs(rho) >> 0,
            #     and it is deservedly scorned"
            #   - Hunter & Schmidt 2004 (p. 197): fixed-effect case (cited by Brannick et al. 2010, eq. 1)
            evar <- (1 - r ^ 2) ^ 2 / (n - 1)

          } else if (var_method == "Allaire2005" && n > 2) {
            # Allaire et al. 2005
            evar <- (1 - r ^ 2) ^ (2 * (n - 2) / n) / (n - 1)
          }

        } else if (method == "spearman") {
          if (var_method == "KendallStuart1973" && n > 2 && is.finite(r)) {
            # Kendall & Stuart 1973
            # Hüsler (2015) script: approximation good for n >= 10 and for n between 10 and
            #   20 better than z-transform
            evar <- (1 - r^2) / (n - 2)

          } else if (var_method == "DavidKendallStuart1951" && n > 0) {
            # "Large sample approximation due to Kendall 1949 and David et al. 1951. ...
            # it does not appear to very very accurate when sample size is as small as 10.
            # ... The Kendall formula (6) does not give the correct value of 1/(n-1) and 0
            # to var(rs) when rho = 0 and 1, respectively." (Fieller et al. 1957): eq. 6
            evar <- 1 / n * (1 - 1.563465 * rp^2 + 0.304743 * rp^4 + 0.155286 * rp^6 +
              0.061552 * rp^8 + 0.022099 * rp^10)

          } else if (var_method == "FiellerHartleyPearson1957" && is.finite(rp) && n > 1) {
            # "A purely empirical adjustment is obtained [from Kendall formula (6)] by
            # substituting n-1 for n as divisor and adding a term +0.019758*rho^12 which
            # reduces the variance to zero when rho = 1" (Fieller et al. 1957): eq. 11
            evar <- 1 / (n - 1) * (1 - 1.563465 * rp^2 + 0.304743 * rp^4 + 0.155286 * rp^6 +
              0.061552 * rp^8 + 0.022099 * rp^10 + 0.019758 * rp^12)

          } else if (var_method == "XuHouHungZou2013") {
            stop("Variance method 'XuHouHungZou2013' (Lemma 3) is not implemented")

          } else if (var_method == "Borkowf2000") {
            stop("Variance method 'Borkowf2000' is not implemented")
          }

        } else if (method == "kendall") {
          if (var_method == "Kendall1938" && n > 2) {
            # Normal approximation good for n > 10; see eq. 6 in Kendall 1938 Biometrika
            # Also used by SuppDists::sKendall
            evar <- 2 * (2 * n + 5) / (9 * n * (n - 1))

          } else if (var_method == "Esscher1924" && is.finite(rp) && n > 2) {
            # Esscher 1924
            pi2 <- pi^2
            evar <- 2 / (n * (n - 1)) * (1 - 4 / pi2 * asin(rp)^2 + 2 * (n - 2) *
              (1 / 9 - 4 / pi2 * asin(rp / 2)^2))
          }
        }

      } else if (tf == "unbiased") {
        # Variance of ru (assuming iid and bivariate normal)
        if (method == "pearson") {
          if (var_method == "Hedges1989" && n > 3 && is.finite(ru)) {
            # "Hedges, 1989 (see eq. 18), but using the exact equation for Q instead of the
            # approximation (see eq. A10)"; as implemented by metafor::escalc(vtype = "UB")
            f <- gsl::hyperg_2F1(1, 1, n / 2, 1 - r^2)
            Q <- 1 - (n - 3) * (1 - r^2) * f / (n - 2)
            evar <- ru^2 - Q
          }
        }

      } else if (tf == "ztransform") {
        # Variance of z-scale r  (assuming iid and bivariate normal)
        if (method == "pearson") {
          if (var_method == "Fisher1921" && n > 3) {
            # Fisher, 1921; as used by metafor::escalc(method = "ZCOR")
            evar <- 1 / (n - 3)
          }

        } else if (method == "spearman") {
          if (var_method == "FiellerHartleyPearson1957" && n > 3) {
            # Fieller, Hartley, & Pearson 1957: eq. 10; only good for abs(r) < 0.8 && n >= 10
            evar <- 1.060 / (n - 3)

          } else if (var_method == "BonettWright2000" && is.finite(r) && n > 3) {
            # Bonett & Wright 2000
            evar <- (1 + r^2 / 2) / (n - 3)

          } else if (var_method == "CarusoCliff1997" && is.finite(z) && n > 2) {
            # Caruso & Cliff 1997
            evar <- 1 / (n - 2) + abs(z) / (6 * n + 4 * sqrt(n))
          }

        } else if (method == "kendall") {
          if (var_method == "FiellerHartleyPearson1957" && n > 4) {
            # Fieller, Hartley, & Pearson 1957: eq. 10; only good for abs(r) < 0.8 && n >= 10
            evar <- 0.437 / (n - 4)
          }
        }
      }

      if (!is.finite(evar) || evar < 0) {
        if (verbose)
          print(paste("'variance estimate' for", method, tf, var_method,
            "was bad with var =", evar, "code set to NA"))
        evar <- NA
      }
    }

    #--- STUDY WEIGHTS
    # Inverse of variance
    winvar = 1 / evar

    # # fragments and # individuals per fragment known: Reed & Frankham 2003
    #   (cited by Vranckx et al. 2012)
    whierarch <- if (is.finite(n) && n > 2 && is.finite(ni)) {
        sqrt((n - 2) * ni)
      } else NA

    c(n, ni, est, evar, winvar, whierarch)
  })


  # Output container
  es_cor <- array(NA, dim = c(N_effects, sf_dims),
    dimnames = c(list(ID_effect = es_dat$ID_effect), sf))

  # Loop over all combinations of data to calculate effect sizes
  for (k1 in seq_len(sf_dims["s1_wcontr"])) {
    iuse1 <- if (sf[["s1_wcontr"]][k1]) iuse0 else iuse0 & is_fragment

    for (k2 in seq_len(sf_dims["s2_wnonnorm"])) {
      iuse2 <- if (sf[["s2_wnonnorm"]][k2]) iuse1 else iuse1 & is_standardized

      for (k3 in seq_len(sf_dims["s3_responses"])) {
        iuse3 <- iuse2 & datl$Response == sf[["s3_responses"]][k3]

        if (any(iuse3)) for (k4 in seq_len(sf_dims["s4_fragsize"])) {
          temp <- strsplit(sf[["s4_fragsize"]][k4], split = "-", fixed = TRUE)[[1]]
          tag_fragsize <- temp[1]
          iuse4 <- if (temp[2] == "published") {
              temp1 <- if (grepl("Fragment_area_", tag_fragsize)) {
                  is_published_farea
                } else if (grepl("Pop_size_", tag_fragsize)) {
                  is_published_psize
                } else NULL
              iuse3 & temp1
            } else iuse3

          if (anyNA(iuse4)) stop("'iuse4' contains NAs")
          if (sum(iuse4) == 0) next

          xdat <- datl[iuse4, c(tag_fragsize, "Value", "Ind_sampled_nPERfragment")]

          for (k5 in seq_len(sf_dims["s5_cormethod"])) {
            for (k6 in seq_len(sf_dims[["s6_cortransform"]])) {

              print(paste(if (sf[[1]][k1]) "withControls" else "withoutControls",
                if (sf[[2]][k2]) "withNonStandardized" else "onlyStandardized",
                sf[[3]][k3], sf[[4]][k4], sf[[5]][k5], sf[[6]][k6], sep = " - "))

              val <- by(xdat, datl$ID_effect[iuse4], calc_effectsize_cor,
                method = sf[["s5_cormethod"]][k5], tf = sf[["s6_cortransform"]][k6],
                var_method = var_options[sf[["s6_cortransform"]][k6], sf[["s5_cormethod"]][k5]])

              temp <- do.call(rbind, val)
              stopifnot(dim(temp)[2] == sf_dims["out_var"])
              temp1 <- array(unlist(temp), dim = dim(temp))
              es_cor[dimnames(temp)[[1]], k1, k2, k3, k4, k5, k6, ] <- temp1
            }
          }
        }
      }
    }
  }

  # Combine effect sizes for He and Sh based on Aguilar et al. 2008: "In cases where
  # heterozygosity was not given (typically in studies using random amplified polymorphic
  # DNA or amplified fragment length polymorhphism), we used molecular variance or gene
  # diversity and analysed these parameters together with expected heterozygosity."
  #' @references Aguilar, R., Quesada, M., Ashworth, L., Herrerias-Diego, Y. & Lobo, J.
  #'   (2008). Genetic consequences of habitat fragmentation in plant populations:
  #'   susceptible signals in plant traits and methodological approaches. Molecular
  #'   Ecology, 17, 5177-88.
  temp <- es_cor[, , , "He", , , , ]
  ids <- is.na(temp)
  temp[ids] <- es_cor[, , , "Sh", , , , ][ids]
  es_cor[, , , "He+Sh", , , , ] <- temp

  # Save data
  save(es_cor, var_options, file = fout_esr_new)
}

