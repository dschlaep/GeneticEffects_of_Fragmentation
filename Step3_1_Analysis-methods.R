
#' Calculate the mean for numeric and the mode for categorical variables
central <- function(x) {
  if (is.factor(x) || is.character(x)) {
    names(which.max(table(x)))
  } else {
    mean(x, na.rm = TRUE)
  }
}

# Calculate the confidence interval of a rma.uni or rma.mv fit
calc_confint <- function(fit) {
  fit[["confint"]] <- if (isTRUE(fit[["has_fit"]])) {
      temp <- try(confint(fit[["fit"]]))
      if (!inherits(temp, "try-error")) temp else NULL
    } else NULL

  fit
}


get_categories <- function(data, moderators) {
  x <- list()
  if (!is.null(data)) {
    x[["names"]] <- moderators
    x[["is_redundant"]] <- if (length(moderators) > 1) {
        !sapply(moderators, function(x)
          isTRUE(var(as.numeric(data[, x]), na.rm = TRUE) > 0))
      } else FALSE
    x[["is_continuous"]] <- !sapply(moderators, function(x)
      isTRUE(is.character(data[, x]) || is.factor(data[, x])))

    x[["cats"]] <- lapply(moderators, function(x) as.character(sort(unique(data[, x]))))
    x[["cats"]][x[["is_continuous"]]] <- NA
    names(x[["cats"]]) <- moderators
    x[["n_cats"]] <- lengths(x[["cats"]])
  }
  x
}

# Percentage of total amount of heterogeneity can be accounted for by including moderators
calc_explained_heterogeneity <- function(fit_nomod, fit_mod) {
  if (fit_nomod[["has_fit"]] && fit_mod[["has_fit"]]) {
    h <- "heterogeneity_total"
    (fit_nomod[["fit"]][[h]] - fit_mod[["fit"]][[h]]) / fit_nomod[["fit"]][[h]]
  } else NULL
}

calc_overall_I2 <- function(m) {
  if (inherits(m, "rma.uni") && !is.null(m[["I2"]]) && !is.na(m[["I2"]])) {
    m[["heterogeneity_total"]] <- m[["tau2"]]
    m[["I2_overall"]] <- m[["I2"]]

  } else {
    # Calculate tau^2 (http://www.metafor-project.org/doku.php/analyses:konstantopoulos2011):
    # "the sum of the variance components can be interpreted as the total amount of
    # heterogeneity in the true effects"
    m[["heterogeneity_total"]] <- sum(m[["sigma2"]], na.rm = TRUE) +
      sum(m[["tau2"]], na.rm = TRUE) + sum(m[["gamma2"]], na.rm = TRUE)

    # Calculate generalized I^2 (Nakagawa & Santos, 2012;
    # http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate):
    # "statistic can be thought of as the overall I2 value that indicates how much of the
    # total variance can be attributed to the total amount of heterogeneity (which is
    # the sum of between- and within-cluster heterogeneity)"
    W <- solve(m[["V"]])
    X <- model.matrix(m)
    P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    m[["sampling_variance_typical"]] <- (m[["k"]] - m[["p"]]) / sum(diag(P))

    m[["I2_overall"]] <- 100 * m[["heterogeneity_total"]] /
      (m[["heterogeneity_total"]] + m[["sampling_variance_typical"]])
  }

  m
}


#' Create short name for interaction effects
name_interaction <- function(mods) {
  paste(substr(mods, 1, 1), collapse = "x")
}


#' Calculate data.frame where columns are aggregated by species
agg_by_sp <- function(d, tree) {
  ids <- d[, "Species"] %in% tree[["tip.label"]]

  if (any(ids)) {
    temp <- d[ids, , drop = FALSE]
    temp2 <- aggregate(temp, by = list(temp[, "Species"]), central)
    temp2[, "E_se"] <- tapply(temp[, "E"], temp[, "Species"], sd)
    # Sort data according to species in subset trees
    temp_ids <- match(tree[["tip.label"]], temp2[, "Species"], nomatch = 0)
    temp2[temp_ids, -1]
  } else {
    NULL
  }
}

#' Estimate phylogenetic signal: Blomberg's K
#--- Estimated Blomberg's K (Blomberg et al., 2003), which measures the strength of the
#phylogenetic signal in phylogenetic trees, in order to understand the potential
#contrasting results between traditional and phy- logenetically independent
#meta-analyses (Chamberlain et al., 2012). In the context of a meta-analysis, K values
#approaching zero imply that closely related species do not share similar effect sizes,
#whereas values of K near or larger than one suggest that closely related species do
#share similar effect sizes (i.e., effect sizes are conserved). This parameter was
#obtained using the R package “phytools” and the function phylosig (Revell, 2012).
calc_BlombergsK <- function(data, yin = "E", yinse = "E_se", tree) {
  temp_mean <- data[, yin]
  names(temp_mean) <- data[, "Species"]

  temp_se <- data[, yinse]
  temp_se[is.na(temp_se)] <- 0
  names(temp_se) <- data[, "Species"]
  if (all(temp_se < sqrt(.Machine$double.eps))) {
    temp_se <- NULL
  }

  phytools::phylosig(tree = tree, x = temp_mean, se = temp_se, method = "K",
    test = TRUE, nsim = 2000)
}

#' Influence measures
#' "calculates the following leave-one-out diagnostics for each study:
#' externally standardized residual, DFFITS value, Cook's distance, covariance ratio, the
#' leave-one-out amount of (residual) heterogeneity, the leave-one-out test statistic for
#' the test of (residual) heterogeneity, DFBETAS value(s)
calc_influence <- function(m) {
  out <- list()

  if (m[["has_fit"]] && inherits(m[["fit"]], "rma.uni")) {
    temp <- try(influence(m[["fit"]]))

    if (!inherits(temp, "try-error")) {
      out[["infl"]] <- list(temp)
      out[["N_infl"]] <- sum(temp[["is.infl"]], na.rm = TRUE)
    }
  }

  out
}

#' Calculate Fail-safe numbers
calc_failsafeN <- function(d, yin, vin) {
  out <- list()

  if (nrow(d) > 0) {
    # The Rosenthal method (sometimes called a ‘file drawer analysis’) calculates the number
    # of studies averaging null results that would have to be added to the given set of
    # observed outcomes to reduce the combined significance level (p-value) to a target
    # alpha level (e.g., .05). The calculation is based on Stouffer's method to combine
    # p-values and is described in Rosenthal (1979).
    out[["Rosenthal"]] <- metafor::fsn(yi = d[, yin], vi = d[, vin], type = "Rosenthal")

    # The Orwin method calculates the number of studies averaging null results that would
    # have to be added to the given set of observed outcomes to reduce the (unweighted)
    # average effect size to a target (unweighted) average effect size. The method is
    # described in Orwin (1983).
    out[["Orwin"]] <- metafor::fsn(yi = d[, yin], vi = d[, vin], type = "Orwin")

    # The Rosenberg method calculates the number of studies averaging null results that
    # would have to be added to the given set of observed outcomes to reduce significance
    # level (p-value) of the (weighted) average effect size (based on a fixed-effects model)
    # to a target alpha level (e.g., .05). The method is described in Rosenberg (2005).
    out[["Rosenberg"]] <- metafor::fsn(yi = d[, yin], vi = d[, vin], type = "Rosenberg")
  }

  out
}



# https://github.com/nicebread/meta-showdown

# threePSM.est(d=MAdat$d, v=MAdat$v, min.pvalues=1, long=TRUE) from
# fiel 'meta-showdown/MA-methods/7-Selection Models.R' on commit
# https://github.com/nicebread/meta-showdown/commit/dcd50fd41f229c9b10ef3c93d3b65868f4071888



# Estimate the 3PSM (a single step-cutpoint) with the weightr package (should be equivalent to the McShane implementation)
# the single cut point is at .025 (one-sided testing)
# The authors suggest to have >= 4 p-values in each interval. If that is not provided, return NA.

#' @param d A vector of meta-analytic effect sizes.
#' @param v A vector of meta-analytic sampling variances
#' @param min.pvalues How many p-values must be present in each bin that the function returns an estimate?

threePSM.est <- function(d, v, min.pvalues=1, long=TRUE) {

  w1 <- tryCatch(
    weightr::weightfunct(d, v, steps = c(0.025, 1), mods = NULL, weights = NULL,
      fe = FALSE, table = TRUE),
    error = function(e) NULL
  )

  res.NA <- data.frame(
    method = "3PSM",
    term = c("tau2", "b0", "pr.nonsig"),
    estimate = NA,
    std.error = NA,
    statistic = NA,
    p.value = NA,
    conf.low = NA,
    conf.high = NA
  )

  if (is.null(w1)) return(res.NA)

  # if <= 3 p-values in an interval: return NA
  p.table <- table(cut(w1$p, breaks=c(0, .025, 1)))
  if (any(p.table < min.pvalues)) {
    return(res.NA)
  } else {
    est <- w1[[2]]$par

    # compute standard errors from hessian
    std.err <- sqrt(abs(diag(solve(w1[[2]]$hessian))))


    res.wide <- data.frame(
      method = "3PSM",
      term = c("tau2", "b0", "pr.nonsig"),
      estimate = est,
      std.error = std.err,
      statistic = est/std.err,
      p.value = pnorm(est/std.err, lower.tail=FALSE)*2,
      conf.low = est + qnorm(.025)*std.err,
      conf.high = est + qnorm(1-.025)*std.err
    )
  }

   return(res.wide)
}




#' Make three attempts to fit with rma.uni
#'
#' - Use REML estimator which "is approximately unbiased and quite efficient"
#' (Viechtbauer 2005).
#' - "Knapp and Hartung (2003) method (knha = TRUE) is an adjustment
#' to the standard errors of the estimated coefficients, which helps to account for the
#' uncertainty in the estimate of tau2 and leads to different reference distributions.
#' Individual coefficients and confidence intervals are then based on the t-distribution
#' with k − p degrees of freedom, while the omnibus test statistic then uses an
#' F-distribution with m and k − p degrees of freedom (p being the total number of model
#' coefficients including the intercept if it is present)." (Viechtbauer 2005)
#' - The "amount of heterogeneity in the true effect is estimated to be tau^2" with
#' "measures for facilitating the interpretation of the estimated amount of heterogeneity":
#'    - "I2 statistic estimates (in percent) how much of the total variability in the
#'    effect size estimates (which is composed of heterogeneity and sampling variability)
#'    can be attributed to heterogeneity among the true effects (tau^2 = 0 therefore
#'    implies I2 = 0%)."
#'      * "In models with moderators, the I2 statistic indicates how much of the
#'        unaccounted variance in the observed effects or outcomes (which is composed of
#'        unaccounted variance in the true effects, that is, residual heterogeneity,
#'        plus sampling variance) can be attributed to residual heterogeneity"
#'    - "The H2 statistic is the ratio of the total amount of variability in the
#'    observed outcomes to the amount of sampling variability (tau^2 = 0 therefore implies
#'    H2 = 1)."
#'    - "It is important to realize, however, that tau^2, I2, and H2 are often estimated
#'    imprecisely, especially when the number of studies is small."
#'
attempt_uni_rma <- function(data, yi, vi, fmods = NULL, intercept = FALSE, wi = NULL, ...) {
  if (NROW(data) > 1) {
    if (is.null(fmods)) {
      intercept <- TRUE
    }

    # First attempt to fit rma model with REML estimator
    m <- try(metafor::rma.uni(yi = yi, vi = vi, intercept = intercept, mods = fmods,
      weights = wi, weighted = TRUE, data = data, method = "REML", test = "knha"),
      silent = TRUE)

    # Second attempt to fit rma model with REML estimator but adjusted steps and more
    # iterations; see http://www.metafor-project.org/doku.php/tips:convergence_problems_rma
    if (inherits(m, "try-error")) {
      m <- try(metafor::rma.uni(yi = yi, vi = vi, intercept = intercept, mods = fmods,
        weights = wi, weighted = TRUE, data = data, method = "REML", test = "knha",
        control = list(maxiter = 1000, stepadj = 0.5)), silent = TRUE)
    }

    # Third attempt to fit rma model but with non-iterative DerSimonian-Laird estimator
    if (inherits(m, "try-error")) {
      m <- try(metafor::rma.uni(yi = yi, vi = vi, intercept = intercept, mods = fmods,
        weights = wi, weighted = TRUE, data = data, method = "DL", test = "knha"),
        silent = TRUE)
    }

    if (inherits(m, "try-error")) {
      print(m)

    } else {
      m <- calc_overall_I2(m)
    }

  } else {
    m <- try(stop("NROW(data) < 2"), silent = TRUE)
  }

  m
}


#' Make three attempts to fit with rma.mv
#'
#' - Use REML estimator which "is approximately unbiased and quite efficient"
#' (Viechtbauer 2005).
#' - Use t-test, "which slightly mimics the Knapp and Hartung (2003) method by using a
#' t-distribution with k-p degrees of freedom for tests of individual coefficients
#' and confidence intervals and an F-distribution with m and k-p degrees of freedom
#' (p being the total number of model coefficients including the intercept if it is
#' present) for the omnibus test statistic."
#'
#' Phylogenetic correction is based on code from OpenMEE which calls metafor::rma.mv:
#'   - https://github.com/gdietz/OpenMEE/blob/master/R/openmeer/R/phyloma.r
#'     -> function 'phylo.meta.analysis'
#'   - https://github.com/gdietz/OpenMEE/blob/523c85b1d4ba9fe3605e3fa9ae8e80c111f1c8fb/python_to_R.py
#'     -> functions '_run_phylo_ma' and 'run_phylo_ma'
#' more infos: https://www.zoology.ubc.ca/~schluter/R/phylogenetic/
attempt_multilevel_rma <- function(data, yi, vi, fmods = NULL, intercept = FALSE, wi = NULL,
  account_repeats = FALSE, account_phylo = FALSE, tree = NULL, evolution = c("BM", "OU"),
  evolution_BM_lambda = 1, evolution_OU_alpha = 0) {

  omax <- 1000

  if (NROW(data) > 1) {
    if (is.null(fmods)) {
      intercept <- TRUE
    }

    # Construct random factors
    frand <- list()
    Nos <- list(
      Marker_cat = table(as.character(data[, "Marker_cat"])),
      fSpecies = table(as.character(data[, "fSpecies"])),
      Marker_x_Species = table(as.character(data[, "Marker_cat"]),
        as.character(data[, "fSpecies"])))

    # variance structure of random 'inner' argument
    # see 'Meta-Regression with All Studies but Different Amounts of (Residual) Heterogeneity'
    # at http://www.metafor-project.org/doku.php/tips:comp_two_independent_estimates
    # where the example uses 'DIAG' (e.g., rho == 0), but my test showed that estimates of
    # rho are +-1 which makes sense since we expect markers to be strongly correlated
    inner_struct <- if (any(Nos[["Marker_cat"]] > 1)) "HCS" else "CS"

    # Only up to two '~ inner | outer' formulas allowed in the 'random' argument
    inner_N <- 0

    if (account_phylo && !is.null(tree)) {
      frand <- c(frand, ~ 1 | phylogenyVariance)
    }

    if (account_repeats) {
      frand <- c(frand, ~ 1 | betweenStudyVariance)
      if (any(Nos[["Marker_cat"]] > 1) && inner_N < 2) {
        # Random-slopes (multivariate model) with no correlation among slope and intercept
        frand <- c(frand, ~ Marker_cat - 1 | betweenStudyVariance)
        inner_N <- inner_N + 1
      }

      if (any(Nos[["fSpecies"]] > 1)) {
        frand <- c(frand, ~ 1 | fSpecies)
        if (any(Nos[["Marker_x_Species"]] > 1) && inner_N < 2) {
          # Random-slopes (multivariate model) with no correlation among slope and intercept
          frand <- c(frand, ~ Marker_cat - 1 | fSpecies)
          inner_N <- inner_N + 1
        }
      }
    }

    # Construct phylogenetic correlation structure
    if (account_phylo && !is.null(tree)) {
      # Check that species in tree and data.frame have the same order
      stopifnot(identical(unique(data[, "Species"]), unique(tree[["tip.label"]])))

      temp_cor <- if (identical(evolution, "BM")) {
          # "The correlation structure from the present model is derived from the Brownian
          # motion model by multiplying the off-diagonal elements (i.e., the covariances) by
          # lambda. The variances are thus the same than for a Brownian motion model."
        	ape::corPagel(value = evolution_BM_lambda, phy = tree, fixed = TRUE)

        } else if (identical(evolution, "OU")) {
          # "Martins and Hansen's (1997) covariance structure [representing the Ornstein-Uhlenbeck model]:
          #     Vij = gamma . exp(-alpha . tij)
          # where tij is the phylogenetic distance between taxa i and j and gamma is a constant."
        	ape::corMartins(value = evolution_OU_alpha, phy = tree, fixed = TRUE)
        }

      temp_sp <- data.frame(tree[["tip.label"]])
      rownames(temp_sp) <- tree[["tip.label"]]

      temp_cor2 <- nlme::Initialize(temp_cor, temp_sp)
      temp_cor3 <- nlme::corMatrix(temp_cor2, corr = TRUE)

      # Expand phylogenetic correlation structure to account for multiple species entries
      # in dataset
      temp_n <- table(data[, "Species"])
      if (any(temp_n > 1)) {
        ids <- rep(seq_along(temp_n), temp_n)
        temp_cor3 <- temp_cor3[ids, ids]
      }

      corR <- list(phylogenyVariance = temp_cor3)

    } else {
      corR <- NULL
    }

    # First attempt with `nlminb` but with a maximum of 1000 instead of 150 iterations
    m <- try(metafor::rma.mv(yi = yi, V = vi, intercept = intercept, mods = fmods,
      random = frand, struct = inner_struct, W = wi, R = corR, data = data,
      method = "REML", test = "t", control =
        list(iter.max = omax, eval.max = omax)),
      silent = TRUE)

if (FALSE) {
  res[["ml2"]][["mod"]][["fit"]] <- m
  res[["ml3"]][["mod"]][["fit"]] <- m
  res[["p"]][["mod"]][["fit"]] <- m
  temp <- summary(res[["uni"]][["mod"]][["fit"]])
  temp <- summary(res[["ml3"]][["mod"]][["fit"]])
  temp <- summary(res[["ml2"]][["mod"]][["fit"]])
  temp <- summary(res[["p"]][["mod"]][["fit"]])
  temp <- m

  metafor::forest(x = temp$b, ci.lb = temp$ci.lb, ci.ub = temp$ci.ub,
    xlab = "Estimates", # clim = c(-1, 1),
    slab = rownames(temp$b))
}

    # Second attempt to fit rma model but with "Nelder-Mead" from 'optim'
    if (inherits(m, "try-error")) {
      m <- try(metafor::rma.mv(yi = yi, V = vi, intercept = intercept, mods = fmods,
        random = frand, struct = inner_struct, W = wi, R = corR, data = data,
        method = "REML", test = "t", control =
          list(optimizer = "optim", optmethod = "Nelder-Mead", maxit = omax)),
        silent = TRUE)
    }

    # Third attempt to fit rma model but with the quasi-Newton method "BFGS" from 'optim'
    if (inherits(m, "try-error")) {
      m <- try(metafor::rma.mv(yi = yi, V = vi, intercept = intercept, mods = fmods,
        random = frand, struct = inner_struct, W = wi, R = corR, data = data,
        method = "REML", test = "t", control =
          list(optimizer = "optim", optmethod = "BFGS", maxit = omax)),
        silent = TRUE)
    }

    # Fourth attempt to fit rma model but with "ucminf"
    if (inherits(m, "try-error")) {
      m <- try(metafor::rma.mv(yi = yi, V = vi, intercept = intercept, mods = fmods,
        random = frand, struct = inner_struct, W = wi, R = corR, data = data,
        method = "ML", test = "t", control =
          list(optimizer = "ucminf", maxit = omax, maxeval = omax)),
        silent = TRUE)
    }

    if (inherits(m, "try-error")) {
      print(m)

    } else {
      m <- calc_overall_I2(m)
    }

  } else {
    m <- try(stop("NROW(data) < 2"), silent = TRUE)
  }

  m
}


has_fit <- function(x) {
  !inherits(x, "try-error") && !is.null(x) && !is.na(x)
}



#' Random/mixed-effect models with Knapp & Hartung adjustment: not accounting for
#' non-independence
fit_uni <- function(mdata, yin, vin, win, intercept = FALSE, interaction = TRUE) {
  # Prepare output containers
  out <- list(nomod = list(), mod = list())

  # Prepare data
  data <- mdata[["all"]][["d"]]

  if (NROW(data) > 1) {
    is_mod_redundant <- mdata[["all"]][["mods"]][["is_redundant"]]
    mods <- mdata[["all"]][["mods"]][["names"]][!is_mod_redundant]
    is_mod_continuous <- mdata[["all"]][["mods"]][["is_continuous"]][!is_mod_redundant]

    if (length(mods) > 1 && interaction) {
      modX <- name_interaction(mods)
      data[, modX] <- factor(apply(data[, mods], 1, paste, collapse = "_x_"))
      fmods <- as.formula(paste("~", modX,
          if (!intercept && !any(is_mod_continuous)) "-1"))

    } else if (length(mods) > 0) {
      fmods <- as.formula(paste("~", paste(mods, collapse = "+"),
          if (!intercept && !any(is_mod_continuous)) "-1"))

    } else {
      fmods <- NULL
    }

    yi <- data[, yin]
    vi <- data[, vin]

    wi <- if (!identical(win, "winvar")) {
        data[, win]
      } else NULL

  }


  # Model without moderator = overall random-effect model
  if (NROW(data) > 1) {
    out[["nomod"]][["fit"]] <- attempt_uni_rma(data = data, yi = yi, vi = vi,
      fmods = NULL, wi = wi)
  }

  out[["nomod"]][["has_fit"]] <- has_fit(out[["nomod"]][["fit"]])

  if (NROW(data) > 1) {
    # Model with moderator and no intercept = mixed-effect model
    #   - no intercept
    #     - level estimates represent the mean estimates of each level of moderator
    #   - intercept (if turned on)
    #     - intercept estimate = mean estimate of reference level of moderator
    #     - level estimates = estimate of difference to reference (i.e., level - reference level)
    # Note: models with/without intercept produce the same tau, I^2, H^2, and QE estimates
    #   and test statistics; test for moderators obviously differs
    out[["mod"]][["fit"]] <- attempt_uni_rma(data = data, yi = yi, vi = vi, fmods = fmods,
      intercept = FALSE, wi = wi)
    out[["mod"]][["has_fit"]] <- has_fit(out[["mod"]][["fit"]])

    #--- Test each moderator
    if (out[["mod"]][["has_fit"]]) {
      x <- mdata[["all"]][["mods"]]
      for (im in seq_along(mods)) {
        if (!is_mod_continuous[im] &&
          any(mdata[["all"]][["mods"]][["n_cats"]][!is_mod_redundant] >= 2)) {

          # All pairwise differences among moderator categories with Holm's method
          temp <- rep(1, length = x[["n_cats"]][mods[im]])
          names(temp) <- x[["cats"]][[mods[im]]]
          ctemp <- multcomp::contrMat(n = temp, type = "Tukey")
          temp <- try(multcomp::glht(out[["mod"]][["fit"]], linfct = ctemp))
          out[["mod"]][["pairdiffs"]][[mods[im]]] <- list()
          if (!inherits(temp, "try-error")) {
            out[["mod"]][["pairdiffs"]][[mods[im]]][["holm"]] <-
              summary(temp, test = adjusted("holm"))
            out[["mod"]][["pairdiffs"]][[mods[im]]][["Wald"]] <-
              anova(out[["mod"]][["fit"]], L = ctemp)
          }
        }

        # Wald-type omnibus test factors with multiple treatment levels
        btt <- grep(mods[im], names(coef(out[["mod"]][["fit"]])))
        if (length(btt) > 0) {
          out[["mod"]][["omnibus"]][[mods[im]]][["Wald"]] <-
            try(anova(out[["mod"]][["fit"]], btt = btt), silent = TRUE)
        }
      }
    }

    # Percentage of total amount of heterogeneity can be accounted for by including moderators
    out[["Htotal_explained"]] <- calc_explained_heterogeneity(out[["nomod"]], out[["mod"]])

  } else {
    out[["mod"]][["has_fit"]] <- FALSE
  }

  out
}



#' Random/mixed-effect models with Knapp & Hartung adjustment: accounting for
#' non-independence in studies, markers, and species phylogenies
fit_multilevel <- function(sdata, yin, vin, win, intercept = FALSE, interaction = TRUE,
  account_repeats = FALSE, account_phylo = FALSE, evolution = c("BM", "OU"),
  tree_branchLengths = FALSE) {

  evolution <- match.arg(evolution)

  # Prepare output containers
  out <- list(nomod = list(), mod = list(), mod_uni = list())

  # Prepare data
  data <- sdata[["d"]]

  if (NROW(data) > 1) {
    tree <- if (tree_branchLengths) {
        sdata[["tree"]]
      } else {
        sdata[["tree_woBL"]]
      }

    is_mod_redundant <- sdata[["mods"]][["is_redundant"]]
    mods <- sdata[["mods"]][["names"]][!is_mod_redundant]
    is_mod_continuous <- sdata[["mods"]][["is_continuous"]][!is_mod_redundant]

    if (length(mods) > 1 && interaction) {
      modX <- name_interaction(mods)
      data[, modX] <- factor(apply(data[, mods], 1, paste, collapse = "_x_"))
      fmods <- as.formula(paste("~", modX,
          if (!intercept && !any(is_mod_continuous)) "-1"))

    } else if (length(mods) > 0) {
      fmods <- as.formula(paste("~", paste(mods, collapse = "+"),
          if (!intercept && !any(is_mod_continuous)) "-1"))

    } else {
      fmods <- NULL
    }

    yi <- data[, yin]
    vi <- data[, vin]

    wi <- if (!identical(win, "winvar")) {
        data[, win]
      } else NULL

  } else {
    yi <- vi <- wi <- NULL
  }

  # Model without moderator = overall random-effect model
  if (NROW(data) > 1) {
    out[["nomod"]][["fit"]] <- attempt_multilevel_rma(data = data, yi = yi,
      vi = vi, fmods = NULL, wi = wi,
      account_repeats = account_repeats, account_phylo = account_phylo, tree = tree,
      evolution = evolution, evolution_BM_lambda = 1, evolution_OU_alpha = 0)
  }

  out[["nomod"]][["has_fit"]] <- has_fit(out[["nomod"]][["fit"]])

  if (NROW(data) > 1) {

    # Model with moderator and no intercept = mixed-effect model
    #   - no intercept
    #     - level estimates represent the mean estimates of each level of moderator
    #   - intercept (if turned on)
    #     - intercept estimate = mean estimate of reference level of moderator
    #     - level estimates = estimate of difference to reference (i.e., level - reference level)
    # Note: models with/without intercept produce the same tau, I^2, H^2, and QE estimates
    #   and test statistics; test for moderators obviously differs
    out[["mod"]][["fit"]] <- attempt_multilevel_rma(data = data, yi = yi, vi = vi,
      fmods = fmods, intercept = FALSE, wi = wi,
      account_repeats = account_repeats, account_phylo = account_phylo, tree = tree,
      evolution = evolution, evolution_BM_lambda = 1, evolution_OU_alpha = 0)
    out[["mod"]][["has_fit"]] <- has_fit(out[["mod"]][["fit"]])

    # Run analysis as-if with rma.uni
    out[["mod_uni"]][["fit"]] <- attempt_multilevel_rma(data = data, yi = yi, vi = vi,
      fmods = fmods, intercept = FALSE, wi = wi,
      account_repeats = FALSE, account_phylo = FALSE, tree = tree,
      evolution = evolution, evolution_BM_lambda = 1, evolution_OU_alpha = 0)
    out[["mod_uni"]][["has_fit"]] <- !inherits(out[["mod_uni"]][["fit"]], "try-error")


    #--- Test each moderator
    if (out[["mod"]][["has_fit"]]) {
      x <- sdata[["mods"]]
      for (im in seq_along(mods)) {
        if (!is_mod_continuous[im] &&
          any(sdata[["mods"]][["n_cats"]][!is_mod_redundant] >= 2)) {

          # All pairwise differences among moderator categories with Holm's method
          temp <- rep(1, length = x[["n_cats"]][mods[im]])
          names(temp) <- x[["cats"]][[mods[im]]]
          ctemp <- multcomp::contrMat(n = temp, type = "Tukey")
          temp <- try(multcomp::glht(out[["mod"]][["fit"]], linfct = ctemp))
          out[["mod"]][["pairdiffs"]][[mods[im]]] <- list()
          if (!inherits(temp, "try-error")) {
            out[["mod"]][["pairdiffs"]][[mods[im]]][["holm"]] <-
              summary(temp, test = adjusted("holm"))
            out[["mod"]][["pairdiffs"]][[mods[im]]][["Wald"]] <-
              anova(out[["mod"]][["fit"]], L = ctemp)
          }
        }

        # Wald-type omnibus test factors with multiple treatment levels
        btt <- grep(mods[im], names(coef(out[["mod"]][["fit"]])))
        if (length(btt) > 0) {
          out[["mod"]][["omnibus"]][[mods[im]]][["Wald"]] <-
            try(anova(out[["mod"]][["fit"]], btt = btt), silent = TRUE)
        }
      }
    }

    # Percentage of total amount of heterogeneity can be accounted for by including moderators
    out[["Htotal_explained"]] <- calc_explained_heterogeneity(out[["nomod"]], out[["mod"]])

  } else {
    out[["mod"]][["has_fit"]] <- FALSE
  }

  out
}


# Create name uniquely describing analysis factors and levels
filetag_ID <- function(interaction_wHabitat3, fragment_sizes.,
  withControlsL., only_wo_controls, temp_withNonStandardizedL, cor_methods.,
  cor_transforms., weight_methods.) {

  temp_withControlsL <- if (only_wo_controls) FALSE else withControlsL.

  stopifnot(lengths(list(fragment_sizes., temp_withControlsL, cor_methods.,
    cor_transforms., weight_methods.)) == 1L)

  temp1 <- if (temp_withControlsL) "withControls" else "withoutControls"
  temp_cor_transforms <- if (cor_methods. == "pearson") {
      cor_transforms.
    } else {
      "ztransform"
    }
  temp2 <- if (temp_withNonStandardizedL) {
      "withNonStandardized"
    } else {
      "onlyStandardized"
    }

  p1 <- paste(fragment_sizes., temp1, temp2, sep = "_")
  p2 <- paste0("COR", cor_methods., "-", temp_cor_transforms)
  p3 <- weight_methods.

  if (interaction_wHabitat3) {
    paste(p1, "modXhabitat", p2, p3, sep = "_")

  } else {
    paste(p1, p2, p3, sep = "_")
  }
}

