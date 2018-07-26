
#--- SETTINGS

anas <- c("all", "aphylo", "pphylo")

lab_responses <- data.frame(
  code = c("mA", "He+Sh", "PLp", "Fis"),
  out = c("A", "He", "PLP", "Fis"),
  hypothesis = c("left", "left", "left", "right"),
  stringsAsFactors = FALSE)
#      temp <- lab_responses[responses[ir] == lab_responses[, "code"], "out"]


lab_fittypes <- data.frame(
  code = c("all_un", "all_mv", "aphylo", "pphylo"),
  out = c("RE", "ML", "MLAP", "MLPP"),
    #inspired by http://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=4
  col = c("orangered", "darkslateblue", "purple2", "springgreen4"),
  stringsAsFactors = FALSE)
#      temp <- lab_fittypes[match(fit_types, lab_fittypes[, "code"]), "out"]

figwidths <- c(3, 6) #inches for Ecological Applications

lab_moderators <- data.frame(
  code = c("Plants_Sex", "Plants_Self_compatibility", "Plants_Mating_system",
    "Plants_Pollination_syndrome", "Plants_Seed_dispersal_vector",
    "EffectAgeCat2_ord_yrs", "Effect_NGenSinceFrag_ord"),
  out = c("Mode of reproduction", "Compatibility system", "Reproductive system",
    "Pollination syndrome", "Seed dispersal",
    "Abs. fragmentation age", "Rel. fragmentation age"),
  stringsAsFactors = FALSE)
#      temp <- lab_moderators[match(moderators, lab_moderators[, "code"]), "out"]

lab_modlevels <- data.frame(
  code = c("sex", "Tf/Ng"),
  out = c("Sexual", "# gen"),
  stringsAsFactors = FALSE)
#      temp <- lab_modlevels[match(modcats, lab_modlevels[, "code"]), "out"]


# Paths
dir_prj1 <- "Prj04_GenetFragment_MetaAnalysis"
dir_ana <- file.path(dir_prj1, "3_Analysis")

#--- Setup
source(file.path(dir_ana, "Step3_0-Analysis-settings.R"))

# Output path
dir_ms_out <- file.path(dir_res_, "Datasets_MS")
dir.create(dir_ms_out, recursive = TRUE, showWarnings = FALSE)


#--- Helper functions
do_ms_loops <- function(fun, args, args_along = NULL, do_targets = TRUE,
  do_interactions = c(TRUE, FALSE)) {

  stopifnot(lengths(args_along) %in% c(0, 5))

  do_responses <- list(`resp-gendiv` = c("mA", "He+Sh", "PLp"), `resp-inbred` = "Fis")

  for (ir in seq_along(do_responses)) for (intH3 in do_interactions) {
    # 0_Miscellaneous
    args1 <- c(args,
        responses = list(do_responses[[ir]]),
        moderators = list(c("Marker_cat", "Effect_latitude_bin")),
        ftag = list(paste("0_Miscellaneous",
          if (do_targets) "target" else "all", names(do_responses)[ir], sep = "-")),
        interaction_wHabitat3 = list(intH3)
      )
    for (k in seq_along(args_along)) {
      args1[[names(args_along)[k]]] <- args_along[[k]][[1L]]
    }
    do.call(fun, args = args1)

    # 1_DecreasedGeneticDiversity && 2_IncreasedInbreeding
    args2 <- c(args,
        responses = list(do_responses[[ir]]),
        moderators = list(moderators_Hs),
        ftag = list(paste("1-2_Diversity&Inbreeding",
          if (do_targets) "target" else "all", names(do_responses)[ir], sep = "-")),
        interaction_wHabitat3 = list(intH3)
      )
    for (k in seq_along(args_along)) {
      args2[[names(args_along)[k]]] <- args_along[[k]][[2L]]
    }
    do.call(fun, args = args2)

    args2b <- c(args,
        responses = list(do_responses[[ir]]),
        moderators = list(moderators_Hs_highres),
        ftag = list(paste("1-2_Diversity&Inbreeding", "highres",
          if (do_targets) "target" else "all", names(do_responses)[ir], sep = "-")),
        interaction_wHabitat3 = list(intH3)
      )
    for (k in seq_along(args_along)) {
      args2b[[names(args_along)[k]]] <- args_along[[k]][[2L]]
    }
    do.call(fun, args = args2b)

    args3a <- c(args,
        responses = list(do_responses[[ir]]),
        moderators = list(moderators_plant_traits1),
        ftag = list(paste("1-2_Diversity&Inbreeding", "PlantTraits1",
          if (do_targets) "target" else "all", names(do_responses)[ir], sep = "-")),
        interaction_wHabitat3 = list(intH3)
      )
    for (k in seq_along(args_along)) {
      args3a[[names(args_along)[k]]] <- args_along[[k]][[3L]]
    }
    do.call(fun, args = args3a)

    args3b <- c(args,
      responses = list(do_responses[[ir]]),
      moderators = list(moderators_plant_traits2),
      ftag = list(paste("1-2_Diversity&Inbreeding", "PlantTraits2",
        if (do_targets) "target" else "all", names(do_responses)[ir], sep = "-")),
      interaction_wHabitat3 = list(intH3)
    )
    for (k in seq_along(args_along)) {
      args3b[[names(args_along)[k]]] <- args_along[[k]][[3L]]
    }
    do.call(fun, args = args3b)

    args4 <- c(args,
        responses = list(do_responses[[ir]]),
        moderators = list(moderators_animal_traits),
        ftag = list(paste("1-2_Diversity&Inbreeding", "AnimalTraits",
          if (do_targets) "target" else "all", names(do_responses)[ir], sep = "-")),
        interaction_wHabitat3 = list(intH3)
      )
    for (k in seq_along(args_along)) {
      args4[[names(args_along)[k]]] <- args_along[[k]][[4L]]
    }
    do.call(fun, args = args4)

    # 3_TimeEffects
    args5 <- c(args,
        responses = list(do_responses[[ir]]),
        moderators = list(moderators_Htime),
        ftag = list(paste("3_TimeEffects",
          if (do_targets) "target" else "all", names(do_responses)[ir], sep = "-")),
        interaction_wHabitat3 = list(intH3)
      )
    for (k in seq_along(args_along)) {
      args5[[names(args_along)[k]]] <- args_along[[k]][[5L]]
    }
    do.call(fun, args = args5)
  }

  invisible(TRUE)
}



# Code from examples of ?toupper
capwords <- function(s, strict = FALSE) {
    cap <- function(s) {
      paste(toupper(substring(s, 1, 1)), {
        s <- substring(s, 2)
        if (strict) tolower(s) else s
      }, sep = "", collapse = " ")
    }
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}


blanks <- function(n) strrep(" ", n)

#' Remove digits if there are any immediately following the third character as in `ExH1`
remove_modcatID <- function(x) {
  m <- gregexpr("[[:alpha:]]x[[:alpha:]][[:digit:]]+", x)
  if (any(sapply(m, function(x) any(x > 0)))) {
    m <- lapply(m, function(x) { # we only want to remove the digits
      if (x > 0) {
        x <- x + 3
        attr(x, "match.length") <- attr(x, "match.length") - 3
      }
      x
    })
    regmatches(x, m) <- ""
  }
  x
}



clean_labels <- function(x, remove_digits = FALSE, break_lines = TRUE,
  remove_interaction = FALSE) {

  labs <- x

  ids <- match(labs, lab_moderators[, "code"], nomatch = 0)
  labs[ids > 0] <- lab_moderators[ids, "out"]

  ids <- match(labs, lab_modlevels[, "code"], nomatch = 0)
  labs[ids > 0] <- lab_modlevels[ids, "out"]

  labs <- gsub("_", " ", labs)

  if (remove_digits) {
    # remove digits at end of words but not those between () and/or [] or numbers after punctuation
    #--- Code based on examples of ?regmatches
    # Create a copy with the matched parts blanked out.
    save_labs <- labs

    # First, match anything between any sets and combinations of parentheses, brackets, or curly braces
    m <- gregexpr("[[({][^][{}()]*[])}]", save_labs)
    if (any(sapply(m, function(x) any(x > 0)))) {
      # Create a copy with the matched parts blanked out.
      regmatches(save_labs, m) <- Map(blanks, lapply(regmatches(save_labs, m), nchar))
    }

    # Second, match numbers following punctuation, but not those at end of a word
    m <- gregexpr("[[:punct:]]{1,2}[[:digit:]]+", save_labs)
    if (any(sapply(m, function(x) any(x > 0)))) {
      # Create a copy with the matched parts blanked out.
      regmatches(save_labs, m) <- Map(blanks, lapply(regmatches(save_labs, m), nchar))
    }

    # Compute the positions of the matches of the digits to remove (note that we cannot call
    # strsplit() on labs with match data from save_labs).
    m <- gregexpr("[\b[:digit:]]+", save_labs)
    if (any(sapply(m, function(x) any(x > 0)))) {
      # And finally extract the non-matched parts and put parts back together.
      labs <- sapply(regmatches(labs, m, invert = TRUE), paste, collapse = "")
    }
  }

  labs <- gsub("-plants", " Plants", labs)
  labs <- gsub("Nonwoody", "Non-woody", labs, ignore.case = TRUE)
  labs <- gsub("cat\\>", "", labs, ignore.case = TRUE)
  labs <- gsub("\\<ord\\>", "", labs, ignore.case = TRUE)
  labs <- gsub("\\<bin\\>", "", labs, ignore.case = TRUE)
  labs <- gsub("\\<effect\\>", "", labs, ignore.case = TRUE)
  labs <- gsub("  ", " ", labs)
  labs <- gsub("dominated", "dom.", labs)

  labs <- capwords(labs)
  labs <- gsub("\\<And\\>", "and", labs, ignore.case = FALSE)
  labs <- gsub("\\<Or\\>", "or", labs, ignore.case = FALSE)
  labs <- gsub("\\<Of\\>", "of", labs, ignore.case = FALSE)
  labs <- gsub("\\<Yrs\\>", "yrs", labs, ignore.case = FALSE)
  labs <- gsub("\\<sc\\>", "SC", labs, ignore.case = TRUE)
  labs <- gsub("\\<si\\>", "SI", labs, ignore.case = TRUE)
  labs <- gsub("\\<Yrs\\>", "yrs", labs, ignore.case = FALSE)
  labs <- gsub(" X ", " x ", labs)

  labs <- gsub("Moist And Wet", "Moist/Wet", labs)
  labs <- gsub("Tf/Ng", "# Gen.", labs)

  # Remove duplicate words
  temp <- strsplit(labs, split = " x ", fixed = TRUE)
  labs <- sapply(temp, function(x) {
        if (length(x) == 2 && (remove_interaction || grepl(x[2], x[1]))) {
          x[1]
        } else {
          paste0(x, collapse = " x ")
        }})

  if (break_lines) {
    labs <- gsub("Animals ", "Animals\n", labs)
    labs <- gsub("Plants ", "Plants\n", labs)

  } else {
    # Remove 'animals' or 'plants' if they occur in each of labs
    if (all(grepl("\\<Animals\\>", labs))) {
      labs <- gsub("\\<Animals\\>", "", labs)
    }
    if (all(grepl("\\<Plants\\>", labs))) {
      labs <- gsub("\\<Plants\\>", "", labs)
    }
  }

  labs <- trimws(labs)
  labs
}



msres_addsensitivity <- function(msdata, sided, interaction_wHabitat3 = FALSE, dir_res,
  ftag) {
  stopifnot(sided %in% c("left", "right"))

  temp_xHabitat3 <- if (interaction_wHabitat3) "modXhabitat" else ""
  tag_fix <- if (interaction_wHabitat3) {
      paste0(ftag, "_", temp_xHabitat3)
    } else {
      ftag
    }
  fsens <- file.path(dir_res, "Datasets_Sensitivity",
    paste0("Data_Sensitivity_", tag_fix, ".rds"))

  if (file.exists(fsens)) {
    sensdata <- readRDS(fsens)
  } else {
    stop(paste("File not found:", shQuote(basename(fsens))))
  }

  temp <- dimnames(msdata[["data"]])
  moderators <- temp[[1]]
  responses <- temp[[2]]
  sensmods <- dimnames(sensdata[["sensitivity"]])[[2]]

  dims_msres <- c(nrow = length(moderators), ncol = length(responses))
  sens <- array(list(), dim = c(2, length(moderators), length(responses)),
    dimnames = list(c("agreement", "target"), moderators, responses))

  idtest <- sapply(sided, function(x)
    switch(EXPR = x, left = "CI_lt0", right = "CI_gt0", NA))
  stopifnot(length(idtest) == length(responses))
  N <- sensdata[["N"]]

  for (ir in seq_along(responses)) {
    for (ig in seq_along(moderators)) {
      idmods <- grep(msdata[["moderatorsX"]][moderators[ig]], sensmods)

      if (length(idmods) > 0) {
        temp <- sensdata[["sensitivity"]][responses[ir], , , idtest[ir]]
        sens["agreement", ig, ir] <- list(temp[idmods, , drop = FALSE] / N)

        temp <- sensdata[["target"]][responses[ir], , , idtest[ir]]
        sens["target", ig, ir] <- list(target = temp[idmods, , drop = FALSE])
      }
    }
  }

  c(msdata, list(sensitivity = sens))
}


# Group meta-analysis output in units for display items
msres_getdatatogether <- function(responses, moderators, interaction_wHabitat3 = FALSE,
  fragment_sizes. = fragment_sizes, withControlsL. = withControlsL, only_wo_controls = TRUE,
  withNonStandardizedL. = withNonStandardizedL, only_useadj_standardized = TRUE,
  cor_methods. = cor_methods, cor_transforms. = cor_transforms,
  weight_methods. = weight_methods, dir_res, dir_out, ftag, ...) {

  temp_withControlsL <- if (only_wo_controls) FALSE else withControlsL.
  temp_xHabitat3 <- if (interaction_wHabitat3) "modXhabitat" else ""

  # Only 'responses' and 'moderators' allowed to vary
  stopifnot(lengths(list(fragment_sizes., temp_withControlsL, cor_methods.,
    cor_transforms., weight_methods.)) == 1L)

  # File identification
  temp_withNonStandardizedL <- !only_useadj_standardized
  tag_fix <- filetag_ID(interaction_wHabitat3, fragment_sizes.,
    withControlsL., only_wo_controls, temp_withNonStandardizedL, cor_methods.,
    cor_transforms., weight_methods.)

  fmsdata <- file.path(dir_out, paste0("MSdata_", ftag, "_", tag_fix, ".rds"))

  if (file.exists(fmsdata)) {
    ms_data <- readRDS(fmsdata)
    do_load <- !all(c("data", "res", "sensitivity", "cats", "n_cats", "n_cats_max",
      "modcats", "tag_fix") %in% names(ms_data))

  } else {
    do_load <- TRUE
  }

  if (do_load) {
    #--- Load analysis outputs
    dims_msres <- c(nrow = length(moderators), ncol = length(responses))
    res <- data <- matrix(list(), nrow = dims_msres[1], ncol = dims_msres[2],
      dimnames = list(moderators, responses))

    missing_N <- 0

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

      ftag_fix <- filetag_ID(interaction_wHabitat3, fragment_sizes.,
        withControlsL., only_wo_controls, temp_withNonStandardizedL, cor_methods.,
        cor_transforms., weight_methods.)

      for (ig in seq_along(moderators)) {
        fname <- paste(responses[ir], moderators[ig], "by", ftag_fix, sep = "_")
        fdata <- file.path(dir_data, paste0(fname, "_data.rds"))
        ftemp <- file.path(dir_resout, paste0(fname, "_output.rds"))

        if (file.exists(fdata)) {
          data[ig, ir] <- list(readRDS(fdata))
        } else {
          missing_N <- missing_N + 1
          print(paste("File not found:", shQuote(basename(fdata))))
        }

        if (file.exists(ftemp)) {
          res[ig, ir] <- list(readRDS(ftemp))
        } else {
          missing_N <- missing_N + 1
          print(paste("File not found:", shQuote(basename(ftemp))))
        }
      }
    }

    if (missing_N < length(moderators) * length(responses)) {
      #--- Moderator levels and numbers
      cats <- array(list(), dim = length(moderators), dimnames = list(moderators))
      n_cats <- array(NA, dim = c(length(moderators), length(responses), length(anas)),
        dimnames = list(moderators, responses, anas))

      for (ir in seq_along(responses)) {
        for (ig in seq_along(moderators)) {
          x <- data[ig, ir][[1]]
          cats[ig] <- list(sort(unique(c(unlist(cats[ig]),
            unlist(sapply(anas, function(k) x[[k]][["mods"]][["cats"]]))))))
          n_cats[ig, ir, ] <- sapply(anas, function(k) x[[k]][["mods"]][["n_cats"]])
        }
      }

      #n_cats_max <- apply(apply(n_cats, 1:2, max), 1, max)
      n_cats_max <- lengths(cats)
      stopifnot(lengths(cats) >= apply(apply(n_cats, 1:2, max), 1, max))

      moderatorsX <- sapply(moderators, function(ig) {
        data[ig, 1][[1]][[1]][["mods"]][["names"]]
      })
      modcats <- paste0(rep(paste0(moderatorsX, seq_along(moderatorsX)), n_cats_max),
        unlist(cats))


      #--- Output object
      ms_data <- list(data = data, res = res,
        cats = cats, n_cats = n_cats, n_cats_max = n_cats_max,
        modcats = modcats, moderators = moderators, moderatorsX = moderatorsX,
        tag_fix = tag_fix)

      saveRDS(ms_data, file = fmsdata)

    } else {
      stop("No data found.")
    }
  }

  invisible(ms_data)
}




has_fit <- function(x) {
  !is.null(x) && !inherits(x, "try-error")
}
