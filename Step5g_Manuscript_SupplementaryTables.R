
#--- SETTINGS
if (!exists("tag_usedate")) tag_usedate <- "20180705" #"20180112"
redo <- FALSE
do_ms <- TRUE
do_targets <- TRUE

do_Appendix3_FigsS1and2 <- TRUE
do_Appendix2_TableS1 <- TRUE
do_Data1_TablesS1and2 <- TRUE
do_Data1_TableS3 <- TRUE
do_Appendix5_TableS2 <- TRUE
do_Appendix5_TableS3 <- TRUE

# Paths
dir_prj1 <- "Prj04_GenetFragment_MetaAnalysis"
dir_ana <- file.path(dir_prj, "3_Analysis")

#--- Setup
source(file.path(dir_ana, "Step5-0_Manuscript-methods.R"))

foutw <- file.path(dir_prj, "2_Data", "5_DataCleaned",
  paste0(tag_usedate, "_ExtractedDataCleaned3_wide.rds"))
dat0 <- readRDS(foutw)

fin <- file.path(dir_prj, "2_Data", "5_DataCleaned",
  paste0(tag_usedate, "_ExtractedDataCleaned3_long.rds"))
datl <- readRDS(fin)

dat_citation <- read.csv(file.path(dir_prj, "2_Data", "1_Literature",
  "20170111_CompleteDB_AllRecords.csv"))

#--------------------------------------
#--- Number of animal and plant species
temp <- unique(dmoderators[, c("Organism_group4", "Species_resolved")])
ttemp <- table(temp[, "Organism_group4"])
print(ttemp)
# animals nonwoody-plants    woody-plants
# 52              40              43


#--------------------------------------
#--- Phylogeny figures

if (do_Appendix3_FigsS1and2) {
  source(file.path(dir_ana, "Step3_2_Analysis-Phylogeny.R"))

  dir_temp <- file.path(dir_res0, "Figs_Phylogeny")
  dir.create(dir_temp, recursive = TRUE, showWarnings = FALSE)

  pdf(height = 5, width = 6, file = file.path(dir_temp, "Fig_Phylogeny_animals_woBL.pdf"))
  ape::plot.phylo(dphyloanimals, use.edge.length = FALSE, node.depth = 2,
    no.margin = TRUE, cex = 0.5, label.offset = 0.25)
  dev.off()

  pdf(height = 8, width = 6, file = file.path(dir_temp, "Fig_Phylogeny_plants_woBL.pdf"))
  ape::plot.phylo(dphyloplants, use.edge.length = FALSE, node.depth = 2,
    no.margin = TRUE, cex = 0.5, label.offset = 0.25)
  dev.off()

  pdf(height = 8, width = 6, file = file.path(dir_temp, "Fig_Phylogeny_plants_wBL.pdf"))
  ape::plot.phylo(dphyloplants,
    no.margin = TRUE, cex = 0.5, label.offset = 5)
  ape::add.scale.bar()
  dev.off()
}


#--------------------------------------
#--- Table Habitat classification

if (do_Appendix2_TableS1) {
  utemp <- unique(dat0[, c("Habitat", "Habitat2", "Habitat3")])
  utemp <- utemp[do.call(order, list(utemp[, "Habitat3"], utemp[, "Habitat2"], utemp[, "Habitat"])), ]

  ftemp0 <- file.path(dir_res0, "Tables")
  dir.create(ftemp0, recursive = TRUE, showWarnings = FALSE)

  write.csv(utemp, file = file.path(ftemp0,
    "Table_Habitat_Input_vs_Classification.csv"))
}



#--------------------------------------
#--- Table Data

if (do_Data1_TablesS1and2) {
  responses_Hall <- c("mA", "He", "Sh", "PLp", "Fis")
  irows <- datl[, "Response"] %in% responses_Hall
  icols <- c(
      # Identification
      "ID_effect", "ID_article", "Paper_Nr", "ID_unit", "FirstAuthor",
      # Fragment data
      "is_fragment", "is_standardized", "Fragment_area_ha", "Fragment_area_rank",
      "Fragment_area_source", "Pop_size_n", "Pop_size_rank", "Pop_size_source",
      # Response
      "Response", "Value", "Ind_sampled_nPERfragment",
      # Moderators:
      # Study project
      "Pub_Year", "Study_type", "Marker", "Marker_cat",
      "Effect_longitude", "Effect_latitude",
      "EffectAgeCat2_ord_yrs",
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
    )

  dat_share <- datl[irows, icols]

  # Data S1 Table S2
  write.csv(dat_share, file = file.path(dir_res0, "Tables", "DataS1_TableS2_PublishedData.csv"),
    row.names = FALSE)

  # Metadata S1
  temp <- data.frame(Field = colnames(dat_share), Explanation = NA)
  write.csv(temp, file = file.path(dir_res0, "Tables", "MetadataS1_Table_PublishedData_ColumnExplanation.csv"),
    row.names = FALSE)


  ids <- match(dat_share[, "ID_article"], dat_citation[, "RecordID"], nomatch = 0)
  stopifnot(ids > 0, dat_citation[ids, "Year"] == dat_share[, "Pub_Year"])

  ids <- sort(unique(ids))
  temp <- apply(dat_citation[ids, c("Authors", "Year", "Title", "Journal", "Volume",
    "Pages")], 1, function(x) paste0(x["Authors"], ". ", x["Year"], ". ", x["Title"],
      ". ", x["Journal"], " ", x["Volume"], ":", x["Pages"], "."))
  dat_references <- data.frame(ID_article = dat_citation[ids, "RecordID"], Reference = temp)

  # Data S1 Table S1
  ftemp0 <- file.path(dir_res0, "Tables")
  dir.create(ftemp0, recursive = TRUE, showWarnings = FALSE)

  write.csv(dat_references, file = file.path(ftemp0,
    "DataS1_TableS1_PublishedDataReferences.csv"), row.names = FALSE)
}


# Data S1 Table S3 (do_targets: Effect sizes of main manuscript)

if (do_Data1_TableS3) {
  # See `std_design` in file `Step3_0_Analysis_settings.R` and
  # do_targets` in file `Step5c_Manuscript_Figure_MainForest.R`
  responses <- std_design[["s3_responses"]]
  moderators <-  c(moderators_H12, moderators_Htime)
  interaction_wHabitat3 <- TRUE
  withControlsL. <- std_design[["s1_wcontr"]]
  only_useadj_standardized <- FALSE
  withNonStandardizedL. <- std_design[["s2_wnonnorm"]]
  fragment_sizes. <- std_design[["s4_fragsize"]]
  cor_methods. <- std_design[["s5_cormethod"]]
  cor_transforms. <- std_design[["s6_cortransform"]]
  weight_methods. <- std_design[["d7_weightmethod"]]


  # Compile data used for analysis back into one table
  mods_used_for_effect <- paste0(moderators, "_inUse")
  vars_table1 <- c("ID_effect", "Species_resolved", "Organism_group4",
    "E", weight_methods., "nf", "ni", "FirstAuthor", "Pub_Year", "Marker")
  vars_table <- c(vars_table1, "Response", mods_used_for_effect)

  table_S1S3 <- data.frame(matrix(NA, nrow = 0, ncol = length(vars_table),
    dimnames = list(NULL, vars_table)))

  for (ir in responses) {
    dtag <- paste0("Response_", ir)
    dir_data <- file.path(dir_res_, "Data", dtag)

    temp_withNonStandardizedL <- if (only_useadj_standardized || ir %in% c("Ar", "mA")) {
      # always use only standardized values for Ar and mA
      FALSE
    } else {
      withNonStandardizedL.
    }

    for (igi in seq_along(moderators)) {
      # File identification
      ftag <- filetag_ID(interaction_wHabitat3, fragment_sizes.,
        withControlsL., only_wo_controls = FALSE, temp_withNonStandardizedL,
        cor_methods., cor_transforms., weight_methods.)
      ftag <- paste(ir, moderators[igi], "by", ftag, sep = "_")

      fdata <- file.path(dir_data, paste0(ftag, "_data.rds"))
      mdata <- readRDS(fdata)

      temp <- mdata[["all"]][["d"]]

      has_vars <- vars_table1 %in% colnames(temp)
      if (any(!has_vars)) {
        stop(paste(vars_table1[!has_vars], collapse = ", "))
      }

      out <- temp[, vars_table1]
      out <- cbind(out, matrix(NA, nrow = nrow(out),
        ncol = 1 + length(mods_used_for_effect),
        dimnames = list(NULL, c("Response", mods_used_for_effect))))
      out[, mods_used_for_effect[igi]] <- TRUE
      out[, "Response"] <- ir

      table_S1S3 <- rbind(table_S1S3, out[, vars_table])
    }
  }

  # Convert factors to character strings
  for (k in seq_len(ncol(table_S1S3))) {
    if (is.factor(table_S1S3[, k])) {
      table_S1S3[, k] <- as.character(table_S1S3[, k])
    }
  }


  # Aggregate data for each effect x response
  other_colnames <- !(colnames(table_S1S3) %in% mods_used_for_effect)
  temp <- by(table_S1S3,
    INDICES = as.list(table_S1S3[, c("ID_effect", "Response")]), function(x) {
      h <- x[, other_colnames]
      hu <- unique(h)
      stopifnot(nrow(hu) == 1)
      a <- x[, mods_used_for_effect]
      ag <- as.integer(apply(a, 2, any, na.rm = TRUE))
      cbind(hu, matrix(ag, nrow = 1))
    })
  temp <- do.call("rbind", temp)
  colnames(temp) <- c(colnames(table_S1S3)[other_colnames], moderators)


  table_S1S3_final <- data.frame(temp, stringsAsFactors = FALSE)

  # Make sure that effects were included consistently per species and moderator
  temp <- by(table_S1S3_final[, moderators],
    INDICES = table_S1S3_final[, "Species_resolved"], function(x) {
      nrow(unique(x))
    })

  ibad <- !(temp == 1)
  if (any(ibad)) {
    print(temp[ibad])
    stop("Effects not consistent.")
  }

  # Beautify table
  table_S1S3_final[, "E"] <- round(table_S1S3_final[, "E"], 3)
  table_S1S3_final[, "whierarch"] <- round(table_S1S3_final[, "whierarch"], 3)
  table_S1S3_final[, "ni"] <- round(table_S1S3_final[, "ni"], 1)
  itemp <- "He+Sh" %in% table_S1S3_final[, "Response"]
  table_S1S3_final[itemp, "Response"] <- "He or Sh"

  ctemp <- colnames(table_S1S3_final)
  ctemp[ctemp == "Species_resolved"] <- "Species"
  ctemp[ctemp == "Organism_group4"] <- "Organism group"
  ctemp[ctemp == "Organism_group4.1"] <- "Organism_group4"
  ctemp[ctemp == "E"] <- "Effect size (z of Pearson cor)"
  ctemp[ctemp == "whierarch"] <- "Hierarchical weight"
  ctemp[ctemp == "nf"] <- "Number of fragments"
  ctemp[ctemp == "ni"] <- "Mean number of individuals per fragment"
  colnames(table_S1S3_final) <- ctemp

  ctemp <- colnames(table_S1S3_final)
  ctemp <- ctemp["ID_effect" != ctemp] # remove `ID_effect`
  # re-order columns:
  new_cols <- c("Organism group", "Species", "Response", "Marker")
  ctemp <- c(new_cols, ctemp[!(ctemp %in% new_cols)])

  ids <- do.call("order", table_S1S3_final[ctemp])
  table_S1S3_final <- table_S1S3_final[ids, ctemp]


  # Write to disk file
  ftemp0 <- file.path(dir_res0, "Tables")
  dir.create(ftemp0, recursive = TRUE, showWarnings = FALSE)

  write.csv(table_S1S3_final, file = file.path(ftemp0,
    "DataS1_TableS3_CalculatedEffectSizesPerSpeciesAndResponse.csv"),
    row.names = FALSE)

  # Number of unique effects
  dim(table_S1S3_final)[1] # 206
}


#--------------------------------------
# Appendix 5: Table S2 (specifications of sensitivity analysis)

if (do_Appendix5_TableS2) {
  # See `full_design` in file `Step3_0_Analysis_settings.R` and `args_full` in
  # `Step5b2_Manuscript_Sensitivity.R`

  sensitivity_design <- data.frame(
    Methods = rep(names(design_arguments[["args_full"]]), lengths(design_arguments[["args_full"]])),
    Levels = unlist(design_arguments[["args_full"]]),
    row.names = NULL
  )

  ftemp0 <- file.path(dir_res0, "Tables")
  dir.create(ftemp0, recursive = TRUE, showWarnings = FALSE)

  write.csv(sensitivity_design, file = file.path(ftemp0,
    "AppendixS5_TableS2_Sensitivity_AnalysisDesign.csv"),
    row.names = FALSE)
}


#--------------------------------------
# Appendix 5: Table S3 (relative influence of analysis decisions)

if (do_Appendix5_TableS3) {
  #-- Determine hypothesis support for each sensitivity component as
  # marginal mean for each response and method decision and averaged across
  # moderators, moderator-levels, and model types
  dir_temp <- file.path(dir_res0, "Tables")
  dir.create(dir_temp, recursive = TRUE, showWarnings = FALSE)

  #--- Inputs
  dir_res <- dir_res_
  dir_sens <- file.path(dir_res, "Datasets_Sensitivity")
  do_targets <- TRUE
  args_full <- design_arguments[["args_full"]]
  responses <- c(responses_H1, responses_H2)
  moderators <- c(moderators_Hs, moderators_Htime)
  interaction_wHabitat3 <- TRUE
  temp_xHabitat3 <- if (interaction_wHabitat3) "modXhabitat" else ""

  #---
  var_design <- names(which(lengths(args_full) > 1))
  var_out <- c("Hypothesis_Support", "Nest", "Ntot")

  temp <- args_full[var_design]
  temp1 <- expand.grid(
    method_levels = paste0(rep(names(temp), lengths(temp)), "__", unlist(temp)),
    responses = responses,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  temp2 <- matrix(0, nrow = nrow(temp1), ncol = length(var_out),
    dimnames = list(NULL, var_out))
  marginal_support <- data.frame(temp1, temp2, stringsAsFactors = FALSE)

  temp0 <- strsplit(marginal_support[, "method_levels"], split = "__")
  marginal_support[, "method"] <- clean_labels(sapply(temp0, function(x) x[[1]]))
  marginal_support[, "levels"] <- clean_labels(sapply(temp0, function(x) x[[2]]))

  sided <- lab_responses[match(responses, lab_responses[, "code"]), "hypothesis"]
  hcols <- sapply(sided, function(x)
    switch(EXPR = x, left = "CI_lt0", right = "CI_gt0", NA))

  for (ir in seq_along(responses)) {
    rname <- if (responses[ir] %in% responses_H1) {
        "resp-gendiv"
      } else if (responses[ir] %in% responses_H2) {
        "resp-inbred"
      } else stop()

    for (ig in seq_along(moderators)) {
      # read sensitivity data
      ftag <- if (moderators[ig] %in% moderators_Hs) {
          paste("1-2_Diversity&Inbreeding",
            if (do_targets) "target" else "all", rname, sep = "-")
        } else if (moderators[ig] %in% moderators_Htime) {
          paste("3_TimeEffects",
            if (do_targets) "target" else "all", rname, sep = "-")
        } else stop()

      tag_fix <- if (interaction_wHabitat3) {
          paste0(ftag, "_", temp_xHabitat3)
        } else {
          ftag
        }

      ftemp <- file.path(dir_sens, paste0("Data_Sensitivity_", tag_fix, ".rds"))
      sens <- readRDS(ftemp)

      # determine indices of moderator
      modcats <- dimnames(sens[["out_res"]])[[2]]
      imods <- if (interaction_wHabitat3) {
          grep(ig, as.integer(substr(modcats, 4, 4)))
        } else {
          grep(moderators[ig], modcats)
        }

      # loop over sensitivity design
      for (k1 in seq_along(var_design)) {
        levels <- args_full[[var_design[k1]]]

        for (k2 in seq_along(levels)) {
          idesign <- which(sens[["design"]][, var_design[k1]] == levels[k2])

          temp_margin <- sens[["out_res"]][responses[ir], imods, idesign, , hcols[ir], drop = FALSE]

          support <- 0
          for (k3 in seq_along(idesign)) {
            tempm <- array(temp_margin[, , k3, , ], dim = dim(temp_margin)[-3],
              dimnames = dimnames(temp_margin)[-3])

            # summed marginal support of hypothesis
            support <- support + apply(tempm, 4, sum, na.rm = TRUE)
          }

          # save in output data.frame
          irow <- marginal_support[, "responses"] == responses[ir] &
            grepl(var_design[k1], marginal_support[, "method_levels"]) &
            grepl(levels[k2], marginal_support[, "method_levels"])
          stopifnot(sum(irow) == 1)

          marginal_support[irow, c("Hypothesis_Support", "Nest", "Ntot")] <-
            marginal_support[irow, c("Hypothesis_Support", "Nest", "Ntot")] +
            c(support, sum(!is.na(temp_margin)), prod(dim(temp_margin)[-5]))
        }
      }
    }
  }

  # Reshape to semi-wide format
  temp <- reshape2::melt(marginal_support)
  out <- reshape2::dcast(temp, method + levels ~ responses + variable)

  ftemp <- file.path(dir_temp, "AppendixS5_TableS3_Sensitivity_MarginalSupport.csv")
  write.csv(out, file = ftemp, row.names = FALSE)


  # Visualize table
  ids <- c(which(out[-1, "method"] != out[-nrow(out), "method"]), nrow(out))
  out2 <- out[0, ]
  addNAs <- out[1, ]
  addNAs[1, ] <- NA

  k0 <- 1
  for (k1 in ids) {
    out2 <- rbind(out2, out[k0:k1, ], addNAs)
    k0 <- k1 + 1
  }
  out2 <- out2[-nrow(out2), ]

  rn <- as.character(apply(out2[, c("method", "levels")], 1, paste0, collapse = ": "))
  rn[rn %in% "NA: NA"] <- NA
  rn <- gsub("-published", "", rn)
  rn <- gsub(".:", ":", rn)
  rn <- gsub("Cor Methods: Kendall", "Kendall's correlation coefficient r", rn)
  rn <- gsub("Cor Methods: Pearson", "Pearson's r*", rn)
  rn <- gsub("Cor Methods: Spearman", "Spearman's r", rn)
  rn <- gsub("Cor Transforms: None", "Untransformed r", rn)
  rn <- gsub("Cor Transforms: Unbiased", "Unbiased r", rn)
  rn <- gsub("Cor Transforms: Ztransform", "z-transformed r*", rn)
  rn <- gsub("Fragment Sizes: Fragment Area Ha", "Fragment size: area (ha)*", rn)
  rn <- gsub("Fragment Sizes: Fragment Area Rank", "--: area rank (ordinal)", rn)
  rn <- gsub("Fragment Sizes: Pop Size N", "--: population size (n)", rn)
  rn <- gsub("Fragment Sizes: Pop Size Rank", "--: population rank (ordinal)", rn)
  rn <- gsub("Weight Methods: Whierarch", "Hierarchical weights*", rn)
  rn <- gsub("Weight Methods: Winvar", "Inverse-variance weights", rn)
  rn <- gsub("WithControlsL: FALSE", "Exclude control areas*", rn)
  rn <- gsub("WithControlsL: TRUE", "Include controls as large fragments", rn)

  out_rel <- as.matrix(data.frame(
    A = out2[, "mA_Hypothesis_Support"] / out2[, "mA_Ntot"],
    He = out2[, "He+Sh_Hypothesis_Support"] / out2[, "He+Sh_Ntot"],
    PLP = out2[, "PLp_Hypothesis_Support"] / out2[, "PLp_Ntot"],
    Fis = out2[, "Fis_Hypothesis_Support"] / out2[, "Fis_Ntot"]
  ))

  ftemp <- file.path(dir_temp, "Fig4_Sensitivity_MarginalSupport.pdf")
  pdf(height = 4, width = 6, file = ftemp)
  par_prev <- par(mar = c(3, 15, 0.5, 0.75), mgp = c(1, 0, 0), tcl = 0.3,
    cex = 0.75)
  pchs <- seq_len(4)
  cols <- viridisLite::viridis(4 + 1)

  ids <- nrow(out_rel):1
  matplot(x = out_rel, y = ids, type = "l", lty = 1, col = cols,
    xlab = "Relative marginal support", ylab = "", axes = FALSE)
  matpoints(x = out_rel, y = ids, pch = pchs, cex = 1.5, col = cols)
  axis(1)
  segments(0, 0, 0, nrow(out_rel))
  text(x = -0.025, y = ids, labels = rn, adj = 1, xpd = NA)
  legend("topright", legend = colnames(out_rel), pch = pchs, col = cols,
    pt.cex = 1.5, lwd = 2)

  par(par_prev)
  dev.off()
}
