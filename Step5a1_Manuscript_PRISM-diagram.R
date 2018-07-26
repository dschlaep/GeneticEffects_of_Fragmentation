##  Moher D, Liberati A, Tetzlaff J, Altman DG, The PRISMA Group (2009). Preferred Reporting
##    Items for Systematic Reviews and Meta-Analyses: The PRISMA Statement.
##    PLoS Med 6(7): e1000097. doi:10.1371/journal.pmed1000097
##
##  Koricheva, J. & Gurevitch, J. (2014). Uses and misuses of meta-analysis in plant
##    ecology. J Ecol, 102, 828-844.


#--- SETTINGS
# Paths
dir_prj1 <- "Prj04_GenetFragment_MetaAnalysis"

dir_ana <- file.path(dir_prj1, "3_Analysis")

#--- Setup
source(file.path(dir_ana, "Step5-0_Manuscript-methods.R"))


# File names
tag <- format(Sys.Date(), "%Y%m%d")
fout <- file.path(dir_out1, paste0(tag, "_PRISM-diagram.rds"))
ffig <- file.path(dir_res0, "Figs_PRISM", paste0("Fig_PRISM-plot_v", tag, ".pdf"))

dir.create(dirname(ffig), recursive = TRUE, showWarnings = FALSE)


#------------------------------------

# Data

set_size <- list(
  identification = list(
    # data in 20160129_CompleteLiteratureSearch.docx
    N_from_DBs = c(
      ISI_WoK = 2299 + 35 + 133,
      Science_Direct = 263+ 33,
      Scopus = 1296 + 108,
      Scielo = 51,
      Agricola_NAL = 1076,
      CiteSeer_x6M = 500,
      ProQuest = 98,
      CAB_Abstracts = 1525,
      DOAJ = 0,
      Google_Scholar = 100
    ),
    # data from References161220_cleaned.tab
    N_from_other = 268
  ),

  screening = list(
    # data from LitDBcomplete_FragmGenet_Cleaned_v3.enlp
    #   groups: v2f_20160217_Reviewers, v3a_Added2016122_New, v3b_Added20170111New
    N_wo_duplicates = N_wo_duplicates <- 3416 + 151 + 158,
    # Title and abstract of all records screened, i.e., N_screened = N_wo_duplicates
    N_screened = N_wo_duplicates,
    #
    N_excluded1 = c(
      Brigitte = 992 - NA,
      Daniel = 992 - 188, # based on 20161102_CompleteDB_AssignedToDaniel.xlsx
      Bruno = 942, # email Feb 6, 2017
      Hans_Peter = 925 # printout Feb 6, 2017
    )
  ),

  eligility = list(
    N_fulltext_assessed = c(
      Brigitte = NA,
      Daniel = 188, # based on 20161102_CompleteDB_AssignedToDaniel.xlsx
      Bruno = 985 - 942, # email Feb 6, 2017
      Hans_Peter = 994 - 925 # printout Feb 6, 2017
    ),
    N_excluded2 = c(
      Brigitte = 992 - NA,
      Daniel = 188 - 94, # based on 20161102_CompleteDB_AssignedToDaniel.xlsx
      Bruno = NA - 40, # email Feb 6, 2017
      Hans_Peter = 33 # printout Feb 6, 2017
    )
  ),

  eligility3 = list(
    N_fulltext_assessed3 = c(
      Bruno = NA, # email Feb 6, 2017
      Hans_Peter = 32 # printout Feb 6, 2017
    ),
    N_excluded3 = c(
      Bruno = NA - 40, # email Feb 6, 2017
      Hans_Peter = 4, #oder 2 # printout Feb 6, 2017
      others = 2 # e.g., Dec 19, 2017: excluded moss/lichen study
    )
  ),

  included = list(
    N_qualitative = NULL,
    N_quantitative = 123
  )
)


#--- Produce a PRISM diagram

library(metagear)

phases <- c(
  paste0("START_PHASE: ", sum(set_size$identification$N_from_DBs), " studies identified ",
    "from database searches:\n",
    paste(unlist(set_size$identification$N_from_DBs), "from",
      gsub("_", " ", names(set_size$identification$N_from_DBs)), collapse = "\n")),
  paste0("START_PHASE: ", set_size$identification$N_from_other,
    " from other sources"),

  paste0(sum(unlist(set_size$identification)), " pooled studies"),
  paste0("EXCLUDE_PHASE: ",
    sum(unlist(set_size$identification)) - set_size$screening$N_wo_duplicates,
    " duplicates excluded"),

  #paste0(set_size$screening$N_screened, " studies with title and abstract screened"),
  #"EXCLUDE_PHASE: # of studies excluded",
  #"# of full-text articles assessed for eligibility",
  #"EXCLUDE_PHASE: # of full-text articles excluded, not fitting eligibility crit "
  # of studies included in qualitative synthesis",
  #"EXCLUDE_PHASE: # studies excluded, incomplete data reported",

  paste0(set_size$screening$N_screened, " studies screened"),
  paste0("EXCLUDE_PHASE: ",
    set_size$screening$N_screened - set_size$included$N_quantitative,
    " studies excluded"),

  paste0(set_size$included$N_quantitative, " studies included in meta-analysis"))



pdf(height = 10, width = 8, file = ffig)
thePlot <- plot_PRISMA(phases, design = c(flatBox = TRUE))
dev.off()
