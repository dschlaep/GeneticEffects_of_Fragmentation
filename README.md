# Genetic effects of anthropogenic habitat fragmentation on remnant animal and plant populations: A meta-analysis.

Daniel R. Schlaepfer, Brigitte Braschler, Hans-Peter Rusterholz, and Bruno Baur

Journal: Ecosphere

Data available from Dryad:
- TABLE S1. References of studies used in meta-analysis
- TABLE S2. Prepared set of extracted data. Data file is in ‘long-format’, i.e., each row contains data for one combination of fragment, genetic measure (A, He, Sh, PLP, and Fis), species, genetic marker, and study
- METADATA S2. Explanation of columns for TABLE S2
- Table S3. Calculated effect sizes for main analysis

Code available from https://github.com/dschlaep/GeneticEffects_of_Fragmentation

Description of how the files are to be used:
- Data S1: Table S2 is a `csv`-representation of underlying data used for the analysis. It was produced by the code file `Step5g_Manuscript_SupplementaryTables.R` with the logical flag `do_Data1_TablesS1and2`. Thus, this table is produced by the variable `dat0_final` as output of the code file `Step1c_Moderators.R` and is loaded in the subsequent steps as `datl`.
- The code file `0_Run_all_analyses.R` is the master file which documents and executes all code; it also identifies the figures and tables of the manuscripts. There are six steps: (i) data preparation, (ii) calculate effect sizes, (iii) fit meta-analysis models, (iv) produce miscellaneous output figures, (v) prepare data for output, and (vi) produce main results, figures, and tables. Each of these steps is controlled by a logical flag in the `0_Run_all_analyses.R` and executes a subset of the other code files.
