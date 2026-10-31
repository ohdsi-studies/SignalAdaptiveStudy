# Site settings for extras/CodeToRun.R. Copy this file to extras/localSettings.R (git-ignored) and fill it in.
# Credentials do not go here: they come from environment variables (see CodeToRun.R).

scratchSchema <- "scratch"         # a schema you can write to (e.g. catalog.schema); cohort tables and temp tables
minCellCount <- 5                  # your site's minimum cell count

databaseName <- "CCAE"             # the database of this run: one of names(cdmSchemas)

# The CDM (and vocabulary) schema of each database, by the short name used in the results
cdmSchemas <- c(
  CCAE = "cdm_ccae",
  MDCD = "cdm_mdcd",
  MDCR = "cdm_mdcr",
  OPTUM_EHR = "cdm_optum_ehr",
  OPTUM_CLAIMS = "cdm_optum_claims"
)

# Worker processes for the models (step 3c): one R process each, one thread each. About the number of cores minus
# one or two; memory is shared out among them automatically. 1: everything in the R session of CodeToRun.R.
workers <- 6

# Optional: run only these tasks in step 3c (task IDs as in studyTasks()); leave out to run all
# onlyTasks <- c("heart_failure_to_COPD_3y", "CKD_to_heart_failure_1y")
