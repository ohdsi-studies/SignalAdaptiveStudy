# Repairs the test AUCs that step 3c recorded as NA with package versions before 2026-11-01 (integer overflow for test
# sets of ~500,000 rows with thousands of events). Each such run's saved model predicts the test set again; the AUC is
# written into result.csv only when the predictions give the stored log loss and Brier score (the original is kept as
# result_na_auc.csv). No refitting and no database connection: it uses the plpData cache. Then writes the zip.
#
#   Rscript extras/production/repairAuc.R            (in the project folder, after git pull and reinstalling)

library(SignalAdaptiveStudy)
source("extras/localSettings.R")
outputFolder <- file.path(getwd(), "SignalAdaptiveStudy_output", databaseName)
# no connection: the data come from the cache (a cache that does not match its task would need one, and fail)
dataSource <- createCdmDataSource(connectionDetails = NULL, cdmDatabaseSchema = cdmSchemas[[databaseName]],
                                  cohortDatabaseSchema = scratchSchema, databaseName = databaseName,
                                  cohortTable = paste0("sar_study_cohort_", tolower(databaseName)))
counts <- execute(dataSource, studyTasks(), outputFolder, mode = "counts")   # the 3b counts (read back)
modelTasks <- addSampleSizes(originalCohortTasks(), counts[counts$feasible, ])
results <- SignalAdaptiveStudy:::repairTestAuc(dataSource, modelTasks, outputFolder,
                                               workers = if (exists("workers")) workers else 1)
message(sprintf("%d runs, %d without a test AUC", nrow(results), sum(is.na(results$auc))))
exportResults(outputFolder, minCellCount = minCellCount)
