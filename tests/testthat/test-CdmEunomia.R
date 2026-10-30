# The production path end to end on Eunomia (a small public CDM; downloaded on first use): cohort generation from the
# shipped SQL, light counts for all 26 original tasks (most have too few events there, which must not stop the run),
# and the export.

test_that("counts run end to end on a CDM (Eunomia)", {
  skip_if_not_installed("Eunomia")
  skip_on_cran()
  details <- tryCatch(Eunomia::getEunomiaConnectionDetails(databaseFile = tempfile(fileext = ".duckdb"), dbms = "duckdb"),
                      error = function(e) skip(paste("Eunomia not available:", conditionMessage(e))))
  outputFolder <- tempfile("sar_eunomia_")
  on.exit(unlink(outputFolder, recursive = TRUE), add = TRUE)
  dataSource <- createCdmDataSource(details, cdmDatabaseSchema = "main", cohortDatabaseSchema = "main",
                                    databaseName = "EUNOMIA")
  tasks <- studyTasks()
  counts <- suppressWarnings(execute(dataSource, tasks, outputFolder, mode = "counts"))
  errors <- list.files(file.path(outputFolder, "counts", "EUNOMIA"), "_error\\.txt$")
  expect_identical(length(errors) + NROW(counts), nrow(tasks)) # every task is either counted or has its error recorded
  if (NROW(counts) > 0) expect_true(all(counts$n == counts$poolN + counts$testN))
  expect_false(any(file.exists(file.path(outputFolder, "data")))) # counts extract no covariates and cache nothing
  zipFile <- exportResults(outputFolder)
  expect_true("counts.csv" %in% utils::unzip(zipFile, list = TRUE)$Name || NROW(counts) == 0)
})
