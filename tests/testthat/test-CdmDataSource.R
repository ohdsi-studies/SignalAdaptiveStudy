# The production data source: cohort definitions shipped with the package, task mapping, and the original task rules
# as PLP population settings (checked at the day boundaries). Against the original SQL on a real CDM:
# extras/production/validateCdmSource.R.

test_that("the cohort definitions cover every task", {
  cohorts <- studyCohortDefinitionSet()
  expect_identical(nrow(cohorts), 26L)
  expect_false(anyDuplicated(cohorts$cohortId) > 0)
  expect_true(all(nzchar(cohorts$json) & grepl("@target_cohort_id", cohorts$sql)))
  tasks <- originalCohortTasks()
  expect_identical(nrow(tasks), 26L)
  expect_true(all(tasks$targetId %in% cohorts$cohortId[cohorts$role == "target"]))
  expect_true(all(tasks$outcomeId %in% cohorts$cohortId[cohorts$role == "outcome"]))
  expect_identical(nrow(originalCohortTasks("IPCI")), 8L)
  catalogue <- catalogueCohortTasks()
  expect_setequal(catalogue$catalogueTask, c("SAR-017", "SAR-175", "SAR-287", "SAR-288"))
  expect_true(all(catalogue$targetId %in% cohorts$cohortId[cohorts$role == "target"]))
  expect_true(all(catalogue$outcomeId %in% cohorts$cohortId[cohorts$role == "outcome"]))
  all <- studyTasks()
  expect_identical(nrow(all), 30L)
  expect_false(anyDuplicated(all$taskId) > 0)
  expect_identical(taskPopulationSettings(all[all$taskId == "SAR-287_MDD_to_bipolar_1y", ])$priorOutcomeLookback, 999999)
  expect_identical(taskPopulationSettings(all[1, ])$priorOutcomeLookback, 366)
})

test_that("population settings follow the original rules at the day boundaries", {
  horizon <- 365
  # one person per case: (days of the outcome relative to index, days of observation after index)
  cases <- data.frame(
    outcomeDay = c(NA, -366, -365, 0, 1, horizon, horizon + 1, NA, 10),
    obsEnd = c(horizon, horizon, horizon, horizon, horizon, horizon, horizon, horizon - 1, horizon - 1),
    kept = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE),
    outcome = c(0, 0, NA, NA, 1, 1, 0, NA, NA))
  n <- nrow(cases)
  cohorts <- data.frame(rowId = seq_len(n), subjectId = seq_len(n), targetId = 1,
                        cohortStartDate = as.Date("2015-01-01"), daysFromObsStart = 400, daysToCohortEnd = 0,
                        daysToObsEnd = cases$obsEnd, ageYear = 50, gender = 8507)
  attr(cohorts, "metaData") <- list(targetId = 1, attrition = data.frame(outcomeId = 2, description = "test",
                                                                          targetCount = n, uniquePeople = n, outcomes = 1))
  has <- !is.na(cases$outcomeDay)
  outcomes <- data.frame(rowId = which(has), outcomeId = 2, outcomeCount = 1, daysToEvent = cases$outcomeDay[has])
  attr(outcomes, "metaData") <- data.frame(outcomeIds = 2)
  covariateData <- Andromeda::andromeda(
    covariates = data.frame(rowId = 1, covariateId = 1, covariateValue = 1),
    covariateRef = data.frame(covariateId = 1, covariateName = "x", analysisId = 1, conceptId = 0),
    analysisRef = data.frame(analysisId = 1, analysisName = "x", domainId = "x", startDay = 0, endDay = 0,
                             isBinary = "Y", missingMeansZero = "Y"))
  class(covariateData) <- "CovariateData"
  plpData <- structure(list(cohorts = cohorts, outcomes = outcomes, covariateData = covariateData,
                            metaData = list(databaseDetails = list(), restrictPlpDataSettings = list())),
                       class = "plpData")
  population <- suppressMessages(PatientLevelPrediction::createStudyPopulation(
    plpData, 2, originalProductionPopulationSettings(horizon)))
  expect_setequal(population$rowId, which(cases$kept))
  expect_identical(as.numeric(population$outcomeCount[order(population$rowId)] > 0), cases$outcome[cases$kept])
})

test_that("large tasks are sampled for the model runs, enough for the largest budget", {
  tasks <- data.frame(taskId = c("small", "large", "rare", "huge_but_rare"), targetId = 1, outcomeId = 2)
  counts <- data.frame(taskId = tasks$taskId, n = c(5e5, 2.8e7, 2.8e7, 3e6), poolN = c(3.75e5, 2.1e7, 2.1e7, 2.25e6),
                       poolEvents = c(400, 3e5, 1.69e4, 900))
  s <- addSampleSizes(tasks, counts)
  expect_true(is.na(s$sampleSize[s$taskId == "small"]))                 # under the cap: not sampled
  expect_identical(s$sampleSize[s$taskId == "large"], 2e6)              # plenty of events at the cap
  rare <- s[s$taskId == "rare", ]
  expect_gt(rare$sampleSize, 2e6)                                      # the cap would give < 1.5 x 1000 events
  expect_gte(rare$expectedPoolEvents, 1500)
  expect_true(is.na(s$sampleSize[s$taskId == "huge_but_rare"]))         # needs more than the whole population
  expect_identical(s$budgets[[which(s$taskId == "small")]], c(100, 250))  # budgets up to the pool events
  expect_identical(s$budgets[[which(s$taskId == "large")]], c(100, 250, 500, 750, 1000))
  expect_equal(s$sampleFraction[s$taskId == "large"], 2e6 / 2.8e7)       # persons sampled: population x fraction
  expect_true(is.na(s$sampleFraction[s$taskId == "small"]))
})

test_that("the person hash matches exact integer arithmetic, in R and in SQL on duckdb and SQLite", {
  ids <- c(1, 2, 3, 41, 123, 2147483646, 2147483647, 2147483648, 9876543210123, -5)
  expected <- c(1347478916, 527213154, 1854431039, 749103141, 59303714, 840526793, 20261031, 1347478916, 1251708369,
                213899515) # Python big integers: ((|id| mod m) * a + seed) mod m
  expect_identical(personHash(ids), expected)
  skip_if_not_installed("DatabaseConnector")
  for (dbms in c("duckdb", "sqlite")) {
    if (dbms == "duckdb") skip_if_not_installed("duckdb")
    if (dbms == "sqlite") skip_if_not_installed("RSQLite")
    connection <- DatabaseConnector::connect(dbms = dbms, server = tempfile())
    sql <- SqlRender::translate(sprintf("SELECT %s AS h FROM (SELECT CAST(%.0f AS BIGINT) AS subject_id) t;",
                                        personHashSql("subject_id", dbms), 9876543210123), targetDialect = dbms)
    expect_equal(DatabaseConnector::querySql(connection, sql)[[1]], expected[9], info = dbms)
    DatabaseConnector::disconnect(connection)
  }
})

test_that("the target sample keeps whole persons, about the fraction, as the hash says", {
  skip_if_not_installed("DatabaseConnector")
  skip_if_not_installed("duckdb")
  connection <- DatabaseConnector::connect(dbms = "duckdb", server = tempfile())
  on.exit(DatabaseConnector::disconnect(connection), add = TRUE)
  set.seed(5)
  ids <- sample.int(1e9, 20000)
  cohort <- data.frame(cohort_definition_id = c(rep(1, 20000), rep(1, 5000), rep(2, 100)),
                       subject_id = c(ids, ids[1:5000], ids[1:100]),
                       cohort_start_date = as.Date("2015-01-01") + c(seq_len(20000), seq_len(5000) + 30000, seq_len(100)),
                       cohort_end_date = as.Date("2030-01-01"))
  DatabaseConnector::insertTable(connection, "main", "cohort", cohort, createTable = TRUE, progressBar = FALSE)
  createTargetSample(connection, "main", "cohort", "cohort_sample", targetId = 1, fraction = 0.1)
  s <- DatabaseConnector::querySql(connection, "SELECT * FROM main.cohort_sample;", snakeCaseToCamelCase = TRUE)
  expect_true(all(s$cohortDefinitionId == 1))
  kept <- unique(ids[personHash(ids) < sampleThreshold(0.1)])
  expect_setequal(unique(s$subjectId), kept)
  expect_equal(length(kept) / 20000, 0.1, tolerance = 0.1)
  expect_identical(nrow(s), sum(cohort$cohort_definition_id == 1 & cohort$subject_id %in% kept)) # every entry of a kept person
})

test_that("the hash SQL translates to every dialect without a bare modulo where it is not supported", {
  for (dbms in SqlRender::listSupportedDialects()$dialect) {
    sql <- SqlRender::translate(sprintf("SELECT %s FROM t;", personHashSql("subject_id", dbms)), targetDialect = dbms)
    if (dbms %in% c("oracle", "iris", "bigquery")) expect_false(grepl("%", sql, fixed = TRUE), info = dbms)
    if (dbms %in% c("oracle", "iris", "bigquery")) expect_true(grepl("MOD(", toupper(sql), fixed = TRUE), info = dbms)
  }
})
