# Unit tests for TaskDiagnosticsV3 (Modules 1-7), privacy censoring, removePriorOutcome,
# restartability, and multi-dialect SqlRender translation.

syntheticV3Cdm <- function(dbms = "duckdb") {
  file <- tempfile(fileext = paste0(".", dbms))
  details <- DatabaseConnector::createConnectionDetails(dbms = dbms, server = file)
  connection <- DatabaseConnector::connect(details)
  on.exit(DatabaseConnector::disconnect(connection))
  d <- function(x) as.Date(x)
  put <- function(name, data) {
    DatabaseConnector::insertTable(connection, databaseSchema = "main", tableName = name, data = data,
                                   createTable = TRUE, dropTableIfExists = TRUE, progressBar = FALSE,
                                   camelCaseToSnakeCase = FALSE)
  }
  n <- 60L
  ids <- seq_len(n)
  # Two care sites: site 1 has events, site 2 has none (tests M5 zero-sites & syAUC)
  put("person", data.frame(
    person_id = ids,
    year_of_birth = rep(c(1950L, 1960L), length.out = n),
    gender_concept_id = ifelse(ids %% 2 == 0, 8532L, 8507L),
    care_site_id = ifelse(ids <= 40L, 1L, 2L)
  ))
  opEnd <- rep("2020-12-31", n)
  opEnd[31:35] <- "2010-03-01" # true LTFU
  opEnd[36:40] <- "2020-11-15" # within 90d of dbLastDate (admin censored if indexed late)
  idxDates <- rep("2010-01-01", n)
  idxDates[36:40] <- "2020-06-01"
  put("observation_period", data.frame(
    person_id = ids,
    observation_period_start_date = d(rep("2000-01-01", n)),
    observation_period_end_date = d(opEnd)
  ))
  put("concept", data.frame(
    concept_id = c(100L, 101L, 200L, 201L, 300L, 301L, 400L, 9201L, 8753L),
    concept_name = c("Target Condition", "Orphan Target Condition", "Outcome Condition",
                     "Drug Oral Tablet", "Drug Clinical Drug Form", "ATC4 Class",
                     "Lab Potassium", "Inpatient Visit", "millimole per liter"),
    domain_id = c("Condition", "Condition", "Condition", "Drug", "Drug", "Drug", "Measurement", "Visit", "Unit"),
    vocabulary_id = c("SNOMED", "SNOMED", "SNOMED", "RxNorm", "RxNorm", "ATC", "LOINC", "Visit", "UCUM"),
    concept_class_id = c("Clinical Finding", "Clinical Finding", "Clinical Finding",
                         "Clinical Drug", "Clinical Drug Form", "ATC 4th", "Lab Test", "Visit", "Unit"),
    standard_concept = c("S", "S", "S", "S", "S", "C", "S", "S", "S"),
    stringsAsFactors = FALSE
  ))
  put("concept_ancestor", data.frame(
    ancestor_concept_id = c(100L, 101L, 200L, 201L, 300L, 300L, 301L, 400L),
    descendant_concept_id = c(100L, 101L, 200L, 201L, 300L, 201L, 201L, 400L)
  ))
  # Conditions: all 60 get target (100) on idxDate; persons 1..30 get outcome (200) on 2010-01-05 (d1-7) or 2010-06-01;
  # persons 1..15 also get orphan concept (101)
  evDates <- c(rep("2010-01-05", 10), rep("2010-06-01", 20))
  conds <- rbind(
    data.frame(person_id = ids, condition_concept_id = 100L, condition_start_date = d(idxDates)),
    data.frame(person_id = 1:30, condition_concept_id = 200L, condition_start_date = d(evDates)),
    data.frame(person_id = 1:15, condition_concept_id = 101L, condition_start_date = d(rep("2010-01-01", 15)))
  )
  put("condition_occurrence", conds)
  # Drugs: persons 1..30 get drug 201 on index & event dates
  put("drug_exposure", rbind(
    data.frame(person_id = 1:30, drug_concept_id = 201L, drug_exposure_start_date = d(idxDates[1:30])),
    data.frame(person_id = 1:30, drug_concept_id = 201L, drug_exposure_start_date = d(evDates))
  ))
  put("procedure_occurrence", data.frame(person_id = 1L, procedure_concept_id = 100L, procedure_date = d("2010-01-01"))[0, ])
  put("measurement", data.frame(
    person_id = 1:10, measurement_concept_id = 400L, measurement_date = d(idxDates[1:10]),
    unit_concept_id = 8753L, value_as_number = seq(4.1, 5.0, by = 0.1)
  ))
  put("visit_occurrence", data.frame(
    person_id = 1:30, visit_concept_id = 9201L,
    visit_start_date = d(evDates), visit_end_date = d(evDates) + 2L
  ))
  # Duplicate death row for person 45 to verify deduplication in #v3_ev
  put("death", data.frame(person_id = c(45L, 45L), death_date = d(c("2010-04-01", "2010-04-02"))))

  # Target cohort (1): all 60 persons; Outcome cohort (2): persons 1..30 in TAR + persons 1..5 ALSO have prior outcome in 2009
  cohort <- rbind(
    data.frame(cohort_definition_id = 1L, subject_id = ids, cohort_start_date = d(idxDates), cohort_end_date = d(idxDates)),
    data.frame(cohort_definition_id = 2L, subject_id = 1:30, cohort_start_date = d(evDates), cohort_end_date = d(evDates)),
    data.frame(cohort_definition_id = 2L, subject_id = 1:5, cohort_start_date = d(rep("2009-01-01", 5)), cohort_end_date = d(rep("2009-01-01", 5)))
  )
  put("cohort", cohort)
  details
}

test_that("TaskDiagnosticsV3 runs all 7 modules, enforces minCellCount, deduplicates death, and honors removePriorOutcome", {
  skip_if_not_installed("duckdb")
  details <- syntheticV3Cdm("duckdb")
  out <- tempfile("sar_v3_")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  tJson <- '{"ConceptSets":[{"id":0,"name":"t","expression":{"items":[{"concept":{"CONCEPT_ID":100,"CONCEPT_NAME":"Target Condition"},"includeDescendants":true,"isExcluded":false}]}},{"id":1,"name":"d","expression":{"items":[{"concept":{"CONCEPT_ID":201,"CONCEPT_NAME":"Drug Oral Tablet"},"includeDescendants":true,"isExcluded":false}]}},{"id":2,"name":"m","expression":{"items":[{"concept":{"CONCEPT_ID":400,"CONCEPT_NAME":"Lab Potassium"},"includeDescendants":true,"isExcluded":false}]}}],"PrimaryCriteria":{"CriteriaList":[{"ConditionOccurrence":{"CodesetId":0}},{"DrugExposure":{"CodesetId":1}},{"Measurement":{"CodesetId":2}}]}}'
  oJson <- '{"ConceptSets":[{"id":0,"name":"o","expression":{"items":[{"concept":{"CONCEPT_ID":200,"CONCEPT_NAME":"Outcome Condition"},"includeDescendants":true,"isExcluded":false}]}}],"PrimaryCriteria":{"CriteriaList":[{"ConditionOccurrence":{"CodesetId":0}}]}}'
  defs <- data.frame(cohortId = c(1L, 2L), json = c(tJson, oJson), sql = c("t", "o"), stringsAsFactors = FALSE)
  tasks <- data.frame(
    taskId = c("rm_prior_true", "rm_prior_false"),
    targetId = c(1L, 1L), outcomeId = c(2L, 2L),
    riskWindowStart = c(1L, 1L), riskWindowEnd = c(365L, 365L),
    removePriorOutcome = c(TRUE, FALSE),
    outcomeKind = c("acute", "acute"),
    stringsAsFactors = FALSE
  )

  zipPath <- runTaskDiagnosticsV3(
    connectionDetails = details, cdmDatabaseSchema = "main",
    cohortDatabaseSchema = "main", cohortTable = "cohort",
    tasks = tasks, cohortDefinitionSet = defs,
    databaseName = "SYN", outputFolder = out,
    minCellCount = 5L, tempEmulationSchema = NULL,
    runOrphanScan = TRUE, minOrphanPersons = 10L
  )
  expect_true(file.exists(zipPath))
  v3Dir <- file.path(out, "taskDiagnosticsV3", "SYN")
  for (f in c("summary_v3.csv", "index_codes_v3.csv", "orphan_concepts_v3.csv",
              "temporal_profile_v3.csv", "visit_context_v3.csv", "drug_route_v3.csv", "lab_units_v3.csv")) {
    expect_true(file.exists(file.path(v3Dir, f)), info = f)
  }

  sm <- utils::read.csv(file.path(v3Dir, "summary_v3.csv"), stringsAsFactors = FALSE)
  expect_equal(nrow(sm), 2L)
  # Duplicate death row for person 45 must NOT inflate population beyond 60 / 55
  expect_equal(sm$population[sm$task_id == "rm_prior_true"], 55)
  expect_equal(sm$events[sm$task_id == "rm_prior_true"], 25)
  expect_equal(sm$population[sm$task_id == "rm_prior_false"], 60)
  expect_equal(sm$events[sm$task_id == "rm_prior_false"], 30)
  # Competing death count = 1 (<5) -> pct_competing_death must be blanked (NA) for privacy
  expect_true(is.na(sm$pct_competing_death[1]))

  # M2 orphan scan catches concept 101 ("Orphan Target Condition")
  orph <- utils::read.csv(file.path(v3Dir, "orphan_concepts_v3.csv"), stringsAsFactors = FALSE)
  expect_true(101L %in% orph$concept_id)

  # M4 visit context sees Inpatient Visit
  vc <- utils::read.csv(file.path(v3Dir, "visit_context_v3.csv"), stringsAsFactors = FALSE)
  expect_true("Inpatient Visit" %in% vc$visit_type)

  # M7 drug route & lab units
  dr <- utils::read.csv(file.path(v3Dir, "drug_route_v3.csv"), stringsAsFactors = FALSE)
  expect_true("Drug Clinical Drug Form" %in% dr$clinical_drug_form)
  lu <- utils::read.csv(file.path(v3Dir, "lab_units_v3.csv"), stringsAsFactors = FALSE)
  expect_true(400L %in% lu$concept_id)
  expect_false(any(c("min_val", "max_val") %in% names(lu)))

  # Privacy check: zero positive counts in 1..4 across every V3 CSV
  for (f in list.files(v3Dir, "_v3\\.csv$", full.names = TRUE)) {
    df <- utils::read.csv(f, stringsAsFactors = FALSE)
    cntCols <- intersect(c("population", "events", "n_persons", "n_cases", "any_n", "n_rows", "person_count", "post_obs_deaths_in_tar"), names(df))
    for (col in cntCols) {
      vals <- df[[col]]
      expect_false(any(!is.na(vals) & vals > 0 & vals < 5), info = paste(basename(f), col))
    }
  }

  # Restartability: second call skips cached tasks (mtime unchanged)
  rds1 <- file.path(v3Dir, "tasks", "rm_prior_true.rds")
  stamp <- file.mtime(rds1)
  Sys.sleep(1.1)
  runTaskDiagnosticsV3(
    connectionDetails = details, cdmDatabaseSchema = "main",
    cohortDatabaseSchema = "main", cohortTable = "cohort",
    tasks = tasks, cohortDefinitionSet = defs,
    databaseName = "SYN", outputFolder = out,
    minCellCount = 5L, tempEmulationSchema = NULL,
    runOrphanScan = TRUE, minOrphanPersons = 10L
  )
  expect_identical(file.mtime(rds1), stamp)
})

test_that("TaskDiagnosticsV3 SQL translates cleanly across major OHDSI dialects", {
  dialects <- c("postgresql", "redshift", "spark", "snowflake", "bigquery", "oracle", "sql server", "duckdb")
  templates <- c(
    "SELECT COUNT(*) AS population, SUM(CASE WHEN op_end >= DATEADD(dd, -90, DATEFROMPARTS(@db_last_y, @db_last_m, @db_last_d)) THEN 1 ELSE 0 END) AS admin_censored FROM #v3_ev;",
    "SELECT subject_id, index_date AS event_dt INTO #v3_samp_t FROM (SELECT subject_id, index_date, ROW_NUMBER() OVER (ORDER BY subject_id) AS rn FROM #v3_ev) x WHERE rn <= 25000;",
    "SELECT concept_id, concept_name FROM (SELECT c.concept_id, c.concept_name, ROW_NUMBER() OVER (ORDER BY COUNT(DISTINCT s.subject_id) DESC, c.concept_id) AS rn FROM #v3_samp_o s INNER JOIN @vocab.concept c ON c.concept_id = s.subject_id GROUP BY c.concept_id, c.concept_name) x WHERE rn <= 10;"
  )
  for (dbms in dialects) {
    for (tpl in templates) {
      r <- SqlRender::render(tpl, vocab = "cdm", db_last_y = 2026L, db_last_m = 6L, db_last_d = 30L, warnOnMissingParameters = FALSE)
      tr <- suppressWarnings(SqlRender::translate(r, targetDialect = dbms, tempEmulationSchema = "scratch"))
      expect_true(nzchar(tr), info = dbms)
      expect_false(grepl("@", tr, fixed = TRUE), info = dbms)
    }
  }
})
