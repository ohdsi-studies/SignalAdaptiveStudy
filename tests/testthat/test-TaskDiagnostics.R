# Task diagnostics on a tiny synthetic CDM where each person is one scenario and every count is known by hand; on duckdb and
# SQLite. The SQL also translates to every dialect.

syntheticCdm <- function(dbms) {
  file <- tempfile(fileext = paste0(".", dbms))
  details <- DatabaseConnector::createConnectionDetails(dbms = dbms, server = file)
  connection <- DatabaseConnector::connect(details)
  on.exit(DatabaseConnector::disconnect(connection))
  d <- function(x) as.Date(x)
  put <- function(name, data) DatabaseConnector::insertTable(connection, databaseSchema = "main", tableName = name, data = data,
                                                             createTable = TRUE, dropTableIfExists = TRUE, progressBar = FALSE, camelCaseToSnakeCase = FALSE)
  ids <- 1:10
  put("person", data.frame(person_id = ids, year_of_birth = 1950, gender_concept_id = ifelse(ids %% 2 == 0, 8532, 8507)))
  start <- rep("2000-01-01", 10); end <- rep("2020-12-31", 10)
  end[4] <- "2010-03-01"; end[5] <- "2010-02-15"; start[7] <- "2009-06-01"
  put("observation_period", data.frame(person_id = ids, observation_period_start_date = d(start), observation_period_end_date = d(end)))
  put("concept", data.frame(concept_id = c(100, 101, 102), domain_id = "Condition"))
  put("concept_ancestor", data.frame(ancestor_concept_id = c(100, 100, 100, 101, 102), descendant_concept_id = c(100, 101, 102, 101, 102)))
  conditions <- data.frame(
    person_id = c(1, 1, 1, 2, 5, 6, 9),
    condition_concept_id = c(100, 100, 100, 101, 100, 100, 102), # 102 is excluded from the concept set
    condition_start_date = d(c("2010-01-06", "2010-04-16", "2005-01-01", "2010-07-20", "2010-01-31", "2005-01-01", "2010-06-01")))
  put("condition_occurrence", conditions)
  put("visit_occurrence", data.frame(person_id = 2, visit_concept_id = 9201, visit_start_date = d("2010-07-19"), visit_end_date = d("2010-07-22")))
  put("death", data.frame(person_id = 10, death_date = d("2010-06-01")))
  # Target (1): each person on 2010-01-01; person 8 outside observation, person 9 a second entry. Outcome (2).
  target <- data.frame(cohort_definition_id = 1, subject_id = c(1:10, 9), cohort_start_date = d(c(rep("2010-01-01", 7), "2021-06-01", "2010-01-01", "2010-01-01", "2012-01-01")))
  outcome <- data.frame(cohort_definition_id = 2, subject_id = c(1, 2, 5, 6, 9),
                        cohort_start_date = d(c("2010-01-06", "2010-07-20", "2010-01-31", "2005-01-01", "2012-01-02")))
  cohort <- rbind(target, outcome)
  cohort$cohort_end_date <- cohort$cohort_start_date
  put("cohort", cohort)
  details
}

outcomeJson <- '{"ConceptSets": [{"id": 0, "name": "o", "expression": {"items": [
  {"concept": {"CONCEPT_ID": 100}, "includeDescendants": true, "isExcluded": false},
  {"concept": {"CONCEPT_ID": 102}, "includeDescendants": false, "isExcluded": true}]}}],
  "PrimaryCriteria": {"CriteriaList": [{"ConditionOccurrence": {"CodesetId": 0}}]}}'
definitions <- data.frame(cohortId = c(1, 2), json = c("{}", outcomeJson), sql = c("target", "outcome"), stringsAsFactors = FALSE)
diagnosticTask <- data.frame(taskId = "t", targetId = 1, outcomeId = 2, riskWindowEnd = 365)

for (dbms in c("duckdb", "sqlite")) {
  test_that(paste("task diagnostics count every scenario on", dbms), {
    skip_if_not_installed(if (dbms == "duckdb") "duckdb" else "RSQLite")
    details <- syntheticCdm(dbms)
    out <- tempfile("sar_diag_")
    on.exit(unlink(out, recursive = TRUE), add = TRUE)
    runTaskDiagnostics(details, "main", "main", "cohort", diagnosticTask, definitions, "SYN", out, minCellCount = 1, tempEmulationSchema = NULL)
    folder <- file.path(out, "taskDiagnostics", "SYN")
    expect_false(file.exists(file.path(folder, "tasks", "t_error.txt")))
    s <- utils::read.csv(file.path(folder, "summary.csv"))
    expect_equal(s$target_entries, 11)    # person 9 twice
    expect_equal(s$target_persons, 10)
    expect_equal(s$in_observation, 9)     # person 8 outside observation
    expect_equal(s$with_washout, 8)       # person 7: 214 days of prior observation
    expect_equal(s$prior_outcome, 1)      # person 6
    expect_equal(s$population, 7)         # 1, 2, 3, 4, 5, 9, 10
    expect_equal(s$events, 3)             # 1, 2, 5; person 9's outcome is after the time at risk of the first entry
    expect_equal(s$non_cases_censored, 1) # person 4
    expect_equal(s$cases_censored, 1)     # person 5
    expect_equal(s$person_days, 5 + 200 + 365 + 59 + 30 + 365 + 365)
    expect_equal(s$events_with_hospital_care, 1) # person 2
    expect_equal(s$deaths_at_risk, 1)     # person 10
    expect_equal(s$records_found, 5)      # person 1: 3 records; persons 2 and 5: 1 each
    expect_equal(s$second_record_365, 1)  # person 1
    expect_equal(s$events_with_earlier_record, 1) # person 1
    expect_equal(s$recurrence_window_cases, 2)    # 1 and 2; person 5 left after 45 days
    timing <- utils::read.csv(file.path(folder, "timing.csv"))
    expect_equal(timing$events[match(c(7, 30, 365), timing$up_to_day)], c(1, 1, 1))
    codes <- utils::read.csv(file.path(folder, "codes.csv"))
    expect_equal(codes$persons[match(c(100, 101), codes$concept_id)], c(2, 1)) # at the event date; 102 is excluded
    strata <- utils::read.csv(file.path(folder, "strata.csv"))
    expect_equal(sum(strata$persons), 7)
    expect_equal(unique(strata$age_group), 60)
    calendar <- utils::read.csv(file.path(folder, "calendar.csv"))
    expect_equal(calendar$persons[calendar$index_year == 2010], 7)
    expect_equal(utils::read.csv(file.path(folder, "database.csv"))$deaths, 1)

    # restart: nothing is recomputed; small counts are censored when asked
    stamp <- file.mtime(file.path(folder, "tasks", "t.rds"))
    Sys.sleep(1.1)
    runTaskDiagnostics(details, "main", "main", "cohort", diagnosticTask, definitions, "SYN", out, minCellCount = 1, tempEmulationSchema = NULL)
    expect_identical(file.mtime(file.path(folder, "tasks", "t.rds")), stamp)
    runTaskDiagnostics(details, "main", "main", "cohort", diagnosticTask, definitions, "SYN", out, minCellCount = 5, tempEmulationSchema = NULL)
    expect_equal(utils::read.csv(file.path(folder, "summary.csv"))$events, -5)
  })
}

test_that("the scorecard turns diagnostics into gates and flags", {
  folder <- file.path(tempfile("sar_card_"), "DB1")
  dir.create(folder, recursive = TRUE)
  on.exit(unlink(dirname(folder), recursive = TRUE), add = TRUE)
  ok <- data.frame(taskId = c("good", "bad", "few"), target_entries = 2e5, target_persons = 1e5, in_observation = 1e5,
                   with_washout = c(9e4, 3e4, 9e4), prior_outcome = 0, prior_outcome_365 = 0, population = c(9e4, 3e4, 9e4),
                   events = c(4000, 4000, 100), non_cases_censored = c(1e4, 2e4, 1e4), cases_censored = 0,
                   person_days = 3e7, outcome_entries = 1e4, events_with_hospital_care = 100, deaths_at_risk = 10,
                   recurrence_window_cases = c(3000, 3000, 80), second_record_365 = c(2000, 100, 40),
                   events_with_earlier_record = c(10, 1000, 1), records_found = c(9000, 3000, 200))
  utils::write.csv(ok, file.path(folder, "summary.csv"), row.names = FALSE)
  timing <- rbind(data.frame(taskId = "good", up_to_day = c(7, 365), events = c(80, 3920)),
                  data.frame(taskId = "bad", up_to_day = c(7, 365), events = c(2000, 2000)),
                  data.frame(taskId = "few", up_to_day = c(7, 365), events = c(2, 98)))
  utils::write.csv(timing, file.path(folder, "timing.csv"), row.names = FALSE)
  years <- 2010:2017
  calendar <- rbind(data.frame(taskId = "good", index_year = years, persons = 1000),
                    data.frame(taskId = "bad", index_year = years, persons = c(1000, 1000, 1000, 5000, 5000, 5000, 5000, 5000)),
                    data.frame(taskId = "few", index_year = years, persons = 1000))
  utils::write.csv(calendar, file.path(folder, "calendar.csv"), row.names = FALSE)
  utils::write.csv(data.frame(persons = 1e6, deaths = 0, first_date = "2005-01-01", last_date = "2020-01-01"), file.path(folder, "database.csv"), row.names = FALSE)
  tasks <- data.frame(taskId = c("good", "bad", "few"), riskWindowEnd = 365, outcomeKind = c("chronic", "chronic", "death"))
  result <- taskScorecard(folder, tasks)
  card <- result$scorecard[match(c("good", "bad", "few"), result$scorecard$taskId), ]
  expect_equal(card$gateEvents, c(TRUE, TRUE, FALSE))
  expect_equal(card$gateDeath, c(TRUE, TRUE, FALSE)) # no death records in the database
  expect_identical(card$flags[1], "")
  for (f in c("attrition", "early-events", "earlier-records", "few-second-records", "calendar-jump", "follow-up")) {
    expect_match(card$flags[2], f, fixed = TRUE)
  }
  expect_equal(card$poolEvents[1], 3000)
  expect_false(any(result$tasks$usable)) # one database only: fewer than the 2 required
  two <- file.path(dirname(folder), "DB2")
  dir.create(two); file.copy(list.files(folder, full.names = TRUE), two)
  summary2 <- ok; summary2$events <- summary2$events * 7
  utils::write.csv(summary2, file.path(two, "summary.csv"), row.names = FALSE)
  both <- taskScorecard(c(folder, two), tasks)$tasks
  expect_equal(both$usable[both$taskId == "good"], TRUE)
  expect_equal(both$incidenceFlag[both$taskId == "good"], TRUE) # 7x between databases
})

test_that("the diagnostic SQL translates to every dialect", {
  statements <- character()
  run <- function(sql) statements <<- c(statements, sql)
  createDiagnosticPopulation(run, diagnosticTask)
  sql <- SqlRender::render(paste(statements, collapse = " "), schema = "s", cohort_table = "c", cdm = "cdm", target_id = 1, outcome_id = 2,
                           risk_start = 1, risk_end = 365, washout = 365, warnOnMissingParameters = FALSE)
  for (dbms in SqlRender::listSupportedDialects()$dialect) {
    translated <- suppressWarnings(SqlRender::translate(sql, dbms, tempEmulationSchema = "scratch"))
    expect_true(nzchar(translated), info = dbms)
    expect_false(grepl("@", translated, fixed = TRUE), info = dbms)
  }
  scalars <- scalarsSql(c(a = "SELECT COUNT(*) AS n FROM x", b = "SELECT COUNT(*) AS n FROM y"))
  expect_match(scalars, "CROSS JOIN", fixed = TRUE)
  expect_match(scalars, "^SELECT s1.n AS a, s2.n AS b FROM")
})
