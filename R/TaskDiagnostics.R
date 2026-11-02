# Task diagnostics (ideas from CohortDiagnostics): whether a candidate prediction task is usable in a database, from
# cohorts, outcomes and care records only; no covariates and no models, so tasks stay blind to model performance. Every
# statement goes through SqlRender, so it runs on all OHDSI dialects. taskScorecard() turns the results into gates and
# flags. What runTaskDiagnostics() writes is aggregate: counts below minCellCount are shown as -minCellCount.

DIAGNOSTIC_DOMAINS <- list(
  Condition = c("condition_occurrence", "condition_concept_id", "condition_start_date"),
  Drug = c("drug_exposure", "drug_concept_id", "drug_exposure_start_date"),
  Procedure = c("procedure_occurrence", "procedure_concept_id", "procedure_date"),
  Measurement = c("measurement", "measurement_concept_id", "measurement_date"),
  Observation = c("observation", "observation_concept_id", "observation_date"),
  Device = c("device_exposure", "device_concept_id", "device_exposure_start_date"),
  Visit = c("visit_occurrence", "visit_concept_id", "visit_start_date"))

# The concept-set items of a Circe cohort's entry events (concept id, include descendants, excluded).
entryConceptItems <- function(json) {
  expression <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  ids <- unlist(lapply(expression$PrimaryCriteria$CriteriaList, function(criterion) lapply(criterion, function(x) x$CodesetId)))
  rows <- lapply(expression$ConceptSets, function(set) {
    if (!set$id %in% ids) return(NULL)
    do.call(rbind, lapply(set$expression$items, function(item) {
      data.frame(concept_id = as.numeric(item$concept$CONCEPT_ID), include_descendants = as.integer(isTRUE(item$includeDescendants)),
                 is_excluded = as.integer(isTRUE(item$isExcluded)))
    }))
  })
  unique(do.call(rbind, rows))
}

# Counts below minCellCount shown as -minCellCount, as in exportResults().
censorCounts <- function(d, minCellCount) {
  for (column in names(d)) {
    if (is.numeric(d[[column]]) && !grepl("^(concept_id|index_year|age_group|gender_concept_id|target_id|other_target_id|up_to_day)$", column)) {
      small <- !is.na(d[[column]]) & d[[column]] > 0 & d[[column]] < minCellCount
      d[[column]][small] <- -minCellCount
    }
  }
  d
}

# Single values in one row, in SQL that needs no FROM-less SELECT (which not every dialect has).
scalarsSql <- function(values) {
  aliases <- paste0("s", seq_along(values))
  paste0("SELECT ", paste(sprintf("%s.n AS %s", aliases, names(values)), collapse = ", "), " FROM ",
         paste(sprintf("(%s) %s", values, aliases), collapse = " CROSS JOIN "))
}

# One task's population under the study's rule (first target entry, days of prior observation, no earlier outcome, at
# least a day at risk; D3 of the design: nobody is removed for leaving before the horizon) and its first event in the
# time at risk: temp tables #diag_pop and #diag_event.
createDiagnosticPopulation <- function(run, task) {
  run("
    DROP TABLE IF EXISTS #diag_first;
    SELECT subject_id, MIN(cohort_start_date) AS index_date
    INTO #diag_first
    FROM @schema.@cohort_table WHERE cohort_definition_id = @target_id
    GROUP BY subject_id;

    DROP TABLE IF EXISTS #diag_obs;
    SELECT f.subject_id, f.index_date,
           DATEDIFF(dd, op.observation_period_start_date, f.index_date) AS prior_days,
           DATEDIFF(dd, f.index_date, op.observation_period_end_date) AS follow_days
    INTO #diag_obs
    FROM #diag_first f
    INNER JOIN @cdm.observation_period op ON op.person_id = f.subject_id
      AND f.index_date >= op.observation_period_start_date AND f.index_date <= op.observation_period_end_date;

    DROP TABLE IF EXISTS #diag_prior;
    SELECT o.subject_id, MIN(DATEDIFF(dd, c.cohort_start_date, o.index_date)) AS days_since_latest
    INTO #diag_prior
    FROM #diag_obs o
    INNER JOIN @schema.@cohort_table c ON c.subject_id = o.subject_id AND c.cohort_definition_id = @outcome_id
      AND c.cohort_start_date < DATEADD(dd, @risk_start, o.index_date)
    GROUP BY o.subject_id;

    DROP TABLE IF EXISTS #diag_pop;
    SELECT o.subject_id, o.index_date, o.follow_days
    INTO #diag_pop
    FROM #diag_obs o
    LEFT JOIN #diag_prior p ON p.subject_id = o.subject_id
    WHERE o.prior_days >= @washout AND p.subject_id IS NULL AND o.follow_days >= @risk_start;

    DROP TABLE IF EXISTS #diag_event;
    SELECT p.subject_id, MIN(c.cohort_start_date) AS event_date
    INTO #diag_event
    FROM #diag_pop p
    INNER JOIN @schema.@cohort_table c ON c.subject_id = p.subject_id AND c.cohort_definition_id = @outcome_id
      AND c.cohort_start_date >= DATEADD(dd, @risk_start, p.index_date)
      AND c.cohort_start_date <= DATEADD(dd, @risk_end, p.index_date)
      AND DATEDIFF(dd, p.index_date, c.cohort_start_date) <= p.follow_days
    GROUP BY p.subject_id;

    DROP TABLE IF EXISTS #diag_risk;
    SELECT p.subject_id, p.index_date, p.follow_days, e.event_date,
           CASE WHEN e.event_date IS NULL THEN 0 ELSE 1 END AS has_event,
           CASE WHEN e.event_date IS NOT NULL THEN DATEDIFF(dd, p.index_date, e.event_date)
                WHEN p.follow_days < @risk_end THEN p.follow_days ELSE @risk_end END AS days_at_risk
    INTO #diag_risk
    FROM #diag_pop p
    LEFT JOIN #diag_event e ON e.subject_id = p.subject_id;")
}

# The cohort's entry-event records of a set of persons, into #diag_records (person, date, concept). Domains come from the
# vocabulary, so a concept set spanning several domains is covered. Returns FALSE when the set has no standard domain.
createEntryRecords <- function(connection, run, query, items, personsTable, tempEmulationSchema) {
  DatabaseConnector::insertTable(connection, tableName = "#diag_items", data = items, dropTableIfExists = TRUE, createTable = TRUE,
                                 tempTable = TRUE, tempEmulationSchema = tempEmulationSchema, progressBar = FALSE, camelCaseToSnakeCase = FALSE)
  run("
    DROP TABLE IF EXISTS #diag_concepts;
    SELECT DISTINCT c.concept_id, c.domain_id
    INTO #diag_concepts
    FROM (
      SELECT a.descendant_concept_id AS concept_id FROM #diag_items i
      INNER JOIN @cdm.concept_ancestor a ON a.ancestor_concept_id = i.concept_id
      WHERE i.include_descendants = 1 AND i.is_excluded = 0
      UNION ALL
      SELECT i.concept_id FROM #diag_items i WHERE i.include_descendants = 0 AND i.is_excluded = 0
    ) included
    INNER JOIN @cdm.concept c ON c.concept_id = included.concept_id
    WHERE c.concept_id NOT IN (
      SELECT a.descendant_concept_id FROM #diag_items i
      INNER JOIN @cdm.concept_ancestor a ON a.ancestor_concept_id = i.concept_id
      WHERE i.include_descendants = 1 AND i.is_excluded = 1
      UNION ALL
      SELECT i.concept_id FROM #diag_items i WHERE i.include_descendants = 0 AND i.is_excluded = 1
    );")
  domains <- intersect(query("SELECT DISTINCT domain_id FROM #diag_concepts")$domain_id, names(DIAGNOSTIC_DOMAINS))
  if (length(domains) == 0) return(FALSE)
  selects <- vapply(DIAGNOSTIC_DOMAINS[domains], function(d) {
    sprintf("SELECT r.person_id AS subject_id, r.%s AS record_date, r.%s AS concept_id
             FROM @cdm.%s r INNER JOIN %s p ON p.subject_id = r.person_id
             INNER JOIN #diag_concepts c ON c.concept_id = r.%s", d[3], d[2], d[1], personsTable, d[2])
  }, character(1))
  run(paste0("DROP TABLE IF EXISTS #diag_records;
              SELECT subject_id, record_date, concept_id INTO #diag_records FROM (",
             paste(selects, collapse = " UNION ALL "), ") records;"))
  TRUE
}

diagnosticQuery <- function(connection, sql) {
  d <- DatabaseConnector::querySql(connection, sql)
  names(d) <- tolower(names(d))
  d[] <- lapply(d, function(x) if (inherits(x, "integer64")) as.numeric(x) else x)
  d
}

taskDiagnostics <- function(connection, run, query, task, cohortDefinitionSet, tempEmulationSchema) {
  createDiagnosticPopulation(run, task)
  summary <- query(scalarsSql(c(
    target_entries = "SELECT COUNT(*) AS n FROM @schema.@cohort_table WHERE cohort_definition_id = @target_id",
    target_persons = "SELECT COUNT(*) AS n FROM #diag_first",
    in_observation = "SELECT COUNT(*) AS n FROM #diag_obs",
    with_washout = "SELECT COUNT(*) AS n FROM #diag_obs WHERE prior_days >= @washout",
    prior_outcome = "SELECT COUNT(*) AS n FROM #diag_prior",
    prior_outcome_365 = "SELECT COUNT(*) AS n FROM #diag_prior WHERE days_since_latest <= 365",
    population = "SELECT COUNT(*) AS n FROM #diag_risk",
    events = "SELECT COUNT(*) AS n FROM #diag_risk WHERE has_event = 1",
    non_cases_censored = "SELECT COUNT(*) AS n FROM #diag_risk WHERE has_event = 0 AND follow_days < @risk_end",
    cases_censored = "SELECT COUNT(*) AS n FROM #diag_risk WHERE has_event = 1 AND follow_days < @risk_end",
    person_days = "SELECT SUM(CAST(days_at_risk AS BIGINT)) AS n FROM #diag_risk",
    outcome_entries = "SELECT COUNT(*) AS n FROM @schema.@cohort_table WHERE cohort_definition_id = @outcome_id")))
  timing <- query("
    SELECT CASE WHEN d <= 7 THEN 7 WHEN d <= 30 THEN 30 WHEN d <= 90 THEN 90 WHEN d <= 180 THEN 180
                WHEN d <= 365 THEN 365 WHEN d <= 730 THEN 730 ELSE 9999 END AS up_to_day, COUNT(*) AS events
    FROM (SELECT DATEDIFF(dd, index_date, event_date) AS d FROM #diag_risk WHERE has_event = 1) e
    GROUP BY CASE WHEN d <= 7 THEN 7 WHEN d <= 30 THEN 30 WHEN d <= 90 THEN 90 WHEN d <= 180 THEN 180
                  WHEN d <= 365 THEN 365 WHEN d <= 730 THEN 730 ELSE 9999 END")
  strata <- query("
    SELECT FLOOR((YEAR(r.index_date) - p.year_of_birth) / 10) * 10 AS age_group, p.gender_concept_id,
           COUNT(*) AS persons, SUM(r.has_event) AS events, SUM(CAST(r.days_at_risk AS BIGINT)) AS person_days
    FROM #diag_risk r INNER JOIN @cdm.person p ON p.person_id = r.subject_id
    GROUP BY FLOOR((YEAR(r.index_date) - p.year_of_birth) / 10) * 10, p.gender_concept_id")
  calendar <- query("
    SELECT YEAR(index_date) AS index_year, COUNT(*) AS persons, SUM(has_event) AS events,
           SUM(CAST(days_at_risk AS BIGINT)) AS person_days
    FROM #diag_risk GROUP BY YEAR(index_date)")
  care <- query("
    SELECT COUNT(DISTINCT r.subject_id) AS events_with_hospital_care
    FROM #diag_risk r INNER JOIN @cdm.visit_occurrence v ON v.person_id = r.subject_id
      AND r.event_date >= v.visit_start_date AND r.event_date <= v.visit_end_date
      AND v.visit_concept_id IN (9201, 262, 9203, 8717, 8870, 8920)
    WHERE r.has_event = 1")
  death <- query("
    SELECT COUNT(DISTINCT r.subject_id) AS deaths_at_risk
    FROM #diag_risk r INNER JOIN @cdm.death d ON d.person_id = r.subject_id
      AND d.death_date > r.index_date AND d.death_date <= DATEADD(dd, @risk_end, r.index_date)")
  summary <- cbind(summary, care, death, recurrence_window_cases = NA_real_, second_record_365 = NA_real_,
                   events_with_earlier_record = NA_real_, records_found = NA_real_)
  codes <- NULL
  json <- cohortDefinitionSet$json[cohortDefinitionSet$cohortId == task$outcomeId]
  items <- if (length(json)) entryConceptItems(json) else NULL
  if (!is.null(items) && summary$events > 0) {
    run("DROP TABLE IF EXISTS #diag_cases; SELECT subject_id, index_date, event_date INTO #diag_cases FROM #diag_risk WHERE has_event = 1;")
    if (createEntryRecords(connection, run, query, items, "#diag_cases", tempEmulationSchema)) {
      evidence <- query("
        SELECT COUNT(*) AS records_found,
               COUNT(DISTINCT CASE WHEN r.record_date > c.event_date AND r.record_date <= DATEADD(dd, 365, c.event_date) THEN r.subject_id END) AS second_record_365,
               COUNT(DISTINCT CASE WHEN r.record_date < c.index_date THEN r.subject_id END) AS events_with_earlier_record
        FROM #diag_cases c INNER JOIN #diag_records r ON r.subject_id = c.subject_id")
      summary$records_found <- evidence$records_found
      summary$second_record_365 <- evidence$second_record_365
      summary$events_with_earlier_record <- evidence$events_with_earlier_record
      summary$recurrence_window_cases <- query("
        SELECT COUNT(*) AS n FROM #diag_cases c INNER JOIN #diag_risk r ON r.subject_id = c.subject_id
        WHERE r.follow_days >= DATEDIFF(dd, c.index_date, c.event_date) + 365")$n
      codes <- query("
        SELECT r.concept_id, COUNT(DISTINCT r.subject_id) AS persons
        FROM #diag_records r INNER JOIN #diag_cases c ON c.subject_id = r.subject_id AND r.record_date = c.event_date
        GROUP BY r.concept_id")
    }
  }
  list(summary = summary, timing = timing, strata = strata, calendar = calendar, codes = codes)
}

#' Task diagnostics
#'
#' Checks whether candidate tasks are usable in a database, from cohorts, outcomes and care records only: no covariates and
#' no models, so tasks stay blind to model performance. For each task, under the study's population rule (first target
#' entry, prior observation, no earlier outcome, nobody removed for leaving before the horizon): attrition, events, follow-up,
#' when events fall in the time at risk, age, sex and calendar patterns, hospital care at the outcome, deaths at risk, whether
#' cases have a second record of the outcome and earlier records of it, and the outcome's codes. Overlap between target
#' cohorts and what the database records are also written. Restartable: a task whose settings have not changed is skipped.
#'
#' @param connectionDetails DatabaseConnector connection details.
#' @param cdmDatabaseSchema Schema of the CDM tables and vocabulary.
#' @param cohortDatabaseSchema Schema of the cohort table.
#' @param cohortTable Cohort table, with the tasks' cohorts already generated.
#' @param tasks Data frame: `taskId`, `targetId`, `outcomeId`, `riskWindowEnd` (days), optionally `riskWindowStart` (default
#'   1), `washoutDays` (default 365).
#' @param cohortDefinitionSet Cohort definitions (`cohortId`, `json`) for the checks on the outcome's codes.
#' @param databaseName Name used in the results.
#' @param outputFolder Folder for the results: `taskDiagnostics/<database>/`.
#' @param minCellCount Counts below this are shown as `-minCellCount`.
#' @param tempEmulationSchema For platforms without temp tables (e.g. Spark, Oracle).
#' @return The folder with the results, invisibly.
#' @export
runTaskDiagnostics <- function(connectionDetails, cdmDatabaseSchema, cohortDatabaseSchema, cohortTable, tasks, cohortDefinitionSet,
                               databaseName, outputFolder, minCellCount = 5,
                               tempEmulationSchema = getOption("sqlRenderTempEmulationSchema")) {
  folder <- file.path(outputFolder, "taskDiagnostics", databaseName)
  dir.create(file.path(folder, "tasks"), recursive = TRUE, showWarnings = FALSE)
  connection <- DatabaseConnector::connect(connectionDetails)
  on.exit(DatabaseConnector::disconnect(connection), add = TRUE)
  defaults <- function(task, column, value) if (column %in% names(task) && !is.na(task[[column]])) task[[column]] else value
  dropTemps <- function(run) {
    run(paste(sprintf("DROP TABLE IF EXISTS #diag_%s;", c("first", "obs", "prior", "pop", "event", "risk", "cases", "items",
                                                          "concepts", "records")), collapse = " "))
  }
  for (i in seq_len(nrow(tasks))) {
    task <- tasks[i, ]
    file <- file.path(folder, "tasks", paste0(task$taskId, ".rds"))
    definition <- cohortDefinitionSet$sql[match(c(task$targetId, task$outcomeId), cohortDefinitionSet$cohortId)]
    settings <- list(database = databaseName, cdm = cdmDatabaseSchema, target = definition[1], outcome = definition[2],
                     riskStart = defaults(task, "riskWindowStart", 1), riskEnd = task$riskWindowEnd,
                     washout = defaults(task, "washoutDays", 365), minCellCount = minCellCount)
    hash <- rlang::hash(settings)
    if (file.exists(file) && identical(readRDS(file)$hash, hash)) next
    parameters <- list(schema = cohortDatabaseSchema, cohort_table = cohortTable, cdm = cdmDatabaseSchema, target_id = task$targetId,
                       outcome_id = task$outcomeId, risk_start = settings$riskStart, risk_end = settings$riskEnd, washout = settings$washout)
    render <- function(sql) do.call(SqlRender::render, c(list(sql), parameters, list(warnOnMissingParameters = FALSE)))
    translate <- function(sql) SqlRender::translate(render(sql), DatabaseConnector::dbms(connection), tempEmulationSchema = tempEmulationSchema)
    run <- function(sql) DatabaseConnector::executeSql(connection, translate(sql), progressBar = FALSE, reportOverallTime = FALSE)
    query <- function(sql) diagnosticQuery(connection, translate(sql))
    start <- Sys.time()
    result <- tryCatch(taskDiagnostics(connection, run, query, task, cohortDefinitionSet, tempEmulationSchema), error = function(e) e)
    if (inherits(result, "error")) {
      writeLines(conditionMessage(result), file.path(folder, "tasks", paste0(task$taskId, "_error.txt")))
      message(databaseName, " ", task$taskId, " diagnostics failed: ", conditionMessage(result))
      next
    }
    unlink(file.path(folder, "tasks", paste0(task$taskId, "_error.txt")))
    tryCatch(dropTemps(run), error = function(e) NULL)
    result <- lapply(result, function(d) if (is.null(d)) NULL else censorCounts(as.data.frame(d), minCellCount))
    saveRDS(c(result, list(hash = hash, taskId = task$taskId, seconds = as.numeric(difftime(Sys.time(), start, units = "secs")))), file)
    message(sprintf("%s %s: population %s, events %s, %.0f s", databaseName, task$taskId, result$summary$population, result$summary$events,
                    as.numeric(difftime(Sys.time(), start, units = "secs"))))
  }
  # Database level: the size and span of the database and whether deaths are recorded; target cohort overlap.
  parameters <- list(schema = cohortDatabaseSchema, cohort_table = cohortTable, cdm = cdmDatabaseSchema,
                     target_ids = paste(unique(tasks$targetId), collapse = ","))
  query <- function(sql) {
    sql <- do.call(SqlRender::render, c(list(sql), parameters, list(warnOnMissingParameters = FALSE)))
    diagnosticQuery(connection, SqlRender::translate(sql, DatabaseConnector::dbms(connection), tempEmulationSchema = tempEmulationSchema))
  }
  database <- query(scalarsSql(c(
    persons = "SELECT COUNT(*) AS n FROM @cdm.person",
    deaths = "SELECT COUNT(*) AS n FROM @cdm.death",
    first_date = "SELECT MIN(observation_period_start_date) AS n FROM @cdm.observation_period",
    last_date = "SELECT MAX(observation_period_end_date) AS n FROM @cdm.observation_period")))
  overlap <- query("
    SELECT a.cohort_definition_id AS target_id, b.cohort_definition_id AS other_target_id, COUNT(*) AS persons_in_both
    FROM (SELECT DISTINCT cohort_definition_id, subject_id FROM @schema.@cohort_table WHERE cohort_definition_id IN (@target_ids)) a
    INNER JOIN (SELECT DISTINCT cohort_definition_id, subject_id FROM @schema.@cohort_table WHERE cohort_definition_id IN (@target_ids)) b
      ON a.subject_id = b.subject_id AND a.cohort_definition_id <= b.cohort_definition_id
    GROUP BY a.cohort_definition_id, b.cohort_definition_id")
  utils::write.csv(database, file.path(folder, "database.csv"), row.names = FALSE)
  utils::write.csv(censorCounts(overlap, minCellCount), file.path(folder, "target_overlap.csv"), row.names = FALSE)
  # One CSV per kind of result, all tasks together.
  results <- lapply(list.files(file.path(folder, "tasks"), "\\.rds$", full.names = TRUE), readRDS)
  for (kind in c("summary", "timing", "strata", "calendar", "codes")) {
    rows <- lapply(results, function(r) if (!is.null(r[[kind]]) && nrow(r[[kind]])) cbind(taskId = r$taskId, r[[kind]]))
    utils::write.csv(do.call(rbind, rows), file.path(folder, paste0(kind, ".csv")), row.names = FALSE)
  }
  invisible(folder)
}

#' Thresholds of the task scorecard
#'
#' Starting values, to calibrate on tasks known to be usable before they decide anything.
#'
#' @return A named list.
#' @export
taskScorecardThresholds <- function() {
  list(minPoolEvents = 1500, minTestEvents = 500, testFraction = 0.25, maxAttritionStep = 0.5, earlyDays = 7, earlyRatio = 5,
       minHorizonForEarly = 31, maxEarlierRecord = 0.05, minSecondRecord = 0.3, maxCalendarRatio = 2, minCalendarEvents = 20,
       maxNonCasesCensored = 0.5, maxIncidenceRatio = 5, minDatabases = 2)
}

#' Task scorecard
#'
#' Gates and flags per task and database from [runTaskDiagnostics()], and across databases. A failed gate rules a task out in
#' that database; a flag asks a reviewer to look. They detect broken or unsuitable tasks; they do not establish that a
#' phenotype is valid.
#'
#' @param folders Folders written by [runTaskDiagnostics()], one per database.
#' @param tasks Data frame: `taskId`, `riskWindowEnd`; optionally `outcomeKind` (`"chronic"`: flag few second records;
#'   `"death"`: requires usable death records).
#' @param thresholds From [taskScorecardThresholds()].
#' @return A list: `scorecard` (one row per task and database) and `tasks` (one row per task).
#' @export
taskScorecard <- function(folders, tasks, thresholds = taskScorecardThresholds()) {
  th <- thresholds
  read <- function(folder, kind) {
    file <- file.path(folder, paste0(kind, ".csv"))
    if (file.exists(file) && file.size(file) > 5) utils::read.csv(file, stringsAsFactors = FALSE) else NULL
  }
  # Censored counts (negative) are small: read as NA in ratios.
  positive <- function(x) ifelse(!is.na(x) & x < 0, NA, x)
  rows <- lapply(folders, function(folder) {
    database <- basename(folder)
    summary <- read(folder, "summary")
    if (is.null(summary)) return(NULL)
    timing <- read(folder, "timing"); calendar <- read(folder, "calendar"); info <- read(folder, "database")
    d <- merge(tasks[, intersect(c("taskId", "riskWindowEnd", "outcomeKind"), names(tasks))], summary, by = "taskId")
    d$database <- database
    events <- positive(d$events); population <- positive(d$population)
    d$poolEvents <- round(events * (1 - th$testFraction)); d$testEvents <- round(events * th$testFraction)
    d$eventRate <- events / population
    d$incidencePer1000Years <- 1000 * events / (positive(d$person_days) / 365.25)
    # step = people removed by one attrition step
    steps <- cbind(lostObservation = 1 - d$in_observation / d$target_persons, lostWashout = 1 - d$with_washout / d$in_observation,
                   lostPriorOutcome = 1 - population / d$with_washout)
    d$largestAttrition <- suppressWarnings(apply(steps, 1, function(x) if (all(is.na(x))) NA else max(x, na.rm = TRUE)))
    d$nonCasesCensoredShare <- positive(d$non_cases_censored) / (population - events)
    d$earlyShare <- vapply(d$taskId, function(id) {
      if (is.null(timing)) return(NA_real_)
      x <- timing[timing$taskId == id, ]
      n <- sum(positive(x$events), na.rm = TRUE)
      if (n == 0) NA_real_ else sum(positive(x$events[x$up_to_day <= th$earlyDays]), na.rm = TRUE) / n
    }, numeric(1))
    d$earlyRatio <- d$earlyShare / pmin(1, th$earlyDays / d$riskWindowEnd)
    d$earlierRecordShare <- positive(d$events_with_earlier_record) / events
    d$secondRecordShare <- positive(d$second_record_365) / positive(d$recurrence_window_cases)
    d$hospitalCareShare <- positive(d$events_with_hospital_care) / events
    d$calendarJump <- vapply(d$taskId, function(id) {
      if (is.null(calendar)) return(NA_real_)
      x <- calendar[calendar$taskId == id, ]
      x <- x[order(x$index_year), ]
      x <- x[!x$index_year %in% 2020:2021 & !is.na(x$persons) & x$persons >= th$minCalendarEvents, ]
      # first and last years are often partial
      if (nrow(x) < 5) return(NA_real_)
      x <- x[-c(1, nrow(x)), ]
      r <- x$persons[-1] / x$persons[-nrow(x)]
      max(r, 1 / r)
    }, numeric(1))
    d$deathsRecorded <- !is.null(info) && !is.na(info$deaths) && info$deaths > 0
    # Gates
    d$gateEvents <- !is.na(d$poolEvents) & !is.na(d$testEvents) & d$poolEvents >= th$minPoolEvents & d$testEvents >= th$minTestEvents
    kind <- if ("outcomeKind" %in% names(d)) d$outcomeKind else rep("", nrow(d))
    d$gateDeath <- kind != "death" | d$deathsRecorded
    # Flags
    flag <- function(condition, name) ifelse(!is.na(condition) & condition, name, "")
    d$flags <- trimws(paste(
      flag(d$largestAttrition > th$maxAttritionStep, "attrition"),
      flag(d$earlyRatio > th$earlyRatio & d$riskWindowEnd >= th$minHorizonForEarly, "early-events"),
      flag(d$earlierRecordShare > th$maxEarlierRecord, "earlier-records"),
      flag(kind == "chronic" & d$secondRecordShare < th$minSecondRecord, "few-second-records"),
      flag(d$calendarJump > th$maxCalendarRatio, "calendar-jump"),
      flag(d$nonCasesCensoredShare > th$maxNonCasesCensored, "follow-up"),
      flag(is.na(d$records_found) | d$records_found == 0, "no-code-records")))
    d$flags <- gsub(" +", " ", d$flags)
    d
  })
  scorecard <- do.call(rbind, rows)
  perTask <- do.call(rbind, lapply(split(scorecard, scorecard$taskId), function(d) {
    usable <- d$gateEvents & d$gateDeath
    rates <- d$incidencePer1000Years[usable & !is.na(d$incidencePer1000Years) & d$incidencePer1000Years > 0]
    data.frame(taskId = d$taskId[1], databases = nrow(d), usableDatabases = sum(usable), usable = sum(usable) >= th$minDatabases,
               incidenceRatio = if (length(rates) > 1) max(rates) / min(rates) else NA_real_,
               incidenceFlag = length(rates) > 1 && max(rates) / min(rates) > th$maxIncidenceRatio,
               flags = paste(sort(unique(unlist(strsplit(d$flags, " ")))), collapse = " "), stringsAsFactors = FALSE)
  }))
  rownames(perTask) <- NULL
  list(scorecard = scorecard, tasks = perTask)
}
