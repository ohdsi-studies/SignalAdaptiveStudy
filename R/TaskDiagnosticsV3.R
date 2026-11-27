# ==============================================================================
# TaskDiagnosticsV3.R
# Production-ready OHDSI SqlRender extension to SignalAdaptiveStudy::TaskDiagnostics
# (Modules 1-7, 100% aggregate-only, minCellCount = 5, multi-dialect OHDSI SQL).
#
# Companion to SignalAdaptiveStudy/R/TaskDiagnostics.R:
#   - Every SQL statement uses SQL Server OHDSI dialect translated via
#     SqlRender::render() + SqlRender::translate(..., dbms, tempEmulationSchema)
#     -> runs cleanly on PostgreSQL, Redshift, Spark/Databricks, Snowflake,
#        BigQuery, SQL Server, Oracle, and DuckDB.
#   - Strictly aggregate-only output (zero patient-level rows; counts 1..(minCellCount-1)
#     censored to -minCellCount; derived percentages blanked when numerators are
#     small; site x year cells stay in memory to compute 5-fold CV
#     auroc_site_plus_year & pct_zero_event_sites without exporting care_site cells).
#   - Restartable per-task .rds cache (hashed on target/outcome SQL + JSON + settings)
#     and zips only aggregate CSVs for safe return from production sites.
#
# Modules:
#   M1: Target AND Outcome top standard index concepts (concept_id, name, vocab)
#   M2: Database-native orphan concept scan (unincluded standard concepts >= 500p)
#   M3: Clinical temporal profile around outcome date ([-365,-31], [-30,-1], d0, [+1,+30])
#       for Conditions & ATC 4th Drug classes
#   M4: Detailed visit context on outcome date (% Inpatient/ED/OP & % on visit_start_date)
#   M5: Site + Calendar Year confounding 5-fold CV AUROC (auroc_site_plus_year,
#       auroc_year_only) and pct_zero_event_sites
#   M6: 3-way censoring (complete non-cases, administrative end-of-DB censoring,
#       true LTFU, competing death) + post-observation-period deaths in TAR
#   M7: Drug route-of-administration audit (Clinical Drug Form via concept_ancestor)
#       + Measurement unit distribution audit (unit_concept_id mix; no MIN/MAX)
# ==============================================================================

DIAGNOSTICS_V3_VERSION <- 3L

V3_DOMAINS <- list(
  Condition   = c("condition_occurrence", "condition_concept_id", "condition_start_date"),
  Drug        = c("drug_exposure",        "drug_concept_id",      "drug_exposure_start_date"),
  Procedure   = c("procedure_occurrence", "procedure_concept_id", "procedure_date"),
  Measurement = c("measurement",          "measurement_concept_id", "measurement_date"),
  Observation = c("observation",          "observation_concept_id", "observation_date"),
  Device      = c("device_exposure",      "device_concept_id",    "device_exposure_start_date"),
  Visit       = c("visit_occurrence",     "visit_concept_id",     "visit_start_date")
)

V3_STOPWORDS <- c(
  "of", "in", "to", "for", "the", "and", "or", "with", "without", "due", "by", "on", "at", "from",
  "disorder", "disease", "finding", "procedure", "measurement", "drug", "exposure", "unspecified",
  "other", "nos", "nec", "acute", "chronic", "primary", "secondary", "history", "initial", "subsequent",
  "stage", "type", "syndrome", "diagnosis", "episode", "patient", "patients", "adults", "adult", "persons",
  "recorded", "newly", "first", "new", "starting", "prior", "no", "non", "after", "visit", "hospital",
  "care", "requiring", "related", "associated", "induced", "left", "right", "bilateral", "general"
)

v3IsSmallCount <- function(x, minCellCount) {
  !is.na(x) & x > 0 & x < minCellCount
}

v3CensorCounts <- function(d, minCellCount = 5L) {
  if (is.null(d) || !is.data.frame(d) || nrow(d) == 0L) return(d)
  skip_pat <- "^(concept_id|unit_concept_id|cohort_id|target_id|outcome_id|index_year|risk_window_start|risk_window_end|horizon_days|n_orphan_candidates)$"
  pct_pat  <- "(pct|rate|ratio|auroc|median|mean|min_val|max_val|p05|p50|p95)"
  censored_rows <- rep(FALSE, nrow(d))
  for (col in names(d)) {
    if (is.numeric(d[[col]]) && !grepl(skip_pat, col) && !grepl(pct_pat, col)) {
      small <- v3IsSmallCount(d[[col]], minCellCount)
      if (any(small)) {
        d[[col]][small] <- -as.numeric(minCellCount)
        if (col %in% c("n_persons", "n_cases", "any_n", "n_rows")) censored_rows <- censored_rows | small
      }
    }
  }
  # Prevent back-calculation of censored counts from derived percentages / means on detail tables
  if (any(censored_rows)) {
    for (col in names(d)) {
      if (is.numeric(d[[col]]) && grepl("^(pct_|mean_)", col)) {
        d[[col]][censored_rows] <- NA_real_
      }
    }
  }
  # On summary table: if events or population is censored, blank derived rates/ratios
  if ("events" %in% names(d)) {
    evSmall <- !is.na(d$events) & d$events < 0
    if (any(evSmall)) {
      for (col in intersect(c("event_rate_pct", "early_d1_7_pct", "early_ratio"), names(d))) {
        d[[col]][evSmall] <- NA_real_
      }
    }
  }
  if ("population" %in% names(d)) {
    popSmall <- !is.na(d$population) & d$population < 0
    if (any(popSmall)) {
      for (col in intersect(c("event_rate_pct", "mean_age", "pct_female", "pct_complete_noncase",
                              "pct_admin_censored", "pct_true_ltfu", "pct_competing_death"), names(d))) {
        d[[col]][popSmall] <- NA_real_
      }
    }
  }
  d
}

v3EntryConceptItems <- function(jsonStr) {
  if (is.null(jsonStr) || length(jsonStr) == 0L || is.na(jsonStr) || !nzchar(jsonStr)) return(NULL)
  expr <- tryCatch(jsonlite::fromJSON(jsonStr, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(expr) || is.null(expr$PrimaryCriteria$CriteriaList) || is.null(expr$ConceptSets)) return(NULL)
  ids <- unlist(lapply(expr$PrimaryCriteria$CriteriaList, function(crit) lapply(crit, function(x) x$CodesetId)))
  if (length(ids) == 0L) return(NULL)
  rows <- lapply(expr$ConceptSets, function(set) {
    if (is.null(set$id) || !set$id %in% ids || is.null(set$expression$items)) return(NULL)
    do.call(rbind, lapply(set$expression$items, function(item) {
      if (is.null(item$concept$CONCEPT_ID)) return(NULL)
      data.frame(
        concept_id = as.numeric(item$concept$CONCEPT_ID),
        concept_name = ifelse(is.null(item$concept$CONCEPT_NAME), "", as.character(item$concept$CONCEPT_NAME)),
        include_descendants = as.integer(isTRUE(item$includeDescendants)),
        is_excluded = as.integer(isTRUE(item$isExcluded)),
        stringsAsFactors = FALSE
      )
    }))
  })
  res <- do.call(rbind, rows)
  if (is.null(res) || nrow(res) == 0L) return(NULL)
  unique(res)
}

v3AutoKeywords <- function(itemsDf, fallbackText = "") {
  toks <- character(0)
  if (!is.null(itemsDf) && nrow(itemsDf) > 0L) {
    inc <- itemsDf[itemsDf$is_excluded == 0L, , drop = FALSE]
    for (nm in inc$concept_name) {
      w <- unlist(regmatches(tolower(nm), gregexpr("[a-z]{5,}", tolower(nm))))
      w <- setdiff(w, V3_STOPWORDS)
      for (x in utils::head(w, 3L)) if (!x %in% toks) toks <- c(toks, x)
      if (length(toks) >= 4L) break
    }
  }
  if (length(toks) == 0L && nzchar(fallbackText)) {
    w <- unlist(regmatches(tolower(fallbackText), gregexpr("[a-z]{5,}", tolower(fallbackText))))
    toks <- utils::head(unique(setdiff(w, V3_STOPWORDS)), 3L)
  }
  utils::head(toks, 4L)
}

# Portable 5-fold CV AUROC from aggregate (care_site_id, idx_yr, pos, neg) counts (runs in R; zero small-cell export)
v3ComputeSiteYearAuroc <- function(cellDf, useSite = TRUE, useYear = TRUE, seed = 20261127L) {
  if (is.null(cellDf) || nrow(cellDf) == 0L) return(NA_real_)
  cellDf <- cellDf[!is.na(cellDf$care_site_id) & !is.na(cellDf$idx_yr) & !is.na(cellDf$pos) & !is.na(cellDf$neg), , drop = FALSE]
  if (nrow(cellDf) == 0L) return(NA_real_)
  if (sum(cellDf$pos) < 20 || sum(cellDf$neg) < 20) return(NA_real_)
  set.seed(seed)
  K <- 5L
  splitCounts <- function(cntVec) {
    mat <- matrix(0L, nrow = length(cntVec), ncol = K)
    for (i in seq_along(cntVec)) {
      n <- as.integer(cntVec[i])
      if (is.na(n) || n <= 0L) next
      base <- n %/% K; rem <- n %% K
      mat[i, ] <- base
      if (rem > 0L) mat[i, sample.int(K, rem)] <- base + 1L
    }
    mat
  }
  posF <- splitCounts(cellDf$pos)
  negF <- splitCounts(cellDf$neg)
  wAuc <- function(sc, pw, nw) {
    keep <- (pw + nw) > 0 & !is.na(sc)
    sc <- sc[keep]; pw <- as.numeric(pw[keep]); nw <- as.numeric(nw[keep])
    if (sum(pw) == 0 || sum(nw) == 0) return(NA_real_)
    ord <- order(sc); sc <- sc[ord]; pw <- pw[ord]; nw <- nw[ord]
    g <- cumsum(c(TRUE, diff(sc) != 0))
    pg <- as.numeric(tapply(pw, g, sum))
    ng <- as.numeric(tapply(nw, g, sum))
    sum(pg * (c(0, cumsum(ng)[-length(ng)]) + 0.5 * ng)) / (sum(pw) * sum(nw))
  }
  logit <- function(p) log(pmin(pmax(p, 1e-6), 1 - 1e-6) / (1 - pmin(pmax(p, 1e-6), 1 - 1e-6)))
  siteKey <- as.character(cellDf$care_site_id)
  yrKey   <- as.character(cellDf$idx_yr)
  aucs <- sapply(seq_len(K), function(k) {
    pTr <- rowSums(posF[, -k, drop = FALSE])
    nTr <- rowSums(negF[, -k, drop = FALSE])
    p0  <- (sum(pTr) + 1) / (sum(pTr + nTr) + 2)
    sc  <- rep(0, nrow(cellDf))
    if (useSite) {
      sr <- (tapply(pTr, siteKey, sum) + 20 * p0) / (tapply(pTr + nTr, siteKey, sum) + 20)
      sc <- sc + logit(as.numeric(sr[siteKey]))
    }
    if (useYear) {
      yr <- (tapply(pTr, yrKey, sum) + 20 * p0) / (tapply(pTr + nTr, yrKey, sum) + 20)
      sc <- sc + logit(as.numeric(yr[yrKey]))
    }
    if (useSite && useYear) sc <- sc - logit(p0)
    wAuc(sc, posF[, k], negF[, k])
  })
  if (all(is.na(aucs))) return(NA_real_)
  round(mean(aucs, na.rm = TRUE), 4)
}

v3Query <- function(connection, sql) {
  d <- DatabaseConnector::querySql(connection, sql)
  names(d) <- tolower(names(d))
  d[] <- lapply(d, function(x) if (inherits(x, "integer64")) as.numeric(x) else x)
  d
}

v3BuildConceptSetTemp <- function(connection, run, query, itemsDf, tempEmulationSchema) {
  run("DROP TABLE IF EXISTS #v3_items; DROP TABLE IF EXISTS #v3_cs_ids;")
  if (is.null(itemsDf) || nrow(itemsDf) == 0L) return(character(0))
  DatabaseConnector::insertTable(
    connection, tableName = "#v3_items",
    data = itemsDf[, c("concept_id", "include_descendants", "is_excluded")],
    dropTableIfExists = TRUE, createTable = TRUE, tempTable = TRUE,
    tempEmulationSchema = tempEmulationSchema, progressBar = FALSE, camelCaseToSnakeCase = FALSE
  )
  run("
    SELECT DISTINCT c.concept_id, c.domain_id
    INTO #v3_cs_ids
    FROM (
      SELECT a.descendant_concept_id AS concept_id FROM #v3_items i
      INNER JOIN @vocab.concept_ancestor a ON a.ancestor_concept_id = i.concept_id
      WHERE i.include_descendants = 1 AND i.is_excluded = 0
      UNION ALL
      SELECT i.concept_id FROM #v3_items i WHERE i.include_descendants = 0 AND i.is_excluded = 0
    ) inc
    INNER JOIN @vocab.concept c ON c.concept_id = inc.concept_id
    WHERE c.concept_id NOT IN (
      SELECT a.descendant_concept_id FROM #v3_items i
      INNER JOIN @vocab.concept_ancestor a ON a.ancestor_concept_id = i.concept_id
      WHERE i.include_descendants = 1 AND i.is_excluded = 1
      UNION ALL
      SELECT i.concept_id FROM #v3_items i WHERE i.include_descendants = 0 AND i.is_excluded = 1
    );
  ")
  query("SELECT DISTINCT domain_id FROM #v3_cs_ids")$domain_id
}

# Build one-time per-database standard concept prevalence table (#v3_concept_counts, >=500 persons) for M2 orphan scans
v3PrepareConceptCounts <- function(run, minOrphanPersons = 500L) {
  run(sprintf("
    DROP TABLE IF EXISTS #v3_concept_counts;
    SELECT domain_id, concept_id, person_count
    INTO #v3_concept_counts
    FROM (
      SELECT 'Condition' AS domain_id, condition_concept_id AS concept_id, COUNT(DISTINCT person_id) AS person_count
      FROM @cdm.condition_occurrence GROUP BY condition_concept_id HAVING COUNT(DISTINCT person_id) >= %d
      UNION ALL
      SELECT 'Drug' AS domain_id, drug_concept_id AS concept_id, COUNT(DISTINCT person_id) AS person_count
      FROM @cdm.drug_exposure GROUP BY drug_concept_id HAVING COUNT(DISTINCT person_id) >= %d
      UNION ALL
      SELECT 'Procedure' AS domain_id, procedure_concept_id AS concept_id, COUNT(DISTINCT person_id) AS person_count
      FROM @cdm.procedure_occurrence GROUP BY procedure_concept_id HAVING COUNT(DISTINCT person_id) >= %d
      UNION ALL
      SELECT 'Measurement' AS domain_id, measurement_concept_id AS concept_id, COUNT(DISTINCT person_id) AS person_count
      FROM @cdm.measurement GROUP BY measurement_concept_id HAVING COUNT(DISTINCT person_id) >= %d
    ) cc;
  ", minOrphanPersons, minOrphanPersons, minOrphanPersons, minOrphanPersons))
}

v3IndexCodes <- function(connection, run, query, taskId, role, cohortId, sampleTbl, itemsDf, tempEmulationSchema) {
  doms <- intersect(v3BuildConceptSetTemp(connection, run, query, itemsDf, tempEmulationSchema), names(V3_DOMAINS))
  if (length(doms) == 0L) return(NULL)
  nTot <- query(sprintf("SELECT COUNT(DISTINCT subject_id) AS n FROM %s", sampleTbl))$n[1]
  if (is.na(nTot) || nTot == 0) return(NULL)
  parts <- vapply(doms, function(d) {
    m <- V3_DOMAINS[[d]]
    sprintf("
      SELECT r.person_id AS subject_id, r.%s AS concept_id
      FROM @cdm.%s r
      INNER JOIN %s s ON s.subject_id = r.person_id AND r.%s = s.event_dt
      INNER JOIN #v3_cs_ids cs ON cs.concept_id = r.%s
    ", m[2], m[1], sampleTbl, m[3], m[2])
  }, character(1))
  df <- query(sprintf("
    SELECT concept_id, concept_name, vocabulary_id, domain_id, n_persons
    FROM (
      SELECT h.concept_id, c.concept_name, c.vocabulary_id, c.domain_id,
             COUNT(DISTINCT h.subject_id) AS n_persons,
             ROW_NUMBER() OVER (ORDER BY COUNT(DISTINCT h.subject_id) DESC, h.concept_id) AS rn
      FROM (%s) h
      INNER JOIN @vocab.concept c ON c.concept_id = h.concept_id
      GROUP BY h.concept_id, c.concept_name, c.vocabulary_id, c.domain_id
    ) x
    WHERE rn <= 10
    ORDER BY n_persons DESC;
  ", paste(parts, collapse = " UNION ALL ")))
  if (nrow(df) == 0L) return(NULL)
  df$pct_of_sample <- round(100.0 * df$n_persons / nTot, 2)
  cbind(task_id = taskId, role = role, cohort_id = cohortId, df, stringsAsFactors = FALSE)
}

v3OrphanScan <- function(connection, run, query, taskId, role, cohortId, itemsDf, fallbackText, tempEmulationSchema) {
  kws <- v3AutoKeywords(itemsDf, fallbackText)
  if (length(kws) == 0L) return(NULL)
  doms <- intersect(v3BuildConceptSetTemp(connection, run, query, itemsDf, tempEmulationSchema), c("Condition", "Drug", "Procedure", "Measurement"))
  if (length(doms) == 0L) return(NULL)
  likeClause <- paste(sprintf("LOWER(c.concept_name) LIKE '%%%s%%'", gsub("'", "''", kws)), collapse = " OR ")
  domIn <- paste(sprintf("'%s'", doms), collapse = ",")
  df <- query(sprintf("
    SELECT concept_id, concept_name, domain_id, vocabulary_id, person_count
    FROM (
      SELECT c.concept_id, c.concept_name, c.domain_id, c.vocabulary_id, cc.person_count,
             ROW_NUMBER() OVER (ORDER BY cc.person_count DESC, c.concept_id) AS rn
      FROM #v3_concept_counts cc
      INNER JOIN @vocab.concept c ON c.concept_id = cc.concept_id
      WHERE cc.domain_id IN (%s)
        AND c.standard_concept = 'S'
        AND (%s)
        AND c.concept_id NOT IN (SELECT concept_id FROM #v3_cs_ids)
    ) x
    WHERE rn <= 8
    ORDER BY person_count DESC;
  ", domIn, likeClause))
  if (nrow(df) == 0L) return(NULL)
  cbind(task_id = taskId, role = role, cohort_id = cohortId, df, stringsAsFactors = FALSE)
}

v3TemporalProfile <- function(query, taskId, sampleTbl, minCellCount = 5L) {
  nSamp <- query(sprintf("SELECT COUNT(*) AS n FROM %s;", sampleTbl))$n[1]
  if (is.na(nSamp) || nSamp < 25) return(NULL)
  qCond <- sprintf("
    SELECT domain, concept_id, concept_name, n_m365_m31, n_m30_m1, n_day0, n_p1_p30, any_n
    FROM (
      SELECT 'Condition' AS domain, e.concept_id, c.concept_name,
             COUNT(DISTINCT CASE WHEN e.win = 'm365_m31' THEN e.subject_id END) AS n_m365_m31,
             COUNT(DISTINCT CASE WHEN e.win = 'm30_m1'   THEN e.subject_id END) AS n_m30_m1,
             COUNT(DISTINCT CASE WHEN e.win = 'day0'     THEN e.subject_id END) AS n_day0,
             COUNT(DISTINCT CASE WHEN e.win = 'p1_p30'   THEN e.subject_id END) AS n_p1_p30,
             COUNT(DISTINCT e.subject_id) AS any_n,
             ROW_NUMBER() OVER (ORDER BY COUNT(DISTINCT CASE WHEN e.win = 'day0' THEN e.subject_id END) DESC, COUNT(DISTINCT e.subject_id) DESC, e.concept_id) AS rn
      FROM (
        SELECT s.subject_id, co.condition_concept_id AS concept_id,
               CASE WHEN co.condition_start_date BETWEEN DATEADD(dd, -365, s.event_dt) AND DATEADD(dd, -31, s.event_dt) THEN 'm365_m31'
                    WHEN co.condition_start_date BETWEEN DATEADD(dd, -30,  s.event_dt) AND DATEADD(dd, -1,  s.event_dt) THEN 'm30_m1'
                    WHEN co.condition_start_date = s.event_dt THEN 'day0'
                    WHEN co.condition_start_date BETWEEN DATEADD(dd, 1,    s.event_dt) AND DATEADD(dd, 30,  s.event_dt) THEN 'p1_p30'
               END AS win
        FROM %s s
        INNER JOIN @cdm.condition_occurrence co ON co.person_id = s.subject_id
         AND co.condition_start_date BETWEEN DATEADD(dd, -365, s.event_dt) AND DATEADD(dd, 30, s.event_dt)
      ) e
      INNER JOIN @vocab.concept c ON c.concept_id = e.concept_id
      WHERE e.win IS NOT NULL
      GROUP BY e.concept_id, c.concept_name
    ) x
    WHERE rn <= 8
    ORDER BY n_day0 DESC, any_n DESC;
  ", sampleTbl)
  qDrug <- sprintf("
    SELECT domain, concept_id, concept_name, n_m365_m31, n_m30_m1, n_day0, n_p1_p30, any_n
    FROM (
      SELECT 'Drug_ATC4' AS domain, e.concept_id, c.concept_name,
             COUNT(DISTINCT CASE WHEN e.win = 'm365_m31' THEN e.subject_id END) AS n_m365_m31,
             COUNT(DISTINCT CASE WHEN e.win = 'm30_m1'   THEN e.subject_id END) AS n_m30_m1,
             COUNT(DISTINCT CASE WHEN e.win = 'day0'     THEN e.subject_id END) AS n_day0,
             COUNT(DISTINCT CASE WHEN e.win = 'p1_p30'   THEN e.subject_id END) AS n_p1_p30,
             COUNT(DISTINCT e.subject_id) AS any_n,
             ROW_NUMBER() OVER (ORDER BY COUNT(DISTINCT CASE WHEN e.win = 'day0' THEN e.subject_id END) DESC, COUNT(DISTINCT e.subject_id) DESC, e.concept_id) AS rn
      FROM (
        SELECT s.subject_id, ca.ancestor_concept_id AS concept_id,
               CASE WHEN d.drug_exposure_start_date BETWEEN DATEADD(dd, -365, s.event_dt) AND DATEADD(dd, -31, s.event_dt) THEN 'm365_m31'
                    WHEN d.drug_exposure_start_date BETWEEN DATEADD(dd, -30,  s.event_dt) AND DATEADD(dd, -1,  s.event_dt) THEN 'm30_m1'
                    WHEN d.drug_exposure_start_date = s.event_dt THEN 'day0'
                    WHEN d.drug_exposure_start_date BETWEEN DATEADD(dd, 1,    s.event_dt) AND DATEADD(dd, 30,  s.event_dt) THEN 'p1_p30'
               END AS win
        FROM %s s
        INNER JOIN @cdm.drug_exposure d ON d.person_id = s.subject_id
         AND d.drug_exposure_start_date BETWEEN DATEADD(dd, -365, s.event_dt) AND DATEADD(dd, 30, s.event_dt)
        INNER JOIN @vocab.concept_ancestor ca ON d.drug_concept_id = ca.descendant_concept_id
        INNER JOIN @vocab.concept ac ON ca.ancestor_concept_id = ac.concept_id
         AND ac.vocabulary_id = 'ATC' AND ac.concept_class_id = 'ATC 4th'
      ) e
      INNER JOIN @vocab.concept c ON c.concept_id = e.concept_id
      WHERE e.win IS NOT NULL
      GROUP BY e.concept_id, c.concept_name
    ) x
    WHERE rn <= 6
    ORDER BY n_day0 DESC, any_n DESC;
  ", sampleTbl)
  df <- rbind(query(qCond), query(qDrug))
  if (nrow(df) == 0L) return(NULL)
  safePct <- function(cnt) {
    pct <- round(100.0 * cnt / nSamp, 1)
    pct[v3IsSmallCount(cnt, minCellCount)] <- NA_real_
    pct
  }
  data.frame(
    task_id      = taskId,
    domain       = df$domain,
    concept_id   = df$concept_id,
    concept_name = df$concept_name,
    pct_m365_m31 = safePct(df$n_m365_m31),
    pct_m30_m1   = safePct(df$n_m30_m1),
    pct_day0     = safePct(df$n_day0),
    pct_p1_p30   = safePct(df$n_p1_p30),
    any_n        = df$any_n,
    stringsAsFactors = FALSE
  )
}

# Module 7: Drug route-of-administration audit (Clinical Drug Form breakdown) + Lab unit distribution audit
v3DrugRouteAudit <- function(connection, run, query, taskId, role, cohortId, sampleTbl, itemsDf, tempEmulationSchema) {
  doms <- v3BuildConceptSetTemp(connection, run, query, itemsDf, tempEmulationSchema)
  if (!"Drug" %in% doms) return(NULL)
  nTot <- query(sprintf("SELECT COUNT(DISTINCT subject_id) AS n FROM %s;", sampleTbl))$n[1]
  if (is.na(nTot) || nTot < 5) return(NULL)
  df <- query(sprintf("
    SELECT concept_id, clinical_drug_form, n_persons
    FROM (
      SELECT fc.concept_id, fc.concept_name AS clinical_drug_form,
             COUNT(DISTINCT s.subject_id) AS n_persons,
             ROW_NUMBER() OVER (ORDER BY COUNT(DISTINCT s.subject_id) DESC, fc.concept_id) AS rn
      FROM %s s
      INNER JOIN @cdm.drug_exposure d ON d.person_id = s.subject_id AND d.drug_exposure_start_date = s.event_dt
      INNER JOIN #v3_cs_ids cs ON cs.concept_id = d.drug_concept_id
      INNER JOIN @vocab.concept_ancestor ca ON ca.descendant_concept_id = d.drug_concept_id
      INNER JOIN @vocab.concept fc ON fc.concept_id = ca.ancestor_concept_id
       AND fc.concept_class_id = 'Clinical Drug Form' AND fc.standard_concept = 'S'
      GROUP BY fc.concept_id, fc.concept_name
    ) x
    WHERE rn <= 8
    ORDER BY n_persons DESC;
  ", sampleTbl))
  if (nrow(df) == 0L) return(NULL)
  df$pct_of_sample <- round(100.0 * df$n_persons / nTot, 2)
  cbind(task_id = taskId, role = role, cohort_id = cohortId, df, stringsAsFactors = FALSE)
}

v3LabUnitAudit <- function(connection, run, query, taskId, role, cohortId, sampleTbl, itemsDf, tempEmulationSchema) {
  doms <- v3BuildConceptSetTemp(connection, run, query, itemsDf, tempEmulationSchema)
  if (!"Measurement" %in% doms) return(NULL)
  # Note: intentionally omits MIN/MAX to prevent extreme single-patient outlier disclosure
  df <- query(sprintf("
    SELECT concept_id, concept_name, unit_concept_id, unit_name, n_rows, n_persons, mean_val
    FROM (
      SELECT m.measurement_concept_id AS concept_id, mc.concept_name,
             COALESCE(m.unit_concept_id, 0) AS unit_concept_id,
             COALESCE(uc.concept_name, 'NO_UNIT') AS unit_name,
             COUNT(*) AS n_rows, COUNT(DISTINCT s.subject_id) AS n_persons,
             AVG(CAST(m.value_as_number AS FLOAT)) AS mean_val,
             ROW_NUMBER() OVER (ORDER BY COUNT(DISTINCT s.subject_id) DESC, m.measurement_concept_id) AS rn
      FROM %s s
      INNER JOIN @cdm.measurement m ON m.person_id = s.subject_id AND m.measurement_date = s.event_dt
      INNER JOIN #v3_cs_ids cs ON cs.concept_id = m.measurement_concept_id
      INNER JOIN @vocab.concept mc ON mc.concept_id = m.measurement_concept_id
      LEFT JOIN @vocab.concept uc ON uc.concept_id = m.unit_concept_id
      GROUP BY m.measurement_concept_id, mc.concept_name, COALESCE(m.unit_concept_id, 0), COALESCE(uc.concept_name, 'NO_UNIT')
    ) x
    WHERE rn <= 8
    ORDER BY n_persons DESC;
  ", sampleTbl))
  if (nrow(df) == 0L) return(NULL)
  df$mean_val <- round(as.numeric(df$mean_val), 2)
  cbind(task_id = taskId, role = role, cohort_id = cohortId, df, stringsAsFactors = FALSE)
}

#' Run TaskDiagnosticsV3 (Modules 1-7) on a CDM database
#'
#' Drop-in OHDSI SqlRender companion to [runTaskDiagnostics()].
#' Exports strictly aggregate-only CSVs (`minCellCount` enforced; site x year cells
#' stay in memory to compute `auroc_site_plus_year` and `pct_zero_event_sites` without
#' small-cell export).
#'
#' @param connectionDetails DatabaseConnector connectionDetails (or supply an open connection via `connection`).
#' @param cdmDatabaseSchema Schema containing OMOP CDM tables.
#' @param vocabDatabaseSchema Schema containing OMOP vocabulary tables (defaults to `cdmDatabaseSchema`).
#' @param cohortDatabaseSchema Schema containing the already-generated cohort table.
#' @param cohortTable Name of the cohort table.
#' @param tasks Data frame with `taskId`, `targetId`, `outcomeId`, `riskWindowEnd`,
#'   and optional `riskWindowStart` (default 1), `washoutDays` (default 365),
#'   `removePriorOutcome` (default TRUE), `outcomeKind` (default ""), `targetDesc`, `outcomeDesc`.
#' @param cohortDefinitionSet Data frame with `cohortId`, `json` (and optional `sql`).
#' @param databaseName Database identifier for output folder & CSVs.
#' @param outputFolder Output directory; writes `taskDiagnosticsV3/<databaseName>/` + zip archive.
#' @param minCellCount Small-cell suppression threshold (default 5; counts `1..(minCellCount-1)` -> `-minCellCount`).
#' @param tempEmulationSchema Temp schema for Spark/Oracle/BigQuery.
#' @param runOrphanScan Logical; set TRUE (default) to build `#v3_concept_counts` once and run M2 orphan scans.
#' @param minOrphanPersons Minimum database person count for `#v3_concept_counts` (default 500).
#' @param connection Optional already-open DatabaseConnector connection (skips connect/disconnect if supplied).
#' @return Path to the aggregate-only zip archive (invisibly).
#' @export
runTaskDiagnosticsV3 <- function(connectionDetails = NULL, cdmDatabaseSchema,
                                 vocabDatabaseSchema = cdmDatabaseSchema,
                                 cohortDatabaseSchema, cohortTable, tasks,
                                 cohortDefinitionSet, databaseName, outputFolder,
                                 minCellCount = 5L,
                                 tempEmulationSchema = getOption("sqlRenderTempEmulationSchema"),
                                 runOrphanScan = TRUE,
                                 minOrphanPersons = 500L,
                                 connection = NULL) {
  folder <- file.path(outputFolder, "taskDiagnosticsV3", databaseName)
  dir.create(file.path(folder, "tasks"), recursive = TRUE, showWarnings = FALSE)
  ownConn <- is.null(connection)
  if (ownConn) {
    connection <- DatabaseConnector::connect(connectionDetails)
    on.exit(DatabaseConnector::disconnect(connection), add = TRUE)
  }
  dbms <- DatabaseConnector::dbms(connection)

  baseParams <- list(cdm = cdmDatabaseSchema, vocab = vocabDatabaseSchema,
                     schema = cohortDatabaseSchema, cohort_table = cohortTable)
  execSql <- function(sql, extra = list()) {
    r <- do.call(SqlRender::render, c(list(sql), baseParams, extra, list(warnOnMissingParameters = FALSE)))
    t <- SqlRender::translate(r, targetDialect = dbms, tempEmulationSchema = tempEmulationSchema)
    DatabaseConnector::executeSql(connection, t, progressBar = FALSE, reportOverallTime = FALSE)
  }
  querySql <- function(sql, extra = list()) {
    r <- do.call(SqlRender::render, c(list(sql), baseParams, extra, list(warnOnMissingParameters = FALSE)))
    t <- SqlRender::translate(r, targetDialect = dbms, tempEmulationSchema = tempEmulationSchema)
    v3Query(connection, t)
  }

  # Accept both SignalAdaptiveStudy camelCase (`taskId`, `targetId`, `outcomeId`, `riskWindowEnd`) and snake_case
  colPick <- function(row, cands, def) {
    for (c in cands) if (c %in% names(row) && !is.na(row[[c]])) return(row[[c]])
    def
  }

  taskMeta <- lapply(seq_len(nrow(tasks)), function(i) {
    tk <- tasks[i, ]
    tid <- as.character(colPick(tk, c("taskId", "task_id"), paste0("task_", i)))
    tId <- as.numeric(colPick(tk, c("targetId", "target_id", "target_cid"), NA_real_))
    oId <- as.numeric(colPick(tk, c("outcomeId", "outcome_id", "outcome_cid"), NA_real_))
    tIdx <- match(tId, cohortDefinitionSet$cohortId)
    oIdx <- match(oId, cohortDefinitionSet$cohortId)
    tJson <- if (!is.na(tIdx) && "json" %in% names(cohortDefinitionSet)) cohortDefinitionSet$json[tIdx] else NA_character_
    oJson <- if (!is.na(oIdx) && "json" %in% names(cohortDefinitionSet)) cohortDefinitionSet$json[oIdx] else NA_character_
    tSql  <- if (!is.na(tIdx) && "sql"  %in% names(cohortDefinitionSet)) cohortDefinitionSet$sql[tIdx]  else NA_character_
    oSql  <- if (!is.na(oIdx) && "sql"  %in% names(cohortDefinitionSet)) cohortDefinitionSet$sql[oIdx]  else NA_character_
    rStart <- as.integer(colPick(tk, c("riskWindowStart", "risk_window_start", "tar_start"), 1L))
    rEnd   <- as.integer(colPick(tk, c("riskWindowEnd", "risk_window_end", "tar_end"), 365L))
    wash   <- as.integer(colPick(tk, c("washoutDays", "washout_days"), 365L))
    rmPri  <- !isFALSE(as.logical(colPick(tk, c("removePriorOutcome", "remove_prior_outcome"), TRUE)))
    okind  <- as.character(colPick(tk, c("outcomeKind", "outcome_kind"), ""))
    tDesc  <- as.character(colPick(tk, c("targetDesc", "target_desc"), ""))
    oDesc  <- as.character(colPick(tk, c("outcomeDesc", "outcome_desc"), ""))
    rdsFile <- file.path(folder, "tasks", paste0(tid, ".rds"))
    cfgHash <- rlang::hash(list(
      version = DIAGNOSTICS_V3_VERSION, database = databaseName, cdm = cdmDatabaseSchema,
      tId = tId, oId = oId, tJson = tJson, oJson = oJson, tSql = tSql, oSql = oSql,
      rStart = rStart, rEnd = rEnd, wash = wash, rmPri = rmPri,
      minCellCount = as.integer(minCellCount), runOrphanScan = isTRUE(runOrphanScan),
      minOrphanPersons = as.integer(minOrphanPersons)
    ))
    cached <- file.exists(rdsFile) && identical(tryCatch(readRDS(rdsFile)$hash, error = function(e) NULL), cfgHash)
    list(tid = tid, tId = tId, oId = oId, tJson = tJson, oJson = oJson,
         rStart = rStart, rEnd = rEnd, wash = wash, rmPri = rmPri,
         okind = okind, tDesc = tDesc, oDesc = oDesc,
         rdsFile = rdsFile, cfgHash = cfgHash, cached = cached)
  })

  needRun <- vapply(taskMeta, function(m) !m$cached, logical(1))
  if (any(needRun)) {
    dbLastDate <- as.Date(querySql("SELECT MAX(observation_period_end_date) AS d FROM @cdm.observation_period")$d[1])
    dbLastY <- as.integer(format(dbLastDate, "%Y"))
    dbLastM <- as.integer(format(dbLastDate, "%m"))
    dbLastD <- as.integer(format(dbLastDate, "%d"))
    if (isTRUE(runOrphanScan)) {
      message(sprintf("[%s] Building one-time database concept_counts table for M2 orphan scans...", databaseName))
      v3PrepareConceptCounts(execSql, minOrphanPersons = as.integer(minOrphanPersons))
      on.exit(try(execSql("DROP TABLE IF EXISTS #v3_concept_counts;"), silent = TRUE), add = TRUE)
    }
  }

  for (m in taskMeta[needRun]) {
    tid <- m$tid; tId <- m$tId; oId <- m$oId
    tJson <- m$tJson; oJson <- m$oJson
    rStart <- m$rStart; rEnd <- m$rEnd; wash <- m$wash; rmPri <- m$rmPri
    okind <- m$okind; tDesc <- m$tDesc; oDesc <- m$oDesc
    rdsFile <- m$rdsFile; cfgHash <- m$cfgHash

    t0 <- Sys.time()
    pExtra <- list(target_id = tId, outcome_id = oId,
                   risk_start = rStart, risk_end = rEnd, washout = wash,
                   db_last_y = dbLastY, db_last_m = dbLastM, db_last_d = dbLastD)
    runT   <- function(sql) execSql(sql, pExtra)
    queryT <- function(sql) querySql(sql, pExtra)

    priorClause <- if (rmPri) {
      "AND NOT EXISTS (SELECT 1 FROM @schema.@cohort_table po WHERE po.cohort_definition_id = @outcome_id AND po.subject_id = f.subject_id AND po.cohort_start_date < DATEADD(dd, @risk_start, f.index_date))"
    } else ""

    res <- tryCatch({
      runT(sprintf("
        DROP TABLE IF EXISTS #v3_first; DROP TABLE IF EXISTS #v3_elig; DROP TABLE IF EXISTS #v3_ev;
        DROP TABLE IF EXISTS #v3_samp_t; DROP TABLE IF EXISTS #v3_samp_o;

        SELECT subject_id, MIN(cohort_start_date) AS index_date
        INTO #v3_first
        FROM @schema.@cohort_table WHERE cohort_definition_id = @target_id
        GROUP BY subject_id;

        SELECT f.subject_id, f.index_date, op.observation_period_end_date AS op_end,
               COALESCE(p.care_site_id, 0) AS care_site_id, p.gender_concept_id,
               YEAR(f.index_date) - p.year_of_birth AS age_at_idx,
               YEAR(f.index_date) AS idx_yr,
               DATEDIFF(dd, f.index_date, op.observation_period_end_date) AS follow_days
        INTO #v3_elig
        FROM #v3_first f
        INNER JOIN @cdm.observation_period op ON op.person_id = f.subject_id
          AND f.index_date >= DATEADD(dd, @washout, op.observation_period_start_date)
          AND f.index_date <= op.observation_period_end_date
        INNER JOIN @cdm.person p ON p.person_id = f.subject_id
        WHERE DATEDIFF(dd, f.index_date, op.observation_period_end_date) >= @risk_start
          %s;

        SELECT e.subject_id, e.index_date, e.op_end, e.care_site_id, e.gender_concept_id,
               e.age_at_idx, e.idx_yr, e.follow_days, o.ev_dt,
               CASE WHEN o.ev_dt IS NOT NULL THEN 1 ELSE 0 END AS has_event,
               CASE WHEN o.ev_dt IS NOT NULL THEN DATEDIFF(dd, DATEADD(dd, @risk_start, e.index_date), o.ev_dt) + 1 ELSE NULL END AS tar_day,
               CASE WHEN o.ev_dt IS NOT NULL THEN DATEDIFF(dd, DATEADD(dd, @risk_start, e.index_date), o.ev_dt) + 1
                    WHEN e.follow_days < @risk_end THEN e.follow_days - @risk_start + 1
                    ELSE @risk_end - @risk_start + 1 END AS tar_days_obs,
               CASE WHEN d.death_date IS NOT NULL AND d.death_date >= DATEADD(dd, @risk_start, e.index_date) AND d.death_date <= DATEADD(dd, @risk_end, e.index_date) THEN 1 ELSE 0 END AS died_in_tar,
               CASE WHEN d.death_date IS NOT NULL AND d.death_date > e.op_end AND d.death_date <= DATEADD(dd, @risk_end, e.index_date) THEN 1 ELSE 0 END AS died_post_obs
        INTO #v3_ev
        FROM #v3_elig e
        LEFT JOIN (
          SELECT e2.subject_id, MIN(c.cohort_start_date) AS ev_dt
          FROM #v3_elig e2
          INNER JOIN @schema.@cohort_table c ON c.cohort_definition_id = @outcome_id AND c.subject_id = e2.subject_id
            AND c.cohort_start_date >= DATEADD(dd, @risk_start, e2.index_date)
            AND c.cohort_start_date <= DATEADD(dd, @risk_end, e2.index_date)
            AND c.cohort_start_date <= e2.op_end
          GROUP BY e2.subject_id
        ) o ON o.subject_id = e.subject_id
        LEFT JOIN (
          SELECT person_id, MIN(death_date) AS death_date
          FROM @cdm.death
          GROUP BY person_id
        ) d ON d.person_id = e.subject_id;

        SELECT subject_id, index_date AS event_dt
        INTO #v3_samp_t
        FROM (SELECT subject_id, index_date, ROW_NUMBER() OVER (ORDER BY subject_id) AS rn FROM #v3_ev) x
        WHERE rn <= 25000;

        SELECT subject_id, ev_dt AS event_dt, index_date
        INTO #v3_samp_o
        FROM (SELECT subject_id, ev_dt, index_date, ROW_NUMBER() OVER (ORDER BY subject_id) AS rn FROM #v3_ev WHERE has_event = 1) x
        WHERE rn <= 5000;
      ", priorClause))

      sm <- queryT("
        SELECT COUNT(*) AS population,
               SUM(has_event) AS events,
               AVG(CAST(age_at_idx AS FLOAT)) AS mean_age,
               SUM(CASE WHEN gender_concept_id = 8532 THEN 1 ELSE 0 END) AS n_female,
               SUM(CASE WHEN has_event = 1 AND tar_day <= 7 THEN 1 ELSE 0 END) AS events_d1_7,
               SUM(CASE WHEN has_event = 0 AND follow_days >= @risk_end THEN 1 ELSE 0 END) AS complete_noncases,
               SUM(CASE WHEN has_event = 0 AND died_in_tar = 1 THEN 1 ELSE 0 END) AS competing_deaths,
               SUM(CASE WHEN has_event = 0 AND died_in_tar = 0 AND follow_days < @risk_end AND op_end >= DATEADD(dd, -90, DATEFROMPARTS(@db_last_y, @db_last_m, @db_last_d)) THEN 1 ELSE 0 END) AS admin_censored,
               SUM(CASE WHEN has_event = 0 AND died_in_tar = 0 AND follow_days < @risk_end AND op_end <  DATEADD(dd, -90, DATEFROMPARTS(@db_last_y, @db_last_m, @db_last_d)) THEN 1 ELSE 0 END) AS true_ltfu,
               SUM(died_post_obs) AS post_obs_deaths_in_tar,
               SUM(CAST(CASE WHEN tar_days_obs > 0 THEN tar_days_obs ELSE 1 END AS BIGINT)) AS person_days
        FROM #v3_ev;
      ")
      pop <- as.numeric(sm$population[1])
      ev  <- as.numeric(ifelse(is.na(sm$events[1]), 0, sm$events[1]))
      evD17 <- as.numeric(ifelse(is.na(sm$events_d1_7[1]), 0, sm$events_d1_7[1]))
      nFem  <- as.numeric(ifelse(is.na(sm$n_female[1]), 0, sm$n_female[1]))
      nComp <- as.numeric(ifelse(is.na(sm$complete_noncases[1]), 0, sm$complete_noncases[1]))
      nAdm  <- as.numeric(ifelse(is.na(sm$admin_censored[1]), 0, sm$admin_censored[1]))
      nLtfu <- as.numeric(ifelse(is.na(sm$true_ltfu[1]), 0, sm$true_ltfu[1]))
      nDead <- as.numeric(ifelse(is.na(sm$competing_deaths[1]), 0, sm$competing_deaths[1]))

      # Module 5: Site + Year confounding (aggregate cell table stays strictly in memory!)
      cells <- queryT("SELECT care_site_id, idx_yr, SUM(has_event) AS pos, COUNT(*) - SUM(has_event) AS neg FROM #v3_ev GROUP BY care_site_id, idx_yr;")
      syAuc <- v3ComputeSiteYearAuroc(cells, useSite = TRUE, useYear = TRUE)
      yrAuc <- v3ComputeSiteYearAuroc(cells, useSite = FALSE, useYear = TRUE)
      siteAgg <- queryT("SELECT care_site_id, COUNT(*) AS n, SUM(has_event) AS ev FROM #v3_ev GROUP BY care_site_id;")
      nSites <- nrow(siteAgg)
      zeroSitePct <- if (nSites > 1 && pop > 0) round(100.0 * sum(siteAgg$n[siteAgg$ev == 0]) / pop, 2) else 0.0

      # Module 4: Visit context on outcome date (pure ANSI ROW_NUMBER, no scalar subquery or TOP N)
      nSampO <- queryT("SELECT COUNT(*) AS n FROM #v3_samp_o;")$n[1]
      vc <- if (ev >= minCellCount && !is.na(nSampO) && nSampO > 0) {
        vcRaw <- queryT("
          SELECT visit_type, n_cases, n_on_vstart, n_rows
          FROM (
            SELECT v.visit_type, COUNT(DISTINCT v.subject_id) AS n_cases,
                   SUM(v.eq_vstart) AS n_on_vstart, COUNT(*) AS n_rows,
                   ROW_NUMBER() OVER (ORDER BY COUNT(DISTINCT v.subject_id) DESC, v.visit_type) AS rn
            FROM (
              SELECT s.subject_id, COALESCE(c.concept_name, 'No_Visit_On_Day') AS visit_type,
                     MAX(CASE WHEN vo.visit_start_date = s.event_dt THEN 1 ELSE 0 END) AS eq_vstart
              FROM #v3_samp_o s
              LEFT JOIN @cdm.visit_occurrence vo ON vo.person_id = s.subject_id
               AND s.event_dt BETWEEN vo.visit_start_date AND vo.visit_end_date
              LEFT JOIN @vocab.concept c ON c.concept_id = vo.visit_concept_id
              GROUP BY s.subject_id, COALESCE(c.concept_name, 'No_Visit_On_Day')
            ) v
            GROUP BY v.visit_type
          ) x
          WHERE rn <= 6
          ORDER BY n_cases DESC;
        ")
        if (nrow(vcRaw) > 0L) {
          pctVstart <- round(100.0 * vcRaw$n_on_vstart / pmax(1, vcRaw$n_rows), 1)
          pctVstart[v3IsSmallCount(vcRaw$n_on_vstart, minCellCount)] <- NA_real_
          data.frame(
            task_id = tid,
            visit_type = vcRaw$visit_type,
            n_cases = vcRaw$n_cases,
            pct_cases = round(100.0 * vcRaw$n_cases / nSampO, 1),
            pct_on_vstart = pctVstart,
            stringsAsFactors = FALSE
          )
        } else NULL
      } else NULL

      # Modules 1, 2, 3, 7
      tItems <- v3EntryConceptItems(tJson)
      oItems <- v3EntryConceptItems(oJson)
      cT <- if (pop >= minCellCount) v3IndexCodes(connection, runT, queryT, tid, "target",  tId,  "#v3_samp_t", tItems, tempEmulationSchema) else NULL
      cO <- if (ev  >= minCellCount) v3IndexCodes(connection, runT, queryT, tid, "outcome", oId,  "#v3_samp_o", oItems, tempEmulationSchema) else NULL
      codesDf <- rbind(cT, cO)

      oT <- if (isTRUE(runOrphanScan)) v3OrphanScan(connection, runT, queryT, tid, "target",  tId,  tItems, tDesc, tempEmulationSchema) else NULL
      oO <- if (isTRUE(runOrphanScan)) v3OrphanScan(connection, runT, queryT, tid, "outcome", oId,  oItems, oDesc, tempEmulationSchema) else NULL
      orphDf <- rbind(oT, oO)
      nOrph <- if (!is.null(orphDf)) nrow(orphDf) else 0L

      tempDf <- if (ev >= 25) v3TemporalProfile(queryT, tid, "#v3_samp_o", minCellCount = minCellCount) else NULL

      drT <- if (pop >= minCellCount) v3DrugRouteAudit(connection, runT, queryT, tid, "target",  tId,  "#v3_samp_t", tItems, tempEmulationSchema) else NULL
      drO <- if (ev  >= minCellCount) v3DrugRouteAudit(connection, runT, queryT, tid, "outcome", oId,  "#v3_samp_o", oItems, tempEmulationSchema) else NULL
      drugRouteDf <- rbind(drT, drO)

      luT <- if (pop >= minCellCount) v3LabUnitAudit(connection, runT, queryT, tid, "target",  tId,  "#v3_samp_t", tItems, tempEmulationSchema) else NULL
      luO <- if (ev  >= minCellCount) v3LabUnitAudit(connection, runT, queryT, tid, "outcome", oId,  "#v3_samp_o", oItems, tempEmulationSchema) else NULL
      labUnitDf <- rbind(luT, luO)

      horizonDays <- rEnd - rStart + 1L
      earlyShare <- if (ev > 0) evD17 / ev else NA_real_
      constShare <- if (pop > 0) min(1, 7 * pop / max(1, as.numeric(sm$person_days[1]))) else NA_real_
      earlyRatio <- if (!is.na(earlyShare) && ev >= 10 && horizonDays >= 31L && !is.na(constShare) && constShare > 0) round(earlyShare / constShare, 2) else NA_real_

      flags <- c()
      if (!is.na(syAuc) && syAuc >= 0.70) flags <- c(flags, "site-year-auroc>=0.70")
      if (zeroSitePct > 25 && ev >= 50)   flags <- c(flags, "zero-sites>25%")
      if (!is.na(earlyRatio) && earlyRatio > 5) flags <- c(flags, "early-events-spike")
      if (nOrph > 0) flags <- c(flags, sprintf("orphans(%d)", nOrph))
      verdict <- if (ev < 100) "LOW_EVENTS"
                 else if (any(grepl("site-year-auroc|zero-sites", flags))) "CONFOUNDING_FLAG"
                 else if (any(grepl("early-events", flags))) "REVIEW_FLAG"
                 else "PASS_V3"

      # Suppress small-count numerators in summary percentages and top-code strings
      earlyPctOut   <- if (!is.na(earlyShare) && !v3IsSmallCount(evD17, minCellCount)) round(100 * earlyShare, 1) else NA_real_
      earlyRatioOut <- if (!is.na(earlyRatio) && !v3IsSmallCount(evD17, minCellCount)) earlyRatio else NA_real_
      safePopPct <- function(num) {
        if (pop <= 0 || v3IsSmallCount(num, minCellCount) || v3IsSmallCount(pop - num, minCellCount)) NA_real_
        else round(100.0 * num / pop, 1)
      }
      topT <- if (!is.null(cT) && nrow(cT) > 0L && !v3IsSmallCount(cT$n_persons[1], minCellCount)) {
        sprintf("%s (%.0f%%)", substr(cT$concept_name[1], 1, 35), cT$pct_of_sample[1])
      } else ""
      topO <- if (!is.null(cO) && nrow(cO) > 0L && !v3IsSmallCount(cO$n_persons[1], minCellCount)) {
        sprintf("%s (%.0f%%)", substr(cO$concept_name[1], 1, 35), cO$pct_of_sample[1])
      } else ""

      sumRow <- data.frame(
        task_id = tid, database = databaseName, outcome_kind = okind,
        risk_window_start = rStart, risk_window_end = rEnd, horizon_days = horizonDays,
        population = pop, events = ev,
        event_rate_pct = if (pop > 0) round(100 * ev / pop, 3) else 0,
        mean_age = round(as.numeric(sm$mean_age[1]), 1),
        pct_female = safePopPct(nFem),
        early_d1_7_pct = earlyPctOut,
        early_ratio = earlyRatioOut,
        pct_complete_noncase = safePopPct(nComp),
        pct_admin_censored   = safePopPct(nAdm),
        pct_true_ltfu        = safePopPct(nLtfu),
        pct_competing_death  = safePopPct(nDead),
        post_obs_deaths_in_tar = sm$post_obs_deaths_in_tar[1],
        auroc_site_plus_year = syAuc, auroc_year_only = yrAuc,
        pct_zero_event_sites = zeroSitePct,
        top_target_code = topT, top_outcome_code = topO,
        n_orphan_candidates = nOrph,
        v3_flags = paste(flags, collapse = "; "), v3_verdict = verdict,
        stringsAsFactors = FALSE
      )
      list(summary = sumRow, index_codes = codesDf, orphan_concepts = orphDf,
           temporal_profile = tempDf, visit_context = vc,
           drug_route = drugRouteDf, lab_units = labUnitDf)
    }, error = function(e) e)

    try(runT("DROP TABLE IF EXISTS #v3_first; DROP TABLE IF EXISTS #v3_elig; DROP TABLE IF EXISTS #v3_ev; DROP TABLE IF EXISTS #v3_samp_t; DROP TABLE IF EXISTS #v3_samp_o; DROP TABLE IF EXISTS #v3_items; DROP TABLE IF EXISTS #v3_cs_ids;"), silent = TRUE)
    if (inherits(res, "error")) {
      writeLines(conditionMessage(res), file.path(folder, "tasks", paste0(tid, "_error.txt")))
      message(sprintf("[%s] %s FAILED: %s", databaseName, tid, conditionMessage(res)))
      next
    }
    unlink(file.path(folder, "tasks", paste0(tid, "_error.txt")))
    res <- lapply(res, function(df) v3CensorCounts(df, minCellCount))
    saveRDS(c(res, list(hash = cfgHash, taskId = tid, seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1))), rdsFile)
    message(sprintf("[%s] %-14s | %-16s | pop=%8.0f ev=%6.0f | syAUC=%5.3f 0site=%4.1f%% | %.1fs",
                    databaseName, tid, res$summary$v3_verdict[1], res$summary$population[1], res$summary$events[1],
                    ifelse(is.na(res$summary$auroc_site_plus_year[1]), -1, res$summary$auroc_site_plus_year[1]),
                    res$summary$pct_zero_event_sites[1], as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }

  # Combine per-task .rds for the requested tasks into aggregate CSVs & zip
  wantFiles <- vapply(taskMeta, function(m) m$rdsFile, character(1))
  rdsFiles <- wantFiles[file.exists(wantFiles)]
  allRes <- lapply(rdsFiles, readRDS)
  csvKinds <- c("summary", "index_codes", "orphan_concepts", "temporal_profile", "visit_context", "drug_route", "lab_units")
  written <- character(0)
  for (kind in csvKinds) {
    parts <- lapply(allRes, function(r) r[[kind]])
    parts <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x) > 0L, parts)
    if (length(parts) > 0L) {
      outFile <- file.path(folder, paste0(kind, "_v3.csv"))
      utils::write.csv(do.call(rbind, parts), outFile, row.names = FALSE)
      written <- c(written, outFile)
    }
  }
  zipPath <- file.path(outputFolder, sprintf("TaskDiagnosticsV3_%s_%s.zip", databaseName, format(Sys.Date(), "%Y%m%d")))
  unlink(zipPath)
  if (length(written) > 0L) {
    if (requireNamespace("zip", quietly = TRUE)) {
      zip::zipr(zipPath, written)
    } else {
      utils::zip(zipPath, written, flags = "-j9X")
    }
  }
  message(sprintf("[%s] Wrote %d aggregate V3 CSVs -> %s", databaseName, length(written), zipPath))
  invisible(zipPath)
}
