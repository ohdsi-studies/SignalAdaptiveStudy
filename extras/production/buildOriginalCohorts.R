# Circe cohort definitions for the original tasks, from the original concept sets (inst/original/definitions_v1.json),
# for the production data source (R/CdmDataSource.R). Writes inst/cohorts/<cohortId>.json and .sql (OHDSI SQL from
# CirceR, so production needs neither CirceR nor Java) and inst/settings/cohorts.csv.
#
# Targets (as the original cohort SQL):
#   - a disease: its first qualifying record, where qualifying means age >= 18 (Circe: year of index - year of birth)
#     and >= 365 days of prior observation; the cohort ends at the end of the observation period.
#   - adults: the first observation period of >= 365 days, entered 365 days after its start, aged >= 18 at entry.
# Outcomes: every record of the concept set (one-day episodes). Washout, time at risk and complete follow-up are PLP
# population settings (originalProductionPopulationSettings()).
#
# CirceR is not in renv; run outside the project folder with the user library:
#   cd /tmp && Rscript --no-init-file <repo>/extras/production/buildOriginalCohorts.R <repo>

args <- commandArgs(trailingOnly = TRUE)
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
repo <- normalizePath(args[1])
definitions <- jsonlite::fromJSON(file.path(repo, "inst/original/definitions_v1.json"), simplifyVector = FALSE)$definitions
tasks <- utils::read.csv(file.path(repo, "inst/settings/originalTasks.csv"))
criteriaName <- c(condition_occurrence = "ConditionOccurrence", procedure_occurrence = "ProcedureOccurrence")

conceptSet <- function(name) {
  ids <- sort(as.numeric(unlist(definitions[[name]]$concept_ids)))
  list(list(id = 0, name = name, expression = list(items = lapply(ids, function(id) list(
    concept = list(CONCEPT_ID = id), isExcluded = FALSE, includeDescendants = FALSE, includeMapped = FALSE)))))
}
cohort <- function(conceptSets, criteria, priorDays, primaryLimit, qualifiedLimit, inclusionRules = list(),
                   endStrategy = NULL) {
  x <- list(
    cdmVersionRange = ">=5.0.0",
    PrimaryCriteria = list(CriteriaList = list(criteria), ObservationWindow = list(PriorDays = priorDays, PostDays = 0),
                           PrimaryCriteriaLimit = list(Type = primaryLimit)),
    ConceptSets = conceptSets, QualifiedLimit = list(Type = qualifiedLimit), ExpressionLimit = list(Type = qualifiedLimit),
    InclusionRules = inclusionRules, CensoringCriteria = list(), CollapseSettings = list(CollapseType = "ERA", EraPad = 0),
    CensorWindow = structure(list(), names = character()))
  if (!is.null(endStrategy)) x$EndStrategy <- endStrategy
  x
}
diseaseTarget <- function(name) {
  domain <- criteriaName[[definitions[[name]]$domain_table]]
  criteria <- structure(list(list(CodesetId = 0, Age = list(Value = 18, Op = "gte"))), names = domain)
  cohort(conceptSet(name), criteria, priorDays = 365, primaryLimit = "First", qualifiedLimit = "First")
}
adultsTarget <- function() {
  # Circe measures PeriodLength on the adjusted period (start + 365 to end): >= 0 is a period of >= 365 days
  criteria <- list(ObservationPeriod = list(
    PeriodLength = list(Value = 0, Op = "gte"),
    DateAdjustment = list(StartWith = "START_DATE", StartOffset = 365, EndWith = "END_DATE", EndOffset = 0)))
  age <- list(name = "aged 18 or over at entry", expression = list(
    Type = "ALL", CriteriaList = list(), DemographicCriteriaList = list(list(Age = list(Value = 18, Op = "gte"))),
    Groups = list()))
  cohort(list(), criteria, priorDays = 0, primaryLimit = "All", qualifiedLimit = "First", inclusionRules = list(age))
}
outcome <- function(name) {
  domain <- criteriaName[[definitions[[name]]$domain_table]]
  criteria <- structure(list(list(CodesetId = 0)), names = domain)
  cohort(conceptSet(name), criteria, priorDays = 0, primaryLimit = "All", qualifiedLimit = "All",
         endStrategy = list(DateOffset = list(DateField = "StartDate", Offset = 0)))
}

targets <- sort(unique(tasks$target)); outcomes <- sort(unique(tasks$outcome))
cohorts <- rbind(
  data.frame(cohortId = 1000 + seq_along(targets), phenotype = targets, role = "target"),
  data.frame(cohortId = 2000 + seq_along(outcomes), phenotype = outcomes, role = "outcome"))
cohorts$cohortName <- paste(cohorts$phenotype, cohorts$role)
dir.create(file.path(repo, "inst/cohorts"), showWarnings = FALSE)
for (i in seq_len(nrow(cohorts))) {
  p <- cohorts$phenotype[i]
  expression <- if (cohorts$role[i] == "outcome") outcome(p) else if (p == "adults") adultsTarget() else diseaseTarget(p)
  json <- jsonlite::toJSON(expression, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  sql <- CirceR::buildCohortQuery(CirceR::cohortExpressionFromJson(as.character(json)),
                                  CirceR::createGenerateOptions(generateStats = FALSE))
  writeLines(json, file.path(repo, "inst/cohorts", sprintf("%d.json", cohorts$cohortId[i])))
  writeLines(sql, file.path(repo, "inst/cohorts", sprintf("%d.sql", cohorts$cohortId[i])))
}
utils::write.csv(cohorts[, c("cohortId", "cohortName", "phenotype", "role")], file.path(repo, "inst/settings/cohorts.csv"),
                 row.names = FALSE)
message("CirceR ", packageVersion("CirceR"), ": ", nrow(cohorts), " cohorts written")
