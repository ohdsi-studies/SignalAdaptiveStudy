# Settings fingerprints. Each cached or finished output records the settings it was made with: a JSON file with the
# hash of all settings together, the hash of each one, and the settings themselves (readable). On a restart the current
# settings are compared with it: a data cache made with other settings is extracted again, counts are recounted, and a
# finished run made with other settings stops execute() with an error, so that results of two designs never mix
# silently. Outputs from before fingerprints have no file and are accepted (noted in the log). Package versions are
# not part of the fingerprint.

FINGERPRINT_FILE <- "fingerprint.json"

# settings: a named list.
fingerprint <- function(settings) {
  list(hash = rlang::hash(settings), fields = lapply(settings, rlang::hash), settings = settings)
}

writeFingerprint <- function(fp, file) ParallelLogger::saveSettingsToJson(fp, file)

# NULL when there is no fingerprint file (an output from before fingerprints); otherwise the names of the settings that
# differ from those recorded (character(0): the same settings).
fingerprintChanges <- function(fp, file) {
  if (!file.exists(file)) return(NULL)
  old <- ParallelLogger::loadSettingsFromJson(file)
  if (identical(as.character(old$hash), fp$hash)) return(character())
  fields <- union(names(fp$fields), names(old$fields))
  changed <- fields[vapply(fields, function(f) !identical(as.character(old$fields[[f]]), as.character(fp$fields[[f]])),
                           logical(1))]
  if (length(changed)) changed else "(all)"
}

# The data a task's outputs depend on (NULL for a data source that cannot say, e.g. an in-memory plpData).
taskDataFingerprint <- function(dataSource, task) {
  if (is.null(dataSource$dataFingerprint)) NULL else dataSource$dataFingerprint(task)
}

splitFingerprint <- function(dataSource) {
  list(label = dataSource$splitLabel, seeds = if (identical(dataSource$splitLabel, "study")) STUDY_SEEDS)
}

# Counts: the cohorts and population rules, not the covariates or the model sample.
countsFingerprint <- function(dataSource, task) {
  data <- taskDataFingerprint(dataSource, task)
  fingerprint(list(data = data[setdiff(names(data), c("covariates", "sample"))],
                   population = dataSource$populationSettings(task), split = splitFingerprint(dataSource)))
}

runFingerprint <- function(dataSource, task, modelSettings, budget, ordering) {
  fingerprint(list(data = taskDataFingerprint(dataSource, task), population = dataSource$populationSettings(task),
                   split = c(splitFingerprint(dataSource), list(budget = budget, ordering = ordering)),
                   preprocessing = studyPreprocessSettings(), method = modelSettings))
}
