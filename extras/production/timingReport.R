# Where the model runs (3c) spend their time: one row per finished run with its timings, sizes and search steps, from
# the run folders (result.csv, the saved model's trainDetails.json and the file times). Read-only, so it can run
# while the study runs. Writes timing_<database>.csv in the output folder: aggregate only (no patient-level data),
# fine to send.
#
#   Rscript extras/production/timingReport.R [outputFolder]     (default: the outputFolder of extras/CodeToRun.R)

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) outputFolder <- args[1] else {
  source("extras/localSettings.R")
  outputFolder <- file.path(getwd(), "SignalAdaptiveStudy_output", databaseName)
}
files <- list.files(file.path(outputFolder, "models"), pattern = "^result\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(files) == 0) stop("no finished runs in ", outputFolder)

runs <- do.call(rbind, lapply(files, function(file) {
  r <- utils::read.csv(file, stringsAsFactors = FALSE)
  model <- file.path(dirname(file), "model")
  details <- tryCatch(ParallelLogger::loadSettingsFromJson(file.path(model, "trainDetails.json")), error = function(e) NULL)
  rows <- function(x) if (is.null(x)) NA else NROW(as.data.frame(x))
  covariates <- tryCatch(nrow(utils::read.csv(file.path(model, "covariateImportance.csv"))), error = function(e) NA)
  data.frame(r[, intersect(c("database", "taskId", "budget", "ordering", "method", "nTrain", "eventsTrain", "nonZero",
                             "fitSeconds", "runSeconds"), names(r))],
             covariates = covariates, searchSteps = rows(details$searchTrace), hyperParameters = rows(details$hyperParamSearch),
             finished = file.mtime(file), stringsAsFactors = FALSE)
}))
runs <- runs[order(runs$finished), ]
# Time between the end of the previous run and the start of this one: loading the task's plpData, the population and
# the split (and, at a task's first run, extracting the data if it was not cached). A gap over 4 hours is taken to be
# a restart (R not running) and left out.
runs$gapSeconds <- c(NA, pmax(0, diff(as.numeric(runs$finished)) - runs$runSeconds[-1]))
restarts <- sum(runs$gapSeconds > 4 * 3600, na.rm = TRUE)
runs$gapSeconds[runs$gapSeconds > 4 * 3600] <- NA
runs$plpOverheadSeconds <- runs$runSeconds - runs$fitSeconds # PLP around the fit: preprocessing, test predictions
runs$finished <- format(runs$finished, "%Y-%m-%d %H:%M:%S")
file <- file.path(outputFolder, sprintf("timing_%s.csv", runs$database[1]))
utils::write.csv(runs, file, row.names = FALSE)

hours <- function(s) round(sum(s, na.rm = TRUE) / 3600, 2)
cat(sprintf("%d runs, %.1f h of runs, %.1f h of gaps between runs (%d gap(s) over 4 h left out as restarts)\n\n",
            nrow(runs), hours(runs$runSeconds), hours(runs$gapSeconds), restarts))
cat("Hours by method (fit / PLP overhead):\n")
print(do.call(rbind, lapply(split(runs, runs$method), function(d)
  data.frame(method = d$method[1], runs = nrow(d), fitHours = hours(d$fitSeconds), overheadHours = hours(d$plpOverheadSeconds),
             medianFitMin = round(stats::median(d$fitSeconds) / 60, 1)))), row.names = FALSE)
cat("\nHours by task:\n")
print(do.call(rbind, lapply(split(runs, runs$taskId), function(d)
  data.frame(taskId = d$taskId[1], runs = nrow(d), maxTrain = max(d$nTrain), runHours = hours(d$runSeconds),
             gapHours = hours(d$gapSeconds)))), row.names = FALSE)
message("written: ", file)
