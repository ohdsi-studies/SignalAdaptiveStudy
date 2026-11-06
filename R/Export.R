#' Export aggregate results
#'
#' Writes the result rows and counts of [execute()] as CSV files, with the step timings, and zips them into one file, with counts below
#' `minCellCount` replaced by `-minCellCount`. Only aggregates leave: no patient-level data, predictions, models or
#' coefficients.
#'
#' @param outputFolder The `outputFolder` of [execute()].
#' @param minCellCount Smallest count that may be shown.
#' @param zipFile Zip file to write.
#' @return The path of the zip file.
#' @export
exportResults <- function(outputFolder, minCellCount = 5,
                          zipFile = file.path(outputFolder, paste0("SignalAdaptiveStudy_results_", format(Sys.Date(), "%Y%m%d"), ".zip"))) {
  folder <- tempfile("sar_export_")
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  censor <- function(d, columns) {
    for (column in intersect(columns, names(d))) {
      small <- !is.na(d[[column]]) & d[[column]] > 0 & d[[column]] < minCellCount
      d[[column]][small] <- -minCellCount
    }
    d
  }
  results <- collectResults(outputFolder)
  utils::write.csv(censor(results, c("nTrain", "eventsTrain", "nTest", "eventsTest")), file.path(folder, "results.csv"), row.names = FALSE)
  countFiles <- list.files(file.path(outputFolder, "counts"), pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  counts <- if (length(countFiles)) do.call(rbind, lapply(countFiles, utils::read.csv, stringsAsFactors = FALSE)) else data.frame()
  utils::write.csv(censor(counts, c("n", "poolN", "poolEvents", "testN", "testEvents")), file.path(folder, "counts.csv"), row.names = FALSE)
  # The duration of each step (cohorts, extraction, cache, population, split, preparation, runs): task names and seconds only
  stepFiles <- list.files(outputFolder, pattern = "^steps_.*\\.csv$", full.names = TRUE)
  if (length(stepFiles)) file.copy(stepFiles, folder)
  writeLines(c(paste("Exported", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), packageVersions(),
               paste("minCellCount", minCellCount, "(smaller counts shown as", -minCellCount, ")")),
             file.path(folder, "export.txt"))
  zipFile <- normalizePath(zipFile, mustWork = FALSE)
  dir.create(dirname(zipFile), recursive = TRUE, showWarnings = FALSE)
  unlink(zipFile)
  zip::zipr(zipFile, list.files(folder, full.names = TRUE))
  zipFile
}
