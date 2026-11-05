# Cohorts and tasks of the PLP benchmark run. Copies the PLP benchmark definitions as they are from the PLPBenchmarks
# package (cohort ID 3000 + the benchmark's ID), and the 64 phenotype features of the PLP vignette "Constrained Predictors" from
# PhenotypeLibrary (cohort ID 100000 + the library's ID). Writes inst/cohorts/<cohortId>.json and .sql (CirceR), adds them to
# inst/settings/cohorts.csv, and writes inst/settings/benchmarkTasks.csv from the benchmark's problem specification and model designs.
#
# CirceR and PhenotypeLibrary are not in renv; run outside the project folder with the user library:
#   cd /tmp && Rscript --no-init-file <repo>/extras/production/addBenchmarkCohorts.R <repo> <PLPBenchmarks checkout>

args <- commandArgs(trailingOnly = TRUE)
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
repo <- normalizePath(args[1]); benchmarks <- normalizePath(args[2])
problems <- utils::read.csv(file.path(benchmarks, "inst/extdata/ProblemSpecification.csv"), stringsAsFactors = FALSE)
designs <- readRDS(file.path(benchmarks, "inst/extdata/BenchmarkModelDesignList.Rds"))
stopifnot(length(designs) == nrow(problems),
          all(vapply(designs, function(d) d$targetId, 1) == problems$targetId),
          all(vapply(designs, function(d) d$outcomeId, 1) == problems$outcomeId),
          all(vapply(designs, function(d) d$populationSettings$riskWindowStart, 1) == problems$TAR_start_day),
          all(vapply(designs, function(d) d$populationSettings$riskWindowEnd, 1) == problems$TAR_end_day))
names <- utils::read.csv(file.path(benchmarks, "inst/settings/CohortsToCreate.csv"))
nameOf <- function(id) names$cohort_name[match(id, names$cohort_id)]
write <- function(cohortId, json) {
  sql <- CirceR::buildCohortQuery(CirceR::cohortExpressionFromJson(json), CirceR::createGenerateOptions(generateStats = FALSE))
  writeLines(json, file.path(repo, "inst/cohorts", sprintf("%d.json", cohortId)))
  writeLines(sql, file.path(repo, "inst/cohorts", sprintf("%d.sql", cohortId)))
}
cohorts <- utils::read.csv(file.path(repo, "inst/settings/cohorts.csv"))
ids <- sort(unique(c(problems$targetId, problems$outcomeId)))
role <- ifelse(ids %in% problems$targetId, "target", "outcome")
new <- data.frame(cohortId = 3000 + ids, cohortName = paste(nameOf(ids), role), phenotype = nameOf(ids), role = role, use = "benchmark")
new <- new[!new$cohortId %in% cohorts$cohortId, ]
for (id in ids) write(3000L + id, paste(readLines(file.path(benchmarks, "inst/cohorts", paste0(id, ".json"))), collapse = "\n"))

library <- as.data.frame(PhenotypeLibrary::getPlCohortDefinitionSet(1152:1215))
stopifnot(nrow(library) == 64)
for (i in seq_len(nrow(library))) write(100000L + library$cohortId[i], library$json[i])
features <- data.frame(cohortId = 100000 + library$cohortId, cohortName = library$cohortName, phenotype = library$cohortName, role = "feature", use = "feature")
cohorts <- rbind(cohorts[cohorts$use == "study" | cohorts$cohortId %in% c(), ], new, features)
utils::write.csv(cohorts[order(cohorts$cohortId), ], file.path(repo, "inst/settings/cohorts.csv"), row.names = FALSE)

tasks <- data.frame(taskId = paste0("benchmark_", problems$analysisName), benchmark = problems$analysisName, target = nameOf(problems$targetId),
                    outcome = nameOf(problems$outcomeId), timeAtRiskDays = problems$TAR_end_day, riskWindowStart = problems$TAR_start_day,
                    firstExposureOnly = vapply(designs, function(d) d$populationSettings$firstExposureOnly, TRUE),
                    targetId = 3000 + problems$targetId, outcomeId = 3000 + problems$outcomeId, populationRule = "benchmark")
utils::write.csv(tasks, file.path(repo, "inst/settings/benchmarkTasks.csv"), row.names = FALSE)
commit <- tryCatch(system2("git", c("-C", benchmarks, "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) NA)
message("CirceR ", packageVersion("CirceR"), "; PLPBenchmarks ", commit, "; PhenotypeLibrary ", packageVersion("PhenotypeLibrary"), ": ", nrow(new), " benchmark and ", nrow(features), " feature cohorts added, ", nrow(tasks), " tasks")
