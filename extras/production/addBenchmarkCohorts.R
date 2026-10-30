# Cohorts of the new catalogue tasks (3b counts only, so the tasks stay unseen): the PLP benchmark definitions of
# SAR-017, SAR-175, SAR-287 and SAR-288 (investigator, 2026-10-30: use the benchmark definitions as they are), copied
# from the PLPBenchmarks package with cohort ID 3000 + the benchmark's ID. Writes inst/cohorts/<cohortId>.json and
# .sql (CirceR) and adds them to inst/settings/cohorts.csv.
#
# CirceR is not in renv; run outside the project folder with the user library:
#   cd /tmp && Rscript --no-init-file <repo>/extras/production/addBenchmarkCohorts.R <repo> <PLPBenchmarks checkout>

args <- commandArgs(trailingOnly = TRUE)
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
repo <- normalizePath(args[1]); benchmarks <- normalizePath(args[2])
benchmarkIds <- c(2, 1, 3, 4, 17, 25, 24, 18) # T2dm_second_line, HF, anemia, colorectal cancer, afib, non-hem-stroke, mdd, bipolar
names <- utils::read.csv(file.path(benchmarks, "inst/settings/CohortsToCreate.csv"))
roles <- c("target", "outcome")[c(1, 2, 1, 2, 1, 2, 1, 2)]
new <- data.frame(cohortId = 3000 + benchmarkIds,
                  cohortName = paste(names$cohort_name[match(benchmarkIds, names$cohort_id)], roles),
                  phenotype = names$cohort_name[match(benchmarkIds, names$cohort_id)], role = roles)
for (i in seq_len(nrow(new))) {
  json <- paste(readLines(file.path(benchmarks, "inst/cohorts", paste0(benchmarkIds[i], ".json"))), collapse = "\n")
  sql <- CirceR::buildCohortQuery(CirceR::cohortExpressionFromJson(json), CirceR::createGenerateOptions(generateStats = FALSE))
  writeLines(json, file.path(repo, "inst/cohorts", sprintf("%d.json", new$cohortId[i])))
  writeLines(sql, file.path(repo, "inst/cohorts", sprintf("%d.sql", new$cohortId[i])))
}
cohorts <- utils::read.csv(file.path(repo, "inst/settings/cohorts.csv"))
cohorts <- rbind(cohorts[!cohorts$cohortId %in% new$cohortId, ], new)
utils::write.csv(cohorts, file.path(repo, "inst/settings/cohorts.csv"), row.names = FALSE)
commit <- tryCatch(system2("git", c("-C", benchmarks, "rev-parse", "--short", "HEAD"), stdout = TRUE), error = function(e) NA)
message("CirceR ", packageVersion("CirceR"), "; PLPBenchmarks ", commit, ": ", nrow(new), " cohorts added")
