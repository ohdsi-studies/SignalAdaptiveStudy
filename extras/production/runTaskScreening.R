# Screens the candidate tasks of the follow-up study in one database: cohorts and counts only, no covariates and no models.
# Writes a zip of aggregates (counts below minCellCount censored) to send back. Safe next to a running study: it uses its own
# cohort table, sar_screen_cohort_<database>. Restartable: rerun after an interruption.
#
#   Rscript extras/production/runTaskScreening.R               (in the project folder; settings as for CodeToRun.R)

# Needed by the Databricks JDBC driver, before Java starts (see CodeToRun.R)
options(java.parameters = "--add-opens=java.base/java.nio=org.apache.arrow.memory.core,ALL-UNNAMED")
library(SignalAdaptiveStudy)

server <- Sys.getenv("DATABASE_SERVER")
connectionString <- if (startsWith(server, "jdbc:")) server else
  sprintf("jdbc:databricks://%s:443/default;transportMode=http;ssl=1;AuthMech=3;httpPath=%s;", server, Sys.getenv("DATABASE_PATH"))
if (!grepl("EnableArrow=", connectionString, ignore.case = TRUE)) {
  connectionString <- paste0(sub(";?$", ";", connectionString), "EnableArrow=1;")
}
connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "spark", connectionString = connectionString, user = "token", password = Sys.getenv("DATABASE_PASSWORD"),
  pathToDriver = Sys.getenv("DATABASE_DRIVER"))

source("extras/localSettings.R")   # scratchSchema, minCellCount, cdmSchemas, databaseName
options(sqlRenderTempEmulationSchema = scratchSchema)
runTaskScreening(connectionDetails, cdmSchemas[[databaseName]], scratchSchema, databaseName,
                 outputFolder = file.path(getwd(), "SignalAdaptiveStudy_screening", databaseName), minCellCount = minCellCount)
