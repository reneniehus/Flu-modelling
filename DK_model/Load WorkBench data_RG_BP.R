#install.packages("RODBC")
library(RODBC)

# Some databases have snapshots (INFL, NCOV, AMC, AMR, FWD), and some don't. Use the snapshot, if it exists.
# Contact your data manager if you have questions.

# VPD.server.and.database <- "Driver=SQL Server;Server=nsql3;Database=VPD"
# VPD.connection <- odbcDriverConnect(VPD.server.and.database)

NCOV_Aggr_Snapshot <- "Driver=SQL Server;Server=nsql3;Database=NCOVAGGR_Snapshot"
NCOV_Aggr.connection <- odbcDriverConnect(NCOV_Aggr_Snapshot)

NCOV_Case_Snapshot <- "Driver=SQL Server;Server=nsql3;Database=NCOV_Snapshot"
NCOV_Case.connection <- odbcDriverConnect(NCOV_Case_Snapshot)


INFL_Snapshot <- "Driver=SQL Server;Server=nsql3;Database=INFL"
INFL.connection <- odbcDriverConnect(INFL_Snapshot)


# list all the databases on the server
sql.code <- "SELECT name FROM sys.databases"
df.workbench.databases <- sqlQuery(INFL.connection, sql.code)

# list of all tables from the database
sql.code <- "SELECT s.name AS Schema_Name, t.name AS Table_Name FROM sys.tables T INNER JOIN sys.schemas S ON S.schema_id = T.schema_id WHERE type='U'"
df.tables <- sqlQuery(INFL.connection, sql.code)

#sql.code <- "SELECT TOP (100000) * FROM [clean].[NCOV]"
# sql.code <- "SELECT TOP (100) * FROM [INFL].[clean].[INFLCLIN_Haggregated] ORDER BY variable"
sql.code <- "SELECT * FROM [INFL].[clean].[INFLCLIN_Haggregated]"
df.inf.clin <- sqlQuery(INFL.connection, sql.code)

#














#######
# load the cleand NCOV aggregated data (from the "NCOVAGGR_Snapshot" database, "NCOVAGGR" table, in the "clean" schema, and filter on the server for 2020 data)
sql.code <- "SELECT * FROM [NCOVAGGR_Snapshot].[clean].[NCOVAGGR] WHERE DateUsedForStatisticsYear > 2020 AND ReportingCountry = 'SI'"
df.NCOV_Aggr <- sqlQuery(NCOV_Aggr.connection, sql.code)

sql.code <- "SELECT TOP (100000) * FROM [clean].[NCOV]"
df.NCOV_Case <- sqlQuery(NCOV_Case.connection, sql.code)
size1000rows <- object.size(df.NCOV_Case)

sql.code <- "SELECT count(*) FROM [clean].[NCOV]"
nr_rows.NCOV_Case <- sqlQuery(NCOV_Case.connection, sql.code)

as.numeric(nr_rows.NCOV_Case)*as.numeric(size1000rows)/100000/1024/1024/1024

sql.code <- "SELECT * FROM [INFL_Snapshot].[clean].[INFLVIR_Haggregated] WHERE DateUsedForStatisticsYear >= 2015"
df.INFLVIR <- sqlQuery(INFL.connection, sql.code)



sql.code <- "SELECT TOP (100) * FROM [INFL_Snapshot].[clean].[INFLSARI_Case]"
# sql.code <- "SELECT * FROM [INFL_Snapshot].[clean].[SARI_Case] WHERE DateUsedForStatisticsYear >= 2015"
df.SAR_Case <- sqlQuery(INFL.connection, sql.code)



sql.code <- "SELECT * FROM [INFL_Snapshot].[clean].[DataSource]"
df.DataSource <- sqlQuery(INFL.connection, sql.code)

sql.code <- "SELECT * FROM [INFL_Snapshot].[clean].[ReportingPeriod]"
df.ReportingPeriod <- sqlQuery(INFL.connection, sql.code)


sql.code <- "SELECT * FROM [NCOVAGGR_Snapshot].[out].[Log_WorkbenchRefresh]"
df.workbench.refresh.date <- sqlQuery(INFL.connection, sql.code)
# Please make sure you're loading the table from the same database as your actual data.
# LogOn is when the refresh actually happened. DataAsOf is the time when the data was extracted from the Data Warehouse. 
# So if LogOn says yesterday, but DataAsOf says September last year, it means you're seeing the data as it existed in TESSy on that date, last year. 
# We can always make a "back-in-time" extraction, so we can "freeze" the data retroactively on any date, at any time. DataAsOf is always the important one.

