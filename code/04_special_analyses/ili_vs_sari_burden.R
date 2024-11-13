# Eeva Broling knows about data set before 2020 on 
# Leah: in 2020 a lot of new SARI systems and thus TESSY tables were started due
# the surveillance needs for COVID

# pre 2020, there was less coherence between countries for hospital/severe respiratory care reporting

# we are looking for severe respiratory cases
# hospitalised respiratory case counts
# sentinel/non-sentinel, we consider both
# age-disaggregation, useful to have not a requirement
# pathogen typing, great to have but not needed

# most of that data should be in NSQL3

source('db/logger.R')
source('db/sql_utils.R')
# Logger is needed for running SQL utils
logger <- forge_logger()(logLevel = 'INFO',
                         fileFormat = "csv",
                         memoryMonitor = T)
logger$info("Logger has been initialised.")
# This is where the SQL Profiles are stored (always needed)
dbDir <- 'db/'

# RESPISEVERE_Haggregated
data <- read_data(table = 'RESPISURV.clean.RESPISEVERE_Haggregated', # directly defining the data table
               connInfo = 'pop' # using pop.PROFILE (defining server and the database)
) 
data$DateUsedForStatisticsYear %>% range() # 2020-2024

# INFLSARI
data <- read_data(table = 'RESPISURV.clean.RESPISURV_Case', # directly defining the data table
                  connInfo = 'NSQL3' # using pop.PROFILE (defining server and the database)
) 


# INFLSARI ? Leah will look into where 
data <- read_data(table = 'INFLSARI', # directly defining the data table
                  connInfo = 'NSQL3' # using pop.PROFILE (defining server and the database)
) 

# use the following query in SQL to find data tables

# USE master;
# GO
# EXEC sp_MSforeachdb 
# 'USE [?]; 
#   SELECT TABLE_SCHEMA, TABLE_NAME, ''?'' AS DatabaseName 
#   FROM INFORMATION_SCHEMA.TABLES 
#   WHERE TABLE_NAME COLLATE SQL_Latin1_General_CP1_CI_AS LIKE ''%SARI%'' COLLATE SQL_Latin1_General_CP1_CI_AS;';



 

data$DateUsedForStatisticsYear %>% range() # 2020 2024
