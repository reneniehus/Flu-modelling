# where is the SARI data?
#NSQL3

source('db/logger.R')
source('db/sql_utils.R')
# Logger is needed for running SQL utils
logger <- forge_logger()(logLevel = 'INFO',
                         fileFormat = "csv",
                         memoryMonitor = T)
logger$info("Logger has been initialised.")
# This is where the SQL Profiles are stored (always needed)
dbDir <- 'db/'

data <- read_data(table = 'clean.RESPISEVERE_Haggregated', # directly defining the data table
               connInfo = 'pop' # using pop.PROFILE (defining server and the database)
) 


pop_data <- read_data(table = 'out.DM_Population_ByCountryEU', # directly defining the data table
                      connInfo = 'pop' # using pop.PROFILE (defining server and the database)
) 
