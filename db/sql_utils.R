#' @title Create SQL connection
#' 
#' @description Returns a db connection object using the specified SQL 
#'   credentials.
#' 
#' @param connInfo string. Name of the profile to be used on this connection.
#' @param read_only boolean. Flag used to optimise the SQL connection when only
#'   reading operations are going to be performed.
#' 
#' @return connection object. An object used to operate ODBC functions such as
#'   RODBC::sqlQuery().
#'    
#' @note This function is intended to be call from **within** the SQL utils
#'   functions, and never outside. Bear in mind that if you do, you'll need to
#'   manually close the connection once you're done
create_sql_connection <- function(connInfo, read_only = F){
  # Load the database profile
  connInfo <- get_database_info(db = connInfo)
  
  # If credentials are not null, they are included within the connString.
  if(length(connInfo$credentials) > 0){
    credString <- glue::glue(";uid={connInfo$credentials$uid};pwd={connInfo$credentials$pwd}")
  } else {
    credString <- ""
  }
  
  connInfo$port <- if(connInfo$port != "") str_c(',', connInfo$port) else ""
  # Write the connection string and set up the connection object
  connString <- glue::glue("Driver={connInfo$sqlDriver};Server={connInfo$server}{connInfo$port};Database={connInfo$database}{credString}")
  
  connection <- RODBC::odbcDriverConnect(connection = connString, DBMSencoding = 'UTF8', believeNRows = F,
                                         rows_at_time = 1, readOnlyOptimize = read_only)
  logger$debug(glue::glue("Connected to database {connInfo$database} on {connInfo$server}."))
  
  return(connection)
}

#' @title Read data from a SQL database
#' 
#' @description Returns a dataframe from a SQL database, or raises an error
#'   if the corresponding server fails to execute the given instruction.
#'   
#' @param query character. String containing a full SELECT query. 
#'   It overrides the table, columns, filter and object_key parameters. Additionally,
#'   it can also be the name of a SQL template file. See parameter path for more 
#'   information. 
#' @param table character. Name of the table to retrieve data from.
#' @param columns character|list. A list of columns to retrieve, or a single string  
#'   containing the columns for the SELECT clause. Default value is "*", all
#'   available columns.
#' @param filters character|list. A string containing the WHERE clause of the SQL
#'   query, or a list of criteria. Object key filter is added afterwards. 
#' @param connInfo string. Name of the profile to be used on this connection.
#' @param additionalParams list. Named list of parameters to be passed to the template
#'   query.
#' @param forceTypeConversion logical. Whether to make R force the types conversion, 
#'   overriding those set by default by RODBC, or leave them as they are. Can be useful
#'   sometimes.
#' @param path character. Path to the directory where queries templates are stored.
#'   Default value is 'templates/sql', but can be overridden without need to specify it
#'   on every call by setting the parameter SQLTemplatePath in .GlobalEnv.
#' 
#' @return table dataFrame. The dataframe containing the data requested from the table.
read_data <- function(query = character(), 
                      table, 
                      columns = "*", 
                      filters = "", 
                      connInfo,
                      additionalParams = list(),
                      forceTypeConversion = F,
                      path = 'templates/sql/'){
  # Assign the base template path in case SQLTemplatePath is not set on .GlobalEnv
  if(exists(x = "SQLTemplatePath", envir = .GlobalEnv)){
    path <- get(x = 'SQLTemplatePath', envir = .GlobalEnv)
  }
  # Check if the user has provided a query or the path to a query template
  if(length(query) != 0){
    if(file.exists(file.path(path, query))){
      # Query is provided in a template
      query <- glue::glue(paste0(readLines(file.path(path, query)), collapse = '\n'))
      logger$debug(query)
    } 
    if(startsWith(query, 'SELECT') | startsWith(query, 'WITH')){ 
      # Query starts with a SELECT or a WITH clause. @WIP we need to sanitize this to prevent INSERT, UPDATE and DELETE
      query <- query
    } else {
      # Query was either invalid or not found.
      stop(glue::glue("Input {query} was not a valid template or SELECT query."))
    }
  } else {
    # Use theq uery builder
    query <- build_sql_query(columns = columns,
                             table = table,
                             filters = filters,
                             sanitizeColumns = stringr::str_detect(string = get_database_info(connInfo)$sqlDriver, 
                                                                   pattern = 'PostgreSQL'))
    logger$debug(query)
  }
  
  withCallingHandlers({
    # Request a connection object
    con <- create_sql_connection(connInfo = connInfo, read_only = T)
    # Execute the query 
    data <- RODBC::sqlQuery(channel = con, query = query, as.is = TRUE, rows_at_time = 1) 
    # Health check
    if(!is.null(data)){
      if(!is.data.frame(data)){
        # Errors are collapsed for better reading
        stop(paste0(data, collapse = '; '))
      } else {
        logger$debug("Query executed correctly.")
        logger$debug(glue::glue("{nrow(data)} lines fetched with {ncol(data)} columns."))
      }
    } else {
      logger$warning("The query did not return a response.")
    }
  },
  error = function(cond){
    logger$error("There was an error when querying the database.")
    logger$error(cond$message)
    if(!stringr::str_detect(string = get_database_info(connInfo)$sqlDriver, pattern = 'SQLite')){
      RODBC::odbcClose(con)
      logger$debug("Connection to server closed.")
    }
    logger$debug("Connection to server closed.")
  }) 
  
  # Apparently SQLite connections don't need to be closed...
  if(!stringr::str_detect(string = get_database_info(connInfo)$sqlDriver, pattern = 'SQLite')){
    RODBC::odbcClose(con)
    logger$debug("Connection to server closed.")
  }
  
  # Restore colnames names for PostgreSQL databases
  if(stringr::str_detect(string = get_database_info(connInfo)$sqlDriver, pattern = 'PostgreSQL')){
    colnames(data) <- stringi::stri_replace_all_fixed(str = colnames(data), 
                                                      pattern = colnames_sanitizer$replacement, 
                                                      replacement = colnames_sanitizer$symbol, 
                                                      vectorize_all = F)
  }
  
  # Assure R returns the right types of data
  if(forceTypeConversion){
    data <- readr::type_convert(data)
  }
  return(data)
}

#' @title Run an arbitrary query against a SQL database
#' 
#' @description Runs any query against a SQL database. Caution is advised.
#'   
#' @param query character. Query to be run.
#' @param conn_info string. Name of the profile to be used on this connection.
#' @param args list. List of arguments in case of using a templatable query.
#' 
#' @note This function requires a working logger object in the upper
#' environment.
run_query = function(query,
                     conn_info,
                     additionalParams = list(),
                     path = 'templates/sql/'){
  # Assign the base template path in case SQLTemplatePath is not set on .GlobalEnv
  if(exists(x = "SQLTemplatePath", envir = .GlobalEnv)){
    path <- get(x = 'SQLTemplatePath', envir = .GlobalEnv)
  }
  withCallingHandlers({
    if(length(query) != 0){
      if(file.exists(file.path(path, query))){
        # Query is provided in a template
        query <- glue::glue(paste0(readLines(file.path(path, query)), collapse = '\n'))
        logger$debug(query)
      } 
    } else {
      stop("No valid query was provided.")
    }
    
    con = create_sql_connection(conn_info)
    res <- sqlQuery(con, glue(query), as.is = TRUE) 
    if(!is.null(res)) if(any(grepl("ERROR", res))){
      stop(paste0(res, collapse = ';'))
    }	
    odbcClose(con)
    logger$debug("Connection to server closed.")
    
  },
  error = function(cond){
    logger$debug(glue(query))
    logger$error("There was an error when querying the database.")
    logger$error(cond$message)
    odbcClose(con)
    logger$debug("Connection to server closed.")
  }) 
  return(res)
}

#' @title Write data to a SQL database
#' 
#' @description Writes a dataset to a SQL table, giving the possibility of
#'   creating it if it doesn't exist.
#'   
#' @param data dataframe. Dataframe containing the data to upload.
#' @param table character. Name of the table to write data to.
#' @param schema character. Name of the schema
#' @param connInfo string. Name of the profile to be used on this connection.
#' @param autoCreate bool. Flag to let SQL Utils create the table if it does
#'   not exist. It will autodetect the types and use approppriate ones.
#'   Currently only PostgreSQL is supported.
#' @param deletePreviousRecords bool. Deletes the previous records wiping the
#'   table clean in the process.
#' @param deleteMatchingRecords named list. When specified, will delete the matching
#'   records from the specified table.
#'
#' @usage write_data(data = myDataset, table = 'dbo.COVIDCases', connectionInfo = 'NSQL3_COVID', 
#'                   deletePreviousRecords = T, deleteMatchinRecords = list(DateUsedForStatisticsISO = "> '2023-W01'"))
#' 
#' @note Make sure the deleteMatchingRecords parameters is set to the right values, respecting SQL formats.
write_data <- function(data, 
                       table, 
                       schema = 'dbo',
                       connInfo,
                       autoCreate = F,
                       deletePreviousRecords = F,
                       deleteMatchingRecords = character()){
  withCallingHandlers({
    # Request a connection object
    con <- create_sql_connection(connInfo = connInfo)
    
    # Table creation
    # Table creation
    if(autoCreate) {
      query = create_new_table(table = table,
                               schema = schema,
                               data = data, 
                               flavor = get_database_info(connInfo)$sqlDriver)
      logger$debug(query)
      res <- RODBC::sqlQuery(con, query, as.is = T)
      logger$debug(res)
      if(!identical(res, character(0)) && grep('ERROR: ', res)){
        stop(paste0(res, collapse = '; '))
      }
    }
    
    # Colnames sanitization
    if( stringr::str_detect(get_database_info(connInfo)$sqlDriver, 'PostgreSQL')){
      colnames(data) <- stringi::stri_replace_all_fixed(str = colnames(data), 
                                                        pattern = colnames_sanitizer$replacement, 
                                                        replacement = colnames_sanitizer$symbol, 
                                                        vectorize_all = F)
    }
    
    # String sanitization (common for all databases)
    if(stringr::str_detect(get_database_info(connInfo)$sqlDriver, 'SQL Server')){
      data = dplyr::ungroup(data) |> 
        dplyr::mutate(
          dplyr::across(everything(), 
                        ~ stringr::str_replace_all(., "'", "''"))) |>
        dplyr::mutate_if(check_numeric_types, ~stringr::str_c("N'", .x, "'")) |>
        dplyr::mutate(
          dplyr::across(
            dplyr::everything(),
            ~ tidyr::replace_na(.,'NULL')
          )
        )
      
    } else {
      data = dplyr::ungroup(data) |> 
        dplyr::mutate(
          dplyr::across(everything(), 
                        ~ stringr::str_replace_all(., "'", "''"))) |>
        dplyr::mutate_if(check_numeric_types, ~stringr::str_c("'", .x, "'")) |>
        dplyr::mutate(
          dplyr::across(
            dplyr::everything(),
            ~ tidyr::replace_na(.,'NULL')
          )
        )
    }    
    
    if(deletePreviousRecords){
      if(length(deleteMatchingRecords) == 0){
        query <- stringi::stri_c("DELETE FROM ", schema, ".", table, ";")
        logger$debug(query)
        RODBC::sqlQuery(con, query, as.is = T)
      } else {
        whereClause <- build_where_clause(filters = deleteMatchingRecords)
        query <- stringi::stri_c("DELETE FROM ", 
                                 schema, ".", table, " ", whereClause,  ";")
        logger$debug(query)
        RODBC::sqlQuery(con, query, as.is = T)
      }
    }
    
    # Calculate the number of steps required to upload all data.
    nSteps <- ceiling(nrow(data) / 1000)
    if(nSteps > 100){
      logger$warning("The requested operation will perform over 100 steps. This may take some time. Halting the process may result in loss of data.")
    }
    for(i in 1:nSteps){
      maxRowNr <- (i * 1000) - 1
      if(maxRowNr > nrow(data)){
        maxRowNr <- nrow(data)
      }
      slicedData <- data[((i-1) * 1000):maxRowNr, ]
      
      query <- glue::glue("{ifelse(grepl('SQL Server', get_database_info(connInfo)$sqlDriver),'--set nocount on', '')}
                INSERT INTO {schema}.{table} ({paste0('\"', colnames(slicedData), '\"', collapse = ', ')}) 
                VALUES 
                {
                    paste0(
                        map(split(slicedData, seq(nrow(slicedData))), 
                        ~paste(
                            '\t(', 
                            paste0(.x, collapse = ', '), 
                            ')', 
                            sep = '')),
                        collapse = ',\n'
                    )
                }
                ")
      
      # Additional sanitization steps
      query <- gsub(pattern = "'NA'", 
                    replacement = "NULL", 
                    x = query)
      query <- gsub(pattern = "NNULL", 
                    replacement = "NULL", 
                    x = query)
      
      res <- RODBC::sqlQuery(channel = con, 
                             query = query, 
                             as.is = TRUE)
      
      if(!is.null(res)) if(any(grepl("ERROR", res))){
        stop(paste0(res, collapse = ';'))
      }	else {
        logger$debug(
          glue::glue("{nrow(slicedData)} ({(i-1) * 1000}:{maxRowNr}) records written to table {table} on step [{i}/{nSteps}]."))
      }
    } 
    # Apparently SQLite connections don't need to be closed...
    if(!stringr::str_detect(string = get_database_info(connInfo)$sqlDriver, pattern = 'SQLite')){
      RODBC::odbcClose(con)
      logger$debug("Connection to server closed.")
    }
    
  },
  error = function(cond){
    logger$error("There was an error when querying the database.")
    logger$error(cond$message)
    
    RODBC::odbcClose(con)
    logger$debug("Connection to server closed.")
  })  
}

#' @title Update data on a SQL database
#' 
#' @description Tries to perform an UPDATE operation against a table on a SQL
#'   database. If you are updating a few columns only and not entire rows, you
#'   will need to specify val_col. Otherwise, you can leave it empty.
#'   
#' @param data dataframe. Dataframe containing the data to upload.
#' @param table character. Name of the table to write data to.
#' @param index character. Name of the column used as a row identifier.
#' @param connInfo string. Name of the profile to be used on this connection.
#' @param val_col list. Specific columns to update. Default value is NULL,
#'   in which case all columns are updated. Note that if val_col is set to NULL,
#'   all columns need to be present within data.
#' 
#' @note This function has not been updated in a year. Caution is advised.
update_data <- function(data,
                        table,
                        index,
                        connInfo,
                        val_col = NULL){
  if(index %in% colnames(data)){
    con <- create_sql_connection(connInfo = connInfo)
    withCallingHandlers({
      # Colnames sanitization
      if(stringr::str_detect(string = get_database_info(connInfo)$sqlDriver, pattern = 'PostgreSQL')){
        colnames(data) <- stringi::stri_replace_all_fixed(str = colnames(data), 
                                                          pattern = colnames_sanitizer$symbol, 
                                                          replacement = colnames_sanitizer$replacement, 
                                                          vectorize_all = F)
      }
      # String sanitization (common for all databases)
      data <- ungroup(data) |> 
        mutate_all(funs(stringr::str_replace_all(string = ., pattern = "'", replacement = "''"))) |>
        mutate_if(check_numeric_types, ~stringr::str_c("'", .x, "'")) |>
        replace(is.na(.), 'NULL')
      # This one is designed for the upload of one or many specific columns, and many, many entries
      if(length(val_col) != 0){ 
        for(query in data |> 
            select(c(val_col, index)) |>
            mutate(across(everything(), ~paste0(cur_column(), "='", ., "'"))) |> 
            rowwise() |>
            transmute(value = paste0(lapply(val_col, function(x) get(x)), collapse = ", \n\t"),
                      filter = paste0(lapply(index, function(x) get(x)), collapse=" AND ")) |>
            transmute(query = glue::glue("UPDATE {table} \nSET {value} \nWHERE {filter}")) |>
            pull(query)){
          # Additional sanitization steps
          query <- gsub(pattern = "'NA'", replacement = "NULL", x = query)
          query <- gsub(pattern = "NNULL", replacement = "NULL", x = query)
          
          logger$debug(query)
          # Run the query
          RODBC::sqlQuery(con, query, as.is = TRUE)
        }
        logger$debug(glue::glue("Successfully updated {nrow(data)} records from {table}."))
      } else { # This one is a bit more straightforward but may require to have the whole collection of columns present in the remote table
        sqlUpdate(channel = con,
                  dat = data,
                  tablename = table,
                  index = index)
      }
      RODBC::odbcClose(con)
      logger$debug("Connection to server closed.")
    },
    error = function(cond){
      logger$error("There was an error when updating records on the database.")
      logger$error(cond$message)
      RODBC::odbcClose(con)
      logger$debug("Connection to server closed.")
    })  
  } else {
    stop(glue::glue("Index {paste0(index, collapse=', ')} not found in dataset colnames."))
  }
}

#' @title Get database info
#' 
#' @description Function to fetch a connection profile. Profiles to be added.
#' 
#' @param db character. Name of the profile to load.
#' @param dbPath character. Directory to look for the profile. Can be globally
#'   overridden by setting dbDir in .GlobalEnv.
#'
#' @return connInfo list. Contains:
#' - sqlDriver character. String containing the reference for the SQL 
#'   driver to use.
#' - server character. String containing the URL of the SQL server to use.
#' - database character. String containing the name of the SQL database
#'   to use.
#' - credentials List. List containing two keys (uid, pwd) holding strings
#'   containing the user name and password respectively to use to connect to
#'   the SQL database.
#'
#' @note Profiles should contain the following JSON structure:
#' {
#' 	"sqlDriver_windows": "",
#' 	"sqlDriver_linux": "",
#' 	"server": "",
#' 	"port": "",
#' 	"database": "",
#' 	"credentials": {
#' 		"uid": "",
#' 		"pwd": ""
#' 	}
#' }
get_database_info <- function(db, dbPath = 'db'){
  if(exists(x = 'dbDir', envir = .GlobalEnv)){
    # The default value assumes you're running this from one of the folders
    dbPath <- get(x = 'dbDir', envir = .GlobalEnv)
  }
  profiles <- dir(path = dbPath, pattern = '*\\.PROFILE', full.names = F) |>
    stringr::str_remove('\\.PROFILE')
  if(db %in% profiles){
    profile <- jsonlite::read_json(stringr::str_c(dbPath, '/', db, '.PROFILE'))
    profile$sqlDriver <- dplyr::if_else(Sys.info()['sysname'] == 'Linux',
                                        profile$sqlDriver_linux,
                                        profile$sqlDriver_windows)
  } else {
    stop(glue::glue("Profile {stringr::str_remove(db, '\\\\.PROFILE')} was not found. Please check it exists."))
  }
  return(profile)
}


#' @title Create a new table
#' 
#' @description Function to create a new table based on the data types of a
#'   given dataframe.
#' 
#' @param tableName character. Name of the table to create.
#' @param schema character
#' @param data dataframe. Dataframe to be used to infer the column types.
#' @param addIdPk bool. Add an id_pk primary key column?.
#' 
#' @return character. String containing the CREATE TABLE IF NOT EXISTS SQL
#'   query.
create_new_table = function(tableName, schema = 'dbo', data, flavor){
  #data = as.data.frame(lapply(data, type.convert))
  columns = data |>
    dplyr::mutate_if(is.factor, as.character) |>
    purrr::map(class) |> 
    data.frame() |> 
    (function(.)
      tidyr::pivot_longer(., cols = colnames(.)))() |>
    dplyr::mutate(value = get_sql_type(value, flavor)) |>
    (function(.)
      base::split((.), seq(nrow(.))))()
  id_pk = glue::glue('\t"id_pk" {case_when(flavor == "SQLite" ~ "INTEGER PRIMARY KEY AUTOINCREMENT,",
                                             flavor == "SQLite3" ~ "INTEGER PRIMARY KEY AUTOINCREMENT,",
                                             flavor == "SQL Server" ~ "UNIQUEIDENTIFIER NOT NULL PRIMARY KEY DEFAULT NEWID(),",
                                             flavor == "ODBC Driver 18 for SQL Server" ~ "UNIQUEIDENTIFIER NOT NULL PRIMARY KEY DEFAULT NEWID(),",
                                             flavor == "ODBC Driver 17 for SQL Server" ~ "UNIQUEIDENTIFIER NOT NULL PRIMARY KEY DEFAULT NEWID(),",
                                             flavor == "PostgreSQL" ~ "SERIAL PRIMARY KEY,",
                                             T ~ "INTEGER PRIMARY KEY AUTOINCREMENT")}')
  
  ifnotexists = glue::glue(
    case_when(flavor == 'SQLite' ~ "CREATE TABLE IF NOT EXISTS {schema}.{tableName}",
              flavor == 'SQLite3' ~ "CREATE TABLE IF NOT EXISTS {schema}.{tableName}",
              flavor == 'SQL Server' ~ "IF NOT EXISTS (SELECT * FROM sysobjects WHERE name='{tableName}' and xtype='U')
                  CREATE TABLE  {schema}.{tableName}",
              flavor == 'ODBC Driver 18 for SQL Server' ~ "IF NOT EXISTS (SELECT * FROM sysobjects WHERE name='{tableName}' and xtype='U')
                  CREATE TABLE  {schema}.{tableName}",
              flavor == 'ODBC Driver 17 for SQL Server' ~ "IF NOT EXISTS (SELECT * FROM sysobjects WHERE name='{tableName}' and xtype='U')
                  CREATE TABLE  {schema}.{tableName}",
              flavor == 'PostgreSQL' ~ "CREATE TABLE IF NOT EXISTS dbo.{tableName}",
              T ~ "CREATE TABLE IF NOT EXISTS dbo.{tableName}"))
  
  glue::glue("
         {ifnotexists} 
         (
           {
            paste0(
              map(columns, ~paste('\t\"', .x$name,'\" ', .x$value, sep = '')),
              collapse = ',\n')
           }
  );")
  
}

#' @title Check for numeric types
#' 
#' @description Checks whether a variable is not a numeric type. Needs a bit of
#'   a rework in the logic to make sense.
#' 
#' @param x object. An object to verify it's not numeric nor double.
#' 
#' @return boolean. T or F depending whether the variable is not a numeric or
#' double type, or is.
check_numeric_types <- function(x){
  check <- !is.numeric(x) & !is.double(x)
  return(check)
}

#' @title Sanitizer reference for column names.
#' 
#' @description PostgreSQL is not usually happy when given non alphanumeric or
#'   underscore values as column names. This lookup take cares of that.
#' 
#' @name colnames_sanitizer
colnames_sanitizer <- data.frame(symbol = c('-', '<', '>', '+', '(', ')', '[', ']', '{', '}', '/', '\\'),
                                 replacement = c('.dash.', '.lesser.', '.greater.', '.plus.', '.oparenth.', '.cparenth.', '.osqbr.', '.csqbr.', '.ocrbr.', '.ccrbr.', '.slash.', '.bslash.'))


#' @title Build a SQL query
#' 
#' @description Returns a fully formed query.
#' 
#' @param columns character|list. List of columns to retrieve.
#' @param table character. Name of the table to retrieve data from.
#' @param filters character|list. WHERE clause fully formed, or named list 
#'   containing the criteria.
#' @param sanitizeColumns logical. Whether to sanitize the name of the columns
#' following the PostgreSQL specifications or not.
#' 
#' @return query character. A string containing the fully formed query.
#'    
#' @note Building a SQL query using this method has limitations. Bear in mind
#'   that the only query that can be built this way is a SELECT query. There are
#'   additional measures to prevent any additional operation (DROP, INSERT,
#'   DELETE, UPDATE, etc...) from being injected into the final query.
build_sql_query <- function(columns, table, filters,
                            sanitizeColumns = F){
  
  # Build the SELECT clause
  # Prevent SQL injections
  prevent_sql_injection(values = columns)
  selectClause <- build_select_clause(columns = columns, sanitizeColumns = sanitizeColumns)
  
  # Build the FROM clause
  # Prevent SQL injections
  prevent_sql_injection(values = table)
  fromClause <- stringi::stri_c('FROM ', table)
  
  # Build the WHERE clause
  # Prevent sql injections
  prevent_sql_injection(values = filters)
  whereClause <- build_where_clause(filters = filters)
  
  # Build the query
  query <- stringi::stri_c(selectClause, ' \n', fromClause, ' \n', whereClause)
  
  return(query)
}

#' @title Build a WHERE clause
#' 
#' @description Returns a WHERE clause.
#' 
#' @param filters character|list. WHERE clause fully formed, or named list 
#'   containing the criteria.
#' 
#' @return whereClause character. A string containing the WHERE clause.
build_where_clause <- function(filters){
  # Check whether filters is a list or a single value
  if(length(filters) > 1){
    # Assert that, if it is a list, it is named
    if(!is.null(names(filters))) {
      # Transform the list into the adequate WHERE clause
      whereClause <- mapply(function(x, i){
        if(length(x) > 1){
          value <- paste0(" IN ('", 
                          paste0(x, collapse = "', '"), 
                          "')", 
                          collapse = "")
        } else if(is.character(x)){
          if(stringr::str_detect(str_sub(string = x, start = 1, end = 2), '>|<|IS')){
            operator <- stri_extract_first_regex(x, '<>|<=|>=|<|>|IS NOT|IS')[[1]]
            criteria <- stri_replace_first_regex(x, operator, '') |> 
              trimws()
            criteria <- case_when(
              criteria == 'NULL' ~ criteria,
              T ~ stringi::stri_c("'", criteria, "'")
            )
            value <- paste0(' ', operator, ' ', criteria)
          } else {
            value <- paste0(" = '", x, "'")
          }
        } else {
          value <- paste0(' = ', x)
        }
        paste0('[', i, ']', value)
      }, filters, names(filters)) |>
        paste0(collapse = ' \n\tAND ') |>
        stringr::str_c('WHERE ', .)
    } else {
      stop("Filters are neither a named character list, nor a WHERE clause.")
    }
  } else {
    # Check if it's an already formed WHERE clause
    if(stringr::str_starts(filters, 'WHERE')){
      whereClause <- filters
      # Check if it's empty
    } else if(filters == ""){
      whereClause <- filters
    } else if(!is.null(names(filters))){
      # Transform the list into the adequate WHERE clause
      whereClause <- mapply(function(x, i){
        if(length(x) > 1){
          value <- paste0(" IN ('", 
                          paste0(x, collapse = "', '"), 
                          "')", 
                          collapse = "")
        } else if(is.character(x)){
          if(stringr::str_detect(stringr::str_sub(string = x, start = 1, end = 2), 
                                 '>|<|IS')){
            operator <- stringi::stri_extract_first_regex(x, '<>|<=|>=|<|>|IS NOT|IS')[[1]]
            criteria <- stringi::stri_replace_first_regex(x, operator, '') |> 
              trimws()
            criteria <- case_when(
              criteria == 'NULL' ~ criteria,
              T ~ stringi::stri_c("'", criteria, "'")
            )
            value <- paste0(' ', operator, ' ', criteria)
          } else {
            value <- paste0(" = '", x, "'")
          }
        } else {
          value <- paste0(' = ', x)
        }
        paste0('[', i, ']', value)
      }, filters, names(filters)) |>
        paste0(collapse = ' \n\tAND ') |>
        (function(.)
          stringr::str_c('WHERE ', .))()
    } else {
      stop("Filters are neither a named character list, nor a WHERE clause.")
    }
  }
  
  return(whereClause)
}

#' @title Build a SELECT clause
#' 
#' @description Returns a SELECT clause.
#' 
#' @param columns character|list. List of columns to retrieve.
#' @param sanitizeColumns logical. Whether to sanitize the name of the columns
#' following the PostgreSQL specifications or not.
#' 
#' @return selectClause character. A string containing the SELECT clause.
build_select_clause <- function(columns, sanitizeColumns){
  # Sanitize column names
  if(sanitizeColumns){
    columns <- stringi::stri_replace_all_fixed(str = columns, 
                                               pattern = colnames_sanitizer$symbol, 
                                               replacement = colnames_sanitizer$replacement, 
                                               vectorize_all = F)
  }
  selectClause <- ifelse(length(columns)>1, 
                         paste('[', 
                               paste0(columns, collapse = '], [', sep = ""), 
                               ']', sep = ""), 
                         columns) |>
    (function(.)
      stringi::stri_c('SELECT ', .))()
  
}

#' @title Prevent SQL injections
#'
#' @description SQL injections are the process in which an external string,
#'   usually input by a third party, ends up alterating the normal or expected
#'   behavior of a query. This can result in a potential harmful situation in
#'   which the perpetrator acquires access to private tables or even performs
#'   operations that alter the database or the whole system.
#'
#' @param values list. Values to be checked against.
#'
#' @note prevent_sql_injection() will raise an error in case a SELECT, DROP,
#'   INSERT, DELETE or UPDATE keywords are found.
prevent_sql_injection <- function(values){
  if(any(stringr::str_detect(string = unlist(values),
                             pattern = "SELECT|DROP|INSERT|DELETE|UPDATE"))){
    values <- paste0(values, collapse = '; ')
    stop(stringr::str_c('The following values contain illegal words: ', values))
  }
}

#' @title SQL types depending on flavor
#' 
#' @description While SQL is a generic set of standards and guiding principles,
#'   certain details such as data type identifiers might not be the same across
#'   implementations. This function returns the right one depending on the SQL
#'   driver of the connection.
#' 
#' @name get_sql_type
get_sql_type <- function(RType, flavor){
  map(.x = RType, .f = function(x){
    switch(which(c(x %in% names(RSQLTypesLookup[[flavor]]),
                   !x %in% names(RSQLTypesLookup[[flavor]]))),
           RSQLTypesLookup[[flavor]][[x]],
           RSQLTypesLookup[[flavor]][['__generic']])}) |>
    unlist() 
}

#' @title Lookup to transform R classes to SQL data types.
#' 
#' @description Used within get_sql_type to obtain the right SQL data type name
#'   depending on SQL flavor.
#' 
#' @name RSQLTypesLookup
RSQLTypesLookup = list(
  'ODBC Driver 18 for SQL Server' = list(
    "character" = 'nvarchar(4000)',
    "complex" = 'float',
    "double" = 'float',
    "expression" = 'nvarchar(4000)',
    "integer" = 'bigint',
    "list" = 'nvarchar(4000)',
    "logical" = 'bit',
    "numeric" = 'float',
    "single" = 'nvarchar(4000)',
    "raw" = 'filestream',
    "POSIXct" = 'datetime',
    "Date" = 'date',
    "__generic" = 'nvarchar(4000)'
  ),
  'ODBC Driver 17 for SQL Server' = list(
    "character" = 'nvarchar(4000)',
    "complex" = 'float',
    "double" = 'float',
    "expression" = 'nvarchar(4000)',
    "integer" = 'bigint',
    "list" = 'nvarchar(4000)',
    "logical" = 'bit',
    "numeric" = 'float',
    "single" = 'nvarchar(4000)',
    "raw" = 'filestream',
    "POSIXct" = 'datetime',
    "Date" = 'date',
    "__generic" = 'nvarchar(4000)'
  ),
  'SQL Server' = list(
    "character" = 'nvarchar(4000)',
    "complex" = 'float',
    "double" = 'float',
    "expression" = 'nvarchar(4000)',
    "integer" = 'bigint',
    "list" = 'nvarchar(4000)',
    "logical" = 'bit',
    "numeric" = 'float',
    "single" = 'nvarchar(4000)',
    "raw" = 'filestream',
    "POSIXct" = 'datetime',
    "Date" = 'date',
    "__generic" = 'nvarchar(4000)'
  ),
  'SQLite' = list(
    "character" = 'nvarchar',
    "complex" = 'nvarchar',
    "double" = 'float',
    "expression" = 'nvarchar',
    "integer" = 'bigint',
    "list" = 'nvarchar',
    "logical" = 'INTEGER',
    "numeric" = 'float',
    "single" = 'nvarchar',
    "raw" = 'blob',
    "POSIXct" = 'datetime',
    "Date" = 'date',
    "__generic" = 'nvarchar'
  ),
  'SQLite3' = list(
    "character" = 'nvarchar',
    "complex" = 'nvarchar',
    "double" = 'float',
    "expression" = 'nvarchar',
    "integer" = 'bigint',
    "list" = 'nvarchar',
    "logical" = 'INTEGER',
    "numeric" = 'float',
    "single" = 'nvarchar',
    "raw" = 'blob',
    "POSIXct" = 'datetime',
    "Date" = 'date',
    "__generic" = 'nvarchar'
  ),
  'PostgreSQL' = list(
    "character" = 'nvarchar(4000)',
    "complex" = 'nvarchar(4000)',
    "double" = 'float',
    "expression" = 'nvarchar(4000)',
    "integer" = 'bigint',
    "list" = 'nvarchar(4000)',
    "logical" = 'boolean',
    "numeric" = 'float',
    "single" = 'nvarchar(4000)',
    "raw" = 'blob',
    "POSIXct" = 'timestamp',
    "Date" = 'date',
    "__generic" = 'nvarchar(4000)'
  )
)
