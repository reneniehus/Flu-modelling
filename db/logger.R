#' @title Logger factory
#' 
#' @description Factorize a Logger class tailored to your specific needs.
#' 
#' @param name character. Name of the class. Default value is Logger.
#' @param ... character. Named values defining the object type of additional
#'   metadata labels. Valid values are "character", "logical" and "numeric". 
#' 
#' @usage Logger <- forge_logger()
#'
#' @returns a Logger class object.
#'
#' @note This function yields a class. In order to use the logger object, it needs
#'   to be properly instantiated. Please check the documentation on Logger for more
#'   information.
forge_logger = function(name = 'Logger',
                        ...) {
  setRefClass(
    name,
    fields = c(
      logLevel = "character",
      fileName = "character",
      consoleOutput = "logical",
      fileOutput = "logical",
      fileFormat = "character",
      memoryMonitor = "logical",
      fieldList = "character",
      quoteEscape = "character",
      list(...)
    ),
    methods = list(
      initialize = function(...) {
        callSuper(...)
        valid_formats <- c('human', 'csv')
        .self$logLevel <-
          ifelse(length(as.character(logLevel)) == 0,
                 "OFF",
                 as.character(logLevel))
        .self$fileName <-
          ifelse(length(as.character(fileName)) == 0,
                 "",
                 as.character(fileName))
        .self$consoleOutput <-
          ifelse(length(as.logical(consoleOutput)) == 0,
                 TRUE,
                 as.logical(consoleOutput))
        .self$fileOutput <-
          ifelse(
            length(as.logical(fileOutput)) == 0,
            ifelse(fileName == "",
                   FALSE,
                   TRUE),
            as.logical(fileOutput)
          )
        .self$fileFormat <-
          ifelse(
            length(as.character(fileFormat)) == 0,
            'human',
            as.character(fileFormat)
          )
        .self$memoryMonitor <-
          ifelse(length(as.logical(memoryMonitor)) == 0,
                 FALSE,
                 as.logical(memoryMonitor))
        .self$fieldList = setdiff(
          names(.self)[stringi::stri_detect(str = names(.self),
                                            regex = "^[a-zA-Z]")],
          c(
            'logLevel',
            'fileName',
            'consoleOutput',
            'fileOutput',
            'fileFormat',
            'memoryMonitor',
            'fieldList',
            'quoteEscape',
            'initFields',
            'initialize',
            'get_level',
            'info',
            'warning',
            'error',
            'debug',
            'fatal',
            'output'
          )
        )
        for (field in .self$fieldList) {
          .self[[field]] <- get(field)
        }
        # Health check
        if(.self$fileFormat == 'csv' & !endsWith(.self$fileName, '.csv')){
          base::warning('Adding csv extension to file name.')
          .self$fileName <- stringr::str_c(.self$fileName, '.csv')
        }
        if(.self$fileFormat == 'human' & !endsWith(.self$fileName, '.log')){
          base::warning('Adding log extension to file name.')
          .self$fileName <- stringr::str_c(.self$fileName, '.log')
        }
        if(!.self$fileFormat %in% valid_formats){
          stop(stringr::str_c("File format ", 
                              .self$fileFormat), " is not valid.")
        }
        if(!dir.exists(dirname(.self$fileName))){
          stop(stringr::str_c("Can't access directory ",
                              dirname(.self$fileName), "."))
        }
        if(file.exists(.self$fileName)){
          base::warning(stringr::str_c('File ', .self$fileName, ' exists already.'))
        }
        if(.self$fileFormat == 'csv'){
          header <- c('Timestamp', 'Memory', 
                      .self$fieldList,'LogLevel', 'Message')
          if(file.exists(.self$fileName)){
            fileHeader <- stringr::str_split(
              readLines(.self$fileName, n = 1),
              ',', simplify = T)
            if(ifelse(length(unlist(as.list(fileHeader))) == length(header), F, !all(unlist(as.list(fileHeader)) == header))){
              stop(stringr::str_c("Target file ", 
                                  .self$fileName, 
                                  " contains a different header, appending is not possible. Use a different file name"))
            } else {
              base::warning("Appending new records to CSV file.")
            }
          } else {
            cat(
              header,
              file = .self$fileName,
              sep = ','
            )
            cat('\n', file = .self$fileName, append = T)
          }
        }
      },
      get_level = function() {
        levels = c("OFF", "FATAL", "ERROR",
                   "WARN", "INFO", "DEBUG")
        return(which(levels %in% .self$logLevel))
      },
      output = function(msg) {
        if (consoleOutput) {
          cat(
            stringr::str_remove_all(
              paste0(msg[1:(length(msg)-2)], collapse = ']['),
              '\\[\\]'),
            '] ',
            msg[length(msg)-1],
            msg[length(msg)],
            sep = ""
          )
        }
        if (fileOutput) {
          msg <- msg[3:length(msg)-1]
          if(.self$fileFormat == 'human'){
            cat(stringr::str_c('[', msg[1:(length(msg)-1)], ']'),
                ' ',
                msg[length(msg)],
                '\n',
                sep = "",
                file = .self$fileName,
                append = T
            )
          } else if(fileFormat == 'csv'){
            cat(
              stringr::str_c('"', stringr::str_replace_all(msg, '"', '""'), '"'),
              file = .self$fileName,
              append = T,
              sep = ','
            )
            cat('\n', file = .self$fileName, append = T)
          }
        }
      },
      debug = function(msg) {
        if (get_level() >= 6) {
          msg <-
            c(
              "\u001b[95m\u001b[1m\u001b[1m[",
              as.character(Sys.time()),
              if (memoryMonitor)
                stringr::str_c(
                  ceiling(sum(.Internal(gc(
                    F, F, F
                  ))[c(3, 4)]) * 1024 * 1024),
                  ' bytes')
              else
                "",
              unlist(
                lapply(.self$fieldList, function(x)
                  get(x))),
              "DEBUG",
              msg,
              "\u001b[39m\u001b[22m\u001b[49m\n"
            )
          output(msg)
        }
      },
      info = function(msg) {
        if (get_level() >= 5) {
          msg <-
            c(
              "\u001b[92m\u001b[1m\u001b[1m[",
              as.character(Sys.time()),
              if (memoryMonitor)
                stringr::str_c(
                  ceiling(sum(.Internal(gc(
                    F, F, F
                  ))[c(3, 4)]) * 1024 * 1024),
                  ' bytes')
              else
                "",
              unlist(
                lapply(.self$fieldList, function(x)
                  get(x))),
              "INFO",
              msg,
              "\u001b[39m\u001b[22m\u001b[49m\n"
            )
          output(msg)
        }
      },
      warning = function(msg) {
        if (get_level() >= 4) {
          msg <-
            c(
              "\u001b[93m\u001b[1m\u001b[1m[",
              as.character(Sys.time()),
              if (memoryMonitor)
                stringr::str_c(
                  ceiling(sum(.Internal(gc(
                    F, F, F
                  ))[c(3, 4)]) * 1024 * 1024),
                  ' bytes')
              else
                "",
              unlist(
                lapply(.self$fieldList, function(x)
                  get(x))),
              "WARNING",
              msg,
              "\u001b[39m\u001b[22m\u001b[49m\n"
            )
          output(msg)
        }
      },
      error = function(msg) {
        if (get_level() >= 3) {
          msg <-
            c(
              "\u001b[91m\u001b[1m\u001b[1m[",
              as.character(Sys.time()),
              if (memoryMonitor)
                stringr::str_c(
                  ceiling(sum(.Internal(gc(
                    F, F, F
                  ))[c(3, 4)]) * 1024 * 1024),
                  ' bytes')
              else
                "",
              unlist(
                lapply(.self$fieldList, function(x)
                  get(x))),
              "ERROR",
              msg,
              "\u001b[39m\u001b[22m\u001b[49m\n"
            )
          output(msg)
        }
      },
      fatal = function(msg) {
        if (get_level() >= 2) {
          msg <-
            c(
              "\u001b[31m\u001b[1m\u001b[40m[",
              as.character(Sys.time()),
              if (memoryMonitor)
                stringr::str_c(
                  ceiling(sum(.Internal(gc(
                    F, F, F
                  ))[c(3, 4)]) * 1024 * 1024),
                  ' bytes')
              else
                "",
              unlist(
                lapply(.self$fieldList, function(x)
                  get(x))),
              "FATAL",
              msg,
              "\u001b[39m\u001b[22m\u001b[49m\n"
            )
          output(msg)
        }
      }
    )
  )
}

#' @title Logger
#' 
#' @description Simple Logger class to generate logger objects. Logger lets you keep
#'   track of different logging messages on up to six different log levels (debug, info,
#'   warn, error, fatal and off).
#'   
#' @param logLevel character. Log level on which the logger object acts. Valid values are
#    DEBUG, INFO, WARN, ERROR, FATAL and OFF (note the upper case).
#' @param fileName character. Path and name of the file to store the logfiles. Note that
#'   the directory needs to exist.
#' @param consoleOutput logical. Whether to output messages to the console or not.
#' @param fileOutput logical. Whether to output messages to a file or not.
#' @param memoryMonitor logical. Whether to include the amount of memory consumed by R or not.
#'   Note that this will impact performance.
#' @param ... implemented by forge. Labels to add context or information to the log messages.
#' 
#' @usage logger <- Logger(logLevel = 'INFO', consoleOutput = T, memoryMonitor = T)
#'
#' @returns a logger object.
#'
#' @note To push a logging message you can access the specific methods using the $ operator
#'   (i.e. logger$info("This is an info-level message")). Note that to modify the logging level 
#'   you can reassign the value of logLevel (i.e. logger$logLevel <- "OFF").
#'   As acknowledged on the README, the logger messages can behave oddly in Windows. Additionally,
#'   you may want to set the logLevel to OFF before rendering a markdown file, as they may get
#'   redirected to the resulting .MD file.
#' @name Logger
NULL