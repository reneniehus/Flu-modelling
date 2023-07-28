#############################
# v. October 2022
# Modified by Adrian Prodan
# Absit omen
#############################

# 1) Config & Initialisation --------------------------------------------------------

list.of.packages <- c("RODBC", "mem")
missing.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(missing.packages)) {
  install.packages(missing.packages)
}

# update.packages(ask=FALSE) #update all packages that have new versions available, do not *ask* for confirmation

library(RODBC)
library(mem)

outputPath <- here("code/flu/MEM model/MEMOutput/")
setwd(outputPath) # if you get this error, the folder does not exist yet and needs to be created: Error in setwd(outputPath) : cannot change working directory

if (length(dir()) != 0) {
  print(
    paste(
      "The folder (",
      outputPath,
      ") needs to be empty in the beginning, because of possible appending problems."
    )
  )
  stop(outputPath, " folder not empty")
}

connectionStringSource <- "Driver=SQL Server;Server=ZSQLCL2.ecdcdmz.europa.eu,4964;Database=TessyDM;Trusted_Connection=TRUE;"
connectionStringTarget <- "Driver=SQL Server;Server=ZSQLCL2.ecdcdmz.europa.eu,4964;Database=TessyDW;Trusted_Connection=TRUE;"

NEWconnectionStringSource <- "Driver=SQL Server;Server=NSQL3;Database=INFL;Trusted_Connection=TRUE;"

channelSource = odbcDriverConnect(connectionStringSource)
channelTarget <- odbcDriverConnect(connectionStringTarget)
NEWchannelTarget <- odbcDriverConnect(NEWconnectionStringSource)

filenameBaseline <- 'Baseline.csv'
filenameTypicalCurve <- 'TypicalCurve.csv'

#country-specific baseline and curve files # Adrian 2022: not sure what these are doing here; the second one is unused..?
filenameBaselineXX <- 'Baseline_XX.csv'
filenameTypicalCurveXX <- 'TypicalCurve_XX.csv'

#countriesWithSurveillanceType <- sqlQuery(channelSource, "exec spMEMGetCountryList 'FI'")
countriesWithSurveillanceType <- sqlQuery(channelSource, "exec spMEMGetCountryList")

countriesILI <- countriesWithSurveillanceType[1][countriesWithSurveillanceType[2] == 'ILI']
countriesARI <- countriesWithSurveillanceType[1][countriesWithSurveillanceType[2] == 'ARI']

# unclear why we try to remove NAs, doesn't look like there are any..?
countriesILI <- countriesILI[!(is.na(countriesILI))]
countriesARI <- countriesARI[!(is.na(countriesARI))]

results <- numeric()
results_by_country <- numeric()

get.number.of.zeroes.and.missings <- function(x){
  sum(!(is.na(x) | x == 0))
}






# 2. Typical Curve Data Frame ---------------------------------------------

# Creating the Typical Curve data frame, including data types and column names;
typicalCurve <- data.frame(ReportingCountryMnemonic = character(),
                           NominatorType = character(), 
                           AbsoluteWeek = integer(),
                           RelativeWeek = integer(),
                           IntensityLower = integer(),
                           Intensity = integer(),
                           IntensityUpper = integer()
)

#writing the column names as fake values in the file, so that we can append later without column names; otherwise colnames would be repeated in the file at every append;
write.table( 
  typicalCurve,
  file = filenameTypicalCurve,
  dec = ".",
  append = FALSE,
  quote = FALSE, # FALSE = do not surround char or factor column values in double quotes (")
  qmethod = "escape", # escape existing double quotes (") in the data with a backslash ( a " becomes \" )
  sep = ",",
  row.names = F,
  col.names = T
) 

# proof that write.table() is better than write.csv() which adds an extra comma in the beginning of the row, probably because the first col is considered to be the row nr
# write.csv.TC.filename <- paste('write.csv.', filenameTypicalCurve, sep = )
# 
# write.csv( #creating the csv file and writing a header (col names) into it; will be appended later;
#   x = typicalCurve,
#   file = write.csv.TC.filename,
#   quote = FALSE
# ) 

 run.model <- function(countries, source) {

# TEMPORARY args for debugging, make sure these are commented if not debugging:
# countries <- countriesILI;
# source <- 'ILI'

  for (i in 1:length(countries)) {
    country <- countries[i]

    sqlCode <-
      paste(
        "EXEC [dbo].[spMEMGetInflClin] @Region = N'",
        country,
        "', @Numerator = N'",
        source,
        "'",
        sep = ''
      )
    
    inputData <- sqlQuery(channelSource, sqlCode)
    
    
    if (is.null(dim(inputData)[2])) { # if no data found for any season (I think...?)
      results <- rbind(results,
                       c(country,
                         source,
                         rep(NA, 10))) # add fake NAs instead of data
    } else {
      i.data = as.data.frame(inputData[, 2:(dim(inputData)[2])]) # remove the first column, which seems to be the week number;
      rownames(i.data) <- as.character(inputData[, 1]) # add the week numbers as row names
      
      # Keep only seasons where we have data for more than 10 weeks; debug code below to see where this applies;
      ind <- apply(i.data, 2, get.number.of.zeroes.and.missings) > 10
      old.i.data <- i.data
      
      i.data <- i.data[ind]
      
      if (!identical(old.i.data, i.data)) {
        print('Debug, pls ignore: filtering based on zeros and missings finally did something!')
        print('Old i.data:')
        print(old.i.data)
        print('New i.data:')
        print(i.data)
        print('End of debug output.')
      }
      
      if (dim(i.data)[2] < 2) { # if fewer than two seasons of data available (I think...)
        results <- rbind(results,
                         c(country,
                           source,
                           rep(NA, 3),
                           dim(i.data)[2],
                           rep(NA, 4))) # then fill the report with NAs
      } else { 
        epi <- memmodel(i.data) # run the model
        # thresholds
        # epi$pre.post.intervals[,3]
        # mean start
        # epi$mean.start
        
        
        # Create and write as a standalone CSV a data frame with the data for the current Country and Source ----
        epi.typical.curve.df <- data.frame(matrix(ncol = 2, nrow = nrow(epi$typ.curve)))
        colnames(epi.typical.curve.df) <- c('RelativeWeek', 'Intensity')
        epi.typical.curve.df[, 1] = rownames(epi$typ.curve)
        epi.typical.curve.df[, 2] = epi$typ.curve[, 2]
        
        write.table( 
          epi.typical.curve.df,
          file = paste0(source, "_", country, ".csv"),
          dec = ".",
          append = FALSE,
          quote = FALSE,
          sep = ",",
          row.names = F,
          col.names = T
        )
        
        # Write out to two CSV files the current country's Typical Curve ----
        
        typicalCurve <- data.frame(
          ReportingCountryMnemonic = country,
          NominatorType            = source,
          AbsoluteWeek             = c(
            paste('W', c(40:52), sep = ''),
            paste('W0', c(1:9), sep = ''),
            paste('W', c(10:(
              length(epi$typ.curve[, 2]) - 13
            )), sep = '')
          ),
          RelativeWeek             = c(1:length(epi$typ.curve[, 2])),
          IntensityLower           = epi$typ.curve[, 1],
          Intensity                = epi$typ.curve[, 2],
          IntensityUpper           = epi$typ.curve[, 3]
        )
        
        write.table( # append the typical curve of the current country to the common Typical Curve CSV file;
          typicalCurve,
          file = filenameTypicalCurve,
          dec = ".",
          append = TRUE,
          quote = FALSE,
          sep = ",",
          row.names = F,
          col.names = F
        )
        
        filenameTypicalCurve_XX <- paste0("TypicalCurve_", country, ".csv") # the name for the country-specific csv file with the typical curve
        write.table( # generate the stand-alone CSV for the current country;
          typicalCurve,
          file = filenameTypicalCurve_XX,
          dec = ".",
          append = FALSE,
          quote = FALSE,
          sep = ",",
          row.names = FALSE,
          col.names = TRUE
        )
        
        #if (saveToDatabase) {
        #	sqlSave(channelTarget, typicalCurve, tablename = "mart.mINFLMEMTypicalCurve_Adrian", rownames = FALSE, append = TRUE)
        #}
        
        # plot the model to a country-specific PNG file  ----
        png(
          file = paste0(source, "_", country, ".png"),
          width = 1280,
          height = 768,
          units = "px",
          pointsize = 12
        )
        plot(epi)
        dev.off()
        
        
        results <- rbind(
          results,
          c(
            country,
            source,
            round(epi$pre.post.intervals[, 3], 2),
            # pre and post-epidemic thresholds
            epi$mean.start[1],
            # avg start of epidemic ( maybe counted from season start..? (Adrian))
            epi$mean.start[2],
            # avg start of epidemic ( maybe counted as actual week nr..? (Adrian))
            epi$mean.length,
            # avg length of epidemic
            dim(i.data)[2],
            # number of seasons
            epi$epi.intervals[1, 4],
            # 40% CI for maximum of epidemic curve
            epi$epi.intervals[2, 4],
            # 90% CI for maximum of epidemic curve
            epi$epi.intervals[3, 4],
            # 97.5% CI for maximum of epidemic curve
            0                                    # IsUsageApproved
          )
        )
        
        
        results_by_country_colnames <- c(
          "ReportingCountryMnemonic",
          "NominatorType",
          "PreEpidemicThreshold",
          "PostEpidemicThreshold",
          "MeanStartWeek(ofSeason)",
          "MeanStartWeek(weekNr)",
          "MeanLength",
          "NumberOfSeasons",
          "CI40",
          "CI90",
          "CI97.5",
          "IsUsageApproved"
        )
        
        results_by_country <-
          rbind(
            results_by_country_colnames,
            c(
              country,
              source,
              round(epi$pre.post.intervals[, 3], 2),
              # pre and post-epidemic thresholds
              epi$mean.start[1],
              # avg start of epidemic ( maybe counted from season start..? (Adrian))
              epi$mean.start[2],
              # avg start of epidemic ( maybe counted as actual week nr..? (Adrian))
              epi$mean.length,
              # avg length of epidemic
              dim(i.data)[2],
              # number of seasons
              epi$epi.intervals[1, 4],
              # 40% CI for maximum of epidemic curve
              epi$epi.intervals[2, 4],
              # 90% CI for maximum of epidemic curve
              epi$epi.intervals[3, 4],
              # 97.5% CI for maximum of epidemic curve
              0                                    # IsUsageApproved
            )
            
          )
      }
    }
    
    # write out the country-specific CSV file with the baseline ----
    
    filenameBaselineXX <- paste0("Baseline_", country, ".csv") #generating the country-speific csv file with the baseline
    
    #print(c(country, "si in plus", results_by_country[2,1]))
    if (country == results_by_country[2, 1])
      #bug fix for when the data for a previous country was added to the baseline file of the current one # Adrian 2022: I no longer remember why this is here
    {
      write.table(
        results_by_country,
        file = filenameBaselineXX,
        dec = ".",
        append = TRUE,
        quote = FALSE,
        sep = ",",
        row.names = F,
        col.names = F
      )
    }
    
  }
  
 return(results)
}

ari <- run.model(countriesARI, 'ARI')
ili <- run.model(countriesILI, 'ILI')

results <- as.data.frame(rbind(cbind(ili),
                               cbind(ari)))

colnames(results) <- c(
  "ReportingCountryMnemonic",
  "NominatorType",
  "PreEpidemicThreshold",
  "PostEpidemicThreshold",
  "MeanStartWeek(ofSeason)",
  "MeanStartWeek(weekNr)",
  "MeanLength",
  "NumberOfSeasons",
  "CI40",
  "CI90",
  "CI97.5",
  "IsUsageApproved"
)

write.table(
  results,
  file = filenameBaseline,
  dec = ".",
  append = FALSE,
  quote = FALSE,
  sep = ",",
  row.names = F,
  col.names = T
)
# write.table( # unclear why this was running twice...
#   results,
#   file = filenameBaseline,
#   dec = ".",
#   append = FALSE,
#   quote = FALSE,
#   sep = ",",
#   row.names = F,
#   col.names = T
# )

odbcCloseAll()
