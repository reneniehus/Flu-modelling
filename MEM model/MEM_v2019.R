####################
# V. 07.2018
####################
list.of.packages <- c("RODBC", "mem")
missing.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(missing.packages))
  install.packages(missing.packages)

#update.packages(ask=FALSE) #update all packages that have new versions available, do not *ask* for confirmation
old.packages()



library(RODBC)
library(mem)

connectionStringSource <-
  "Driver=SQL Server;Server=ZSQLCL2.ecdcdmz.europa.eu,4964;Database=TessyDM;Trusted_Connection=TRUE;"
connectionStringTarget <-
  "Driver=SQL Server;Server=ZSQLCL2.ecdcdmz.europa.eu,4964;Database=TessyDW;Trusted_Connection=TRUE;"

outputPath <- "./Output/"
# outputPath <- "C:/MEMOutput/"
setwd(outputPath)
if (length(dir()) != 0)
{
  print(
    paste(
      "The folder (",
      outputPath,
      ") needs to be empty in the beginning, because of possible appending problems."
    )
  )
  stop(outputPath, " folder not empty")
}

excelFilenameBaseline <- paste(outputPath, "Baseline.csv", sep = "")
excelFilenameTypicalCurve <-
  paste(outputPath, "TypicalCurve.csv", sep = "")

#country-specific baseline and curve files
excelFilenameBaselineXX <-
  paste(outputPath, "Baseline_XX.csv", sep = "")
excelFilenameTypicalCurveXX <-
  paste(outputPath, "TypicalCurve_XX.csv", sep = "")

# channelSource <- odbcConnect(dataSource)
channelSource = odbcDriverConnect(connectionStringSource)
channelTarget <- odbcDriverConnect(connectionStringTarget)


#countriesWithSurveillanceType <- sqlQuery(channelSource, "exec spMEMGetCountryList 'FI'")
countriesWithSurveillanceType <-
  sqlQuery(channelSource, "exec spMEMGetCountryList")

countriesILI <-
  countriesWithSurveillanceType[1][countriesWithSurveillanceType[2] == 'ILI']
countriesARI <-
  countriesWithSurveillanceType[1][countriesWithSurveillanceType[2] == 'ARI']

countriesILI <- countriesILI[!(is.na(countriesILI))]
countriesARI <- countriesARI[!(is.na(countriesARI))]

results <- numeric()
results_by_country <- numeric()

number.of.zeroes.and.missings <- function(x)
  sum(!(is.na(x) | x == 0))


#***********************************************************************************************
{
  #faking the column names as actual values and starting a new file with just these column names, for the merged curves file
  typicalCurve <- data.frame(
    ReportingCountryMnemonic = 'ReportingCountryMnemonic',
    NominatorType            = 'NominatorType',
    AbsoluteWeek             = 'AbsoluteWeek',
    RelativeWeek             = 'RelativeWeek',
    IntensityLower           = 'IntensityLower',
    Intensity                = 'Intensity',
    IntensityUpper           = 'IntensityUpper'
  )
  
  write.table(
    typicalCurve,
    file = excelFilenameTypicalCurve,
    dec = ".",
    append = FALSE,
    quote = FALSE,
    sep = ",",
    row.names = F,
    col.names = F
  ) #writing the fake values in the file
  }
#***********************************************************************************************


proceed <- function(countries, source) {
  for (i in 1:length(countries)) {
    country <- countries[i]
    #sql <- paste("EXEC spMEMGetInflClin '", country, "', '", source, "'", sep = '') #*******
    sql <-
      paste(
        "EXEC [dbo].[spMEMGetInflClin] @Region = N'",
        country,
        "', @Numerator = N'",
        source,
        "'",
        sep = ''
      ) #*******
    
    inputData <- sqlQuery(channelSource, sql)
    
    
    if (is.null(dim(inputData)[2])) {
      results <- rbind(results,
                       c(country,
                         source,
                         rep(NA, 10)))
    } else {
      i.data = as.data.frame(inputData[, 2:(dim(inputData)[2])])
      rownames(i.data) <- as.character(inputData[, 1])
      
      ind <- apply(i.data, 2, number.of.zeroes.and.missings) > 10
      i.data <- i.data[ind]
      
      if (dim(i.data)[2] < 2) {
        results <- rbind(results,
                         c(country,
                           source,
                           rep(NA, 3),
                           dim(i.data)[2],
                           rep(NA, 4)))
      } else {
        epi <- memmodel(i.data)
        # thresholds
        # epi$pre.post.intervals[,3]
        # mean start
        # epi$mean.start
        
        {
          #faking the column names for the ILI_XX.csv and ARI_XX.csv files
          epi_typ_curve_ColumnNames <- data.frame(RelativeWeek 	= 'RelativeWeek',
                                                  Intensity 		= 'Intensity')
          write.table(
            epi_typ_curve_ColumnNames,
            file = paste(outputPath, source, "_", country, ".csv", sep = ""),
            dec = ".",
            append = FALSE,
            quote = FALSE,
            sep = ",",
            row.names = F,
            col.names = F
          )
        }
        
        write.table(
          epi$typ.curve[, 2],
          file = paste(outputPath, source, "_", country, ".csv", sep = ""),
          dec = ".",
          append = TRUE,
          quote = FALSE,
          sep = ",",
          row.names = T,
          col.names = F
        )
        
        #***********************************************************************************************
        {
          #faking the column names as actual values and starting a new file with just these column names, for the country-specific curve
          
          excelFilenameTypicalCurve_XX <-
            paste(outputPath, "TypicalCurve_", country, ".csv", sep = "") #generating the country-speific csv file with the typical curve
          
          typicalCurve <- data.frame(
            ReportingCountryMnemonic = 'ReportingCountryMnemonic',
            NominatorType            = 'NominatorType',
            AbsoluteWeek             = 'AbsoluteWeek',
            RelativeWeek             = 'RelativeWeek',
            IntensityLower           = 'IntensityLower',
            Intensity                = 'Intensity',
            IntensityUpper           = 'IntensityUpper'
          )
          
          write.table(
            typicalCurve,
            file = excelFilenameTypicalCurve_XX,
            dec = ".",
            append = TRUE,
            quote = FALSE,
            sep = ",",
            row.names = F,
            col.names = F
          ) #printing the fake values
        }
        #***********************************************************************************************
        
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
        
        write.table(
          typicalCurve,
          file = excelFilenameTypicalCurve,
          dec = ".",
          append = TRUE,
          quote = FALSE,
          sep = ",",
          row.names = F,
          col.names = F
        )
        
        write.table(
          typicalCurve,
          file = excelFilenameTypicalCurve_XX,
          dec = ".",
          append = TRUE,
          quote = FALSE,
          sep = ",",
          row.names = F,
          col.names = F
        )
        
        
        #if (saveToDatabase) {
        #	sqlSave(channelTarget, typicalCurve, tablename = "mart.mINFLMEMTypicalCurve_Adrian", rownames = FALSE, append = TRUE)
        #}
        
        #bmp(file=paste(outputPath, source, "_", country, ".bmp", sep=""))
        png(
          file = paste(outputPath, source, "_", country, ".png", sep = ""),
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
    
    excelFilenameBaselineXX <-
      paste(outputPath, "Baseline_", country, ".csv", sep = "") #generating the country-speific csv file with the baseline
    
    #print(c(country, "si in plus", results_by_country[2,1]))
    if (country == results_by_country[2, 1])
      #bug fix for when the data for a previous country was added to the baseline file of the current one
    {
      write.table(
        results_by_country,
        file = excelFilenameBaselineXX,
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

ari <- proceed(countriesARI, 'ARI')
ili <- proceed(countriesILI, 'ILI')

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
  file = excelFilenameBaseline,
  dec = ".",
  append = FALSE,
  quote = FALSE,
  sep = ",",
  row.names = F,
  col.names = T
)
write.table(
  results,
  file = excelFilenameBaseline,
  dec = ".",
  append = FALSE,
  quote = FALSE,
  sep = ",",
  row.names = F,
  col.names = T
)

odbcCloseAll()
