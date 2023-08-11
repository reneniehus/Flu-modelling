#Reshape INFLCLIN and INFLWAGGR data
#Access to ZSQLCL5.ECDCDMZ.EUROPA.EU is required



reshape_INFL_data<-function(SQLite_files = T,
                            save_files = F,
                            destination_folder = getwd())
{
  library(RODBC)
  library(dplyr)
  library(DBI)
  library(glue)
  library(fst)
  library(tidyverse)
  
  
  
  #Set-up connection to ECDC DW (ZSQLCL5.ECDCDMZ.EUROPA.EU)
  connStr <-"Driver=SQL Server;Server=ZSQLCL5.ECDCDMZ.EUROPA.EU;Database=INFL;Trusted_Connection=TRUE;"
  mycon <- odbcDriverConnect(connStr)
  
  
  
  #Reshape Denominator data within fClinicalWeeklyOverview data into long-format
  ClinicalWeeklyDenominators <- sqlQuery(mycon,"SELECT RecordId, SubjectCode,LocationCode, ReportingCountry, SeasonCode,
                   TimeCode, NumberOfILICases, NumberOfARICases, TotalNumberOfReportingPhysicians AS NumberOfReportingPhysicians
                          FROM [INFL].[infl].[fClinicalWeeklyOverview]
                          WHERE ValidTo is NULL",
                   stringsAsFactors = FALSE,
                   na.strings="") %>%
    mutate(SurvType="STL") %>%
    pivot_longer(cols = starts_with("NumberOf"),
                 names_to = "Indicator",
                 names_prefix = "NumberOf",
                 values_to = "Value",
                 values_drop_na = TRUE)
  
  
  
  #Reshape Qualitative data within fClinicalWeeklyOverview data into long-format
  WeeklyQualitative <- sqlQuery(mycon,"SELECT RecordId, SubjectCode,LocationCode, ReportingCountry, SeasonCode,
                   TimeCode, GeographicSpreadCode AS Ind_GeographicSpread, IntensityCode AS Ind_Intensity,
                   ImpactCode AS Ind_Impact, TrendCode AS Ind_Trend
                          FROM [INFL].[infl].[fClinicalWeeklyOverview]
                          WHERE ValidTo is NULL",
                   stringsAsFactors = FALSE,
                   na.strings="") %>%
    mutate(SurvType="STL") %>%
    pivot_longer(cols = starts_with("Ind_"),
                 names_to = "Indicator",
                 names_prefix = "Ind_",
                 values_to = "Value",
                 values_drop_na = TRUE)
  
  
  
  
  #Prepare fClinicalWeeklyByAge data
  ClinicalWeeklyByAge <- sqlQuery(mycon,"SELECT RecordId, SubjectCode,LocationCode, ReportingCountry, SeasonCode,
                   TimeCode, AgeClassificationCode AS Age, ClinicalSyndromeCode AS Indicator, NumberOfCases,
                   PopulationDemographyAsReported AS Denominator
                          FROM [INFL].[infl].[fClinicalWeeklyByAge]
                          WHERE ValidTo is NULL",
                   stringsAsFactors = FALSE,
                   na.strings="") %>% mutate(SurvType="STL") %>%
    mutate(Age=case_when(Age=="Age05_14" ~ "05-14",
                         Age=="Age00_04" ~ "00-04",
                         Age=="Age15_64" ~ "15-64",
                         Age=="Age65+" ~ "65+",
                         Age=="AgeUnk" ~ "Unk"))
  
  
  
  #Prepare VirologicalDetectionsWeekly data combining denominator and numerature data
  VirologicalDetectionsWeeklyDenominators <- sqlQuery(mycon,"SELECT RecordId, SubjectCode,LocationCode, ReportingCountry, SeasonCode,
         TimeCode,SurveillanceSystemTypeCode AS SurvType,DiseaseCode, TotalNumberOfSpecimensAnalysed AS Denominator
                          FROM [INFL].[infl].[fVirologicalWeeklyOverview]
                          WHERE ValidTo is NULL",
         stringsAsFactors = FALSE,
         na.strings="")
  
  
  
  
  VirologicalDetectionsWeekly <- sqlQuery(mycon,"SELECT RecordId, SubjectCode,LocationCode, ReportingCountry, SeasonCode,
         TimeCode,SurveillanceSystemTypeCode AS SurvType,DiseaseCode, PathogenTypeCode, PathogenSubtypeCode,
         NumberOfSpecimensDetected AS Value
                          FROM [INFL].[infl].[fVirologicalWeeklyByVirusType]
                          WHERE ValidTo is NULL",
         stringsAsFactors = FALSE,
         na.strings="") %>%
    left_join(.,VirologicalDetectionsWeeklyDenominators, by=c("RecordId", "SubjectCode","LocationCode", "ReportingCountry", "SeasonCode",
                                                              "TimeCode","SurvType","DiseaseCode")) %>%
    mutate(SubjectCode="INFLVIR")
  
  
  
  # MEM Theshold data
  Thresholds <- sqlQuery(mycon,"SELECT LocationCode, SeasonCode, ThresholdCode, ClinicalSyndromeCode, ValueThreshold FROM [INFL].[common].[fThresholds]
                              WHERE ValidTo is NULL AND IsUsageApproved = 1",
                         stringsAsFactors = FALSE,
                         na.strings="")
  
  
  
  
  
  if(SQLite_files){
    con <- dbConnect(RSQLite::SQLite(), ":memory:")
    dbListTables(con)
    dbWriteTable(con, "ClinicalWeeklyDenominators", ClinicalWeeklyDenominators)
    dbWriteTable(con, "WeeklyQualitative", WeeklyQualitative)
    dbWriteTable(con, "ClinicalWeeklyByAge", ClinicalWeeklyByAge)
    dbWriteTable(con, "VirologicalDetectionsWeekly", VirologicalDetectionsWeekly)
    dbWriteTable(con, "MEMThresholds", Thresholds)
    return(con)
  }
  
  
  
  if(!SQLite_files){
    return(list(ClinicalWeeklyDenominators=ClinicalWeeklyDenominators,
                WeeklyQualitative=WeeklyQualitative,
                ClinicalWeeklyByAge=ClinicalWeeklyByAge,
                VirologicalDetectionsWeekly=VirologicalDetectionsWeekly,
                Thresholds=Thresholds))
  }
  
  
  
  if(save_files){
    destination_folder_final = glue('{destination_folder}/data/')
    dir.create(destination_folder_final)
    data_fst(ClinicalWeeklyDenominators, glue('{destination_folder_final}/ClinicalWeeklyDenominators.fst'))
    data_fst(WeeklyQualitative, glue('{destination_folder_final}/WeeklyQualitative.fst'))
    data_fst(ClinicalWeeklyByAge, glue('{destination_folder_final}/ClinicalWeeklyByAge.fst'))
    data_fst(Thresholds, glue('{destination_folder_final}/Thresholds.fst'))
    
    
    
  }
}


X <- reshape_INFL_data(SQLite_files = F)
X$ClinicalWeeklyByAge -> Y
Y$date = as.Date(paste(sub("W","",X$ClinicalWeeklyByAge$TimeCode), 1, sep="-"), "%Y-%U-%u")

X$ClinicalWeeklyByAge %>% 
  group_by(ReportingCountry, date, Indicator, SeasonCode) %>%
  summarise(value = sum(NumberOfCases)) %>%
  ungroup() %>%
  filter(ReportingCountry=="AT",Indicator=="ILI", date>"2009-01-01") %>%  #View()
  ggplot(aes(x=date, y=value)) + 
  geom_line() +
  geom_hline(yintercept = c(250,500,750))
  #facet_wrap(~SeasonCode)


X$Thresholds %>% #View()
  filter(LocationCode=="AT", ClinicalSyndromeCode=="ILI")

         