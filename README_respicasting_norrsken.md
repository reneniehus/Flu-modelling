Document describing in very brief how to run ECDC's norrsken model

As desired in the Flu-forecasting repo, everything should be run from the 
00_main.R script in the code folder. 

The user can disregard most lines of that script and move all the way down to the 
# ---- |-Respicasting: Norrsken ---- section. 

That sections sources scripts are can be run independently, meaning they will
each load all required libraries and create the files that they needed. 
Internet connection is needed to load the data. 

Example: run the line that does source("./code/respicasting_ILI.R")
this script will 
-load the libraries
-define settings and parameters (fit_rerun = T should be the default to refit
the stan model, while fit_rerun = F loads existing fits if run before on the 
local machine and it much faster)
-loop through green (log-transform) and blue (sqrt-transform) model
-loading the github data for the respective indicator
-filtering out countries
-looping through countries
-preparing country-specific dataframe for stan model
-(optional plotting)
-fit the model
-(optional output plotting)
-wrangle the output to be in submission format
-(optional model side-by-side plotting)

The results files should then be moved into the correct folders in the repository: 
.../main/model-output/ECDC-norrsken_(blue/green)


Model fits (large files) are stored locally only in a "Big data" folder that the 
user needs to create

Model forecasts are saved under Flu-modelling/output/ in folders with respective hub name 

The model meta-data files can be found in: Flu-modelling/metadata/