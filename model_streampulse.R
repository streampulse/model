# You need to have an internet connection for this code to function.
# REQUIRED packages
library(coda)
library(dplyr)
library(httr)
library(jsonlite)
library(R2jags)
library(streamMetabolizer)
library(tidyr)

# REQUIRED code
# These source our StreamPULSE functions from GitHub
# In the future, turn these codes into a package and just source the package...


rm(list=ls()); cat('\014')
source('~/git/streampulse/model/sp_functions.R')
# sp_functions = GET("https://raw.githubusercontent.com/streampulse/model/master/sp_functions.R")
# eval(parse(text=content(sp_functions, as="text", encoding="UTF-8")), envir=.GlobalEnv)

# Model type for streamMetabolizer
# We recommend the Bayesian model ("bayes"), but you can also fit "mle",
# which runs much faster.
model_type = "bayes"

# Which modeling framework to use
# "streamMetabolizer" is default; can also use "BASE"
model_name = "streamMetabolizer"

# Select site and date range
site_code = "NC_Eno" # a combination of the regionID and siteID
start_date = "2016-11-01"
end_date = "2017-01-01"

# Download data from streampulse
streampulse_data = retrieve_data(sitecode=site_code,
    startdate=start_date, enddate=end_date)

# Format data for metabolism modeling
fitdata = prep_metabolism(d=streampulse_data, type=model_type,
    model=model_name, fillgaps=TRUE)

# Fit metabolism model
#  - if using streamMetabolizer, returns a metab model object
#  - if using BASE, returns a temporary directory where results are saved
modelfit = fit_metabolism(fitdata)

# Gather metabolism predictions
predictions = predict_metabolism(modelfit)


