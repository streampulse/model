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
sp_functions <- GET("https://raw.githubusercontent.com/streampulse/model/master/sp_functions.R")
eval(parse(text = content(sp_functions, as="text", encoding="UTF-8")), envir= .GlobalEnv)

# Model type for streamMetabolizer
# We recommend the Bayesian model, but you can also fit "mle", which runs much faster.
model_type <- "bayes"
# Which modeling framework to use
# "streamMetabolizer" is default, can also use "BASE"
model_name <- "BASE"

# Download data from streampulse and prepare for metabolism modeling
fitdata <- prep_metabolism(sitecode = "NC_Eno",
    startdate = "2016-01-01", enddate = "2017-01-01",
    type = model_type, model = model_name, fillgaps = TRUE)

# Fit metabolism model
#  - if using streamMetabolizer, returns a metab model object
#  - if using BASE, returns a temporary directory where results are saved
modelfit <- fit_metabolism(fitdata)

# Gather metabolism predictions
predictions <- predict_metabolism(modelfit)
