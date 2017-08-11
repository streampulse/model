# You need to have an internet connection for this code to function.
# REQUIRED packages
library(coda)
library(dplyr)
library(httr)
library(jsonlite)
library(R2jags)
library(streamMetabolizer)
library(tidyr)
library(zoo)

# REQUIRED code
# This sources our StreamPULSE functions from GitHub.
spfunctions <- GET("https://raw.githubusercontent.com/berdaniera/StreamPulsev2/master/sp_functions.R")
eval(parse(text = content(spfunctions,as="text",encoding="UTF-8")), envir= .GlobalEnv)

# Model type for streamMetabolizer
# We recommend the Bayesian model, but you can also fit "mle", which runs much faster.
model_type <- "bayes"
# Which modeling framework to use
# "streamMetabolizer" is default, can also use "BASE"
model_name <- "streamMetabolizer"

# Get StreamPULSE data for metabolism modeling
fitdata <- sp_data_metab(sitecode = "NC_Eno",
    startdate = "2016-07-01", enddate = "2017-07-01",
    type = model_type, model = model_name, fillgaps=TRUE)

# Fit models
predictions <- fit_metabolism(fitdata, model_name, model_type)
