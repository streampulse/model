# Model
Streampulse metabolism modeling code

# Workflow
Fitting metabolism models with streampulse follows two steps:

1. Downloading data from the streampulse platform
2. Fitting the chosen model
  * `model_name = "streamMetabolizer"` runs [streamMetabolizer](https://github.com/USGS-R/streamMetabolizer)
    - `model_type="bayes"` fits a Bayesian model (recommended)
    - `model_type="mle"` fits a maximum likelihood model
  * `model_name = "BASE"` runs [BASE](https://github.com/dgiling/BASE)

# Code
```r
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
# These source our StreamPULSE functions from GitHub
# In the future, turn these codes into a package and just source the package...
sp_functions <- GET("https://raw.githubusercontent.com/streampulse/model/master/sp_functions.R")
eval(parse(text = content(sp_functions, as="text", encoding="UTF-8")), envir= .GlobalEnv)
gapfill_functions <- GET("https://raw.githubusercontent.com/streampulse/model/master/gapfill_functions.R")
eval(parse(text = content(gapfill_functions, as="text", encoding="UTF-8")), envir= .GlobalEnv)
BASE_functions <- GET("https://raw.githubusercontent.com/streampulse/model/master/BASE_functions.R")
eval(parse(text = content(BASE_functions, as="text", encoding="UTF-8")), envir= .GlobalEnv)

# Which modeling framework to use
# "streamMetabolizer" is default, can also use "BASE"
model_name <- "streamMetabolizer"
# Model type for streamMetabolizer
# We recommend the Bayesian model, but you can also fit "mle", which runs much faster.
model_type <- "bayes"

# Get StreamPULSE data for metabolism modeling
fitdata <- sp_data_metab(sitecode = "NC_Eno",
    startdate = "2016-01-01", enddate = "2017-01-01",
    type = model_type, model = model_name, fillgaps = TRUE)

# Fit models
predictions <- fit_metabolism(fitdata, model_name, model_type)
```
