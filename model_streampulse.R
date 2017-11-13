install.packages("streamMetabolizer", dependencies=TRUE,
    repos=c("https://owi.usgs.gov/R","https://cran.rstudio.com")) #stable
remove.packages('streamMetabolizer')
library(devtools)
devtools::install_github("USGS-R/streamMetabolizer", ref="develop") #dev

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
model_type <- "mle"
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




par(mfrow=c(1,2))

# pstan = predictions
ymin = min(c(pstan$ER,pstan$GPP), na.rm=TRUE)
ymax = max(c(pstan$ER,pstan$GPP), na.rm=TRUE)
plot(pstan$GPP, type='l', ylim=c(ymin,ymax))
lines(pstan$ER, col='red')
# saveRDS(pjags, '~/git/streampulse/model/temp/pstan.rds')

# pjags = predictions
ymin = min(c(pjags$ER.mean,pjags$GPP.mean), na.rm=TRUE)
ymax = max(c(pjags$ER.mean,pjags$GPP.mean), na.rm=TRUE)
plot(pjags$GPP.mean, type='l', ylim=c(ymin,ymax))
lines(pjags$ER.mean, col='red')
# saveRDS(pjags, '~/git/streampulse/model/temp/pjags.rds')

# x = prep_metabolism('NC_Eno', '2016-11-13', '2017-09-10')
# x = prep_metabolism('OC', '2016-01-01', '2017-03-01')
# sitecode='OC'; startdate='2017-06-13'; enddate='2017-09-10'
# variables=''
