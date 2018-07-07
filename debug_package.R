library(dplyr)
library(accelerometry)
library(imputeTS)
library(zoo)
library(tidyr)
library(tibble)
library(geoknife)
library(readr)
library(geosphere)
library(streamMetabolizer)

detach('package:StreamPULSE', unload=TRUE)
setwd('~/git/streampulse/model/StreamPULSE/R/')
fs = list.files()
for(i in fs) source(i)

#now go to hax/R/R_packages/package_write.R if ready to test for reals
library(StreamPULSE)
library(streamMetabolizer)


#testing model comparison and API upload infrastructure ####
x = request_data('NC_Eno', '2017-01-19', '2017-02-05',
    token='901bb597c597e3237d00')

x$specs

y = prep_metabolism(x)
y$specs

o = fit_metabolism(y)

# saveRDS(fff, '~/Desktop/test_fit.rds')
# saveRDS(ppp, '~/Desktop/test_pred.rds')

model_fit = readRDS('~/Desktop/test_fit.rds')
predictions = readRDS('~/Desktop/test_pred.rds')
mod_startyr = '2017'
d = y
output = list(predictions=predictions, fit=model_fit)
deets = extract_model_details(model_fit, predictions, d$specs)

output = modelfit
predictions = modelfit$predictions
model_fit = modelfit$fit
d = fitdata
deets = StreamPULSE:::extract_model_details(model_fit, predictions, d$specs)
