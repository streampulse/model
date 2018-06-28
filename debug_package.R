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

#now go to hax/R/package_write.R if necessary

x = request_data('NC_Eno', '2017-02-01', '2017-02-15',
    token='901bb597c597e3237d00')

x$specs

y = prep_metabolism(x)
y$specs

o = fit_metabolism(y)
