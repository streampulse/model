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


#debug request_data####
sitecode='NC_Eno'; startdate='2017-03-26'; enddate='2017-04-05'; token=NULL
variables = c('DO_mgL','DOsat_pct','satDO_mgL','WaterPres_kPa',
    'Depth_m','WaterTemp_C','Light_PAR','AirPres_kPa')#,'Discharge_m3s')

#debug prep_metabolism ####
d=streampulse_data; model="streamMetabolizer"; type="bayes"
interval='15 min'; rm_flagged=list('Bad Data', 'Questionable')
fillgaps='interpolation'; maxhours=3
zq_curve=list(sensor_height=NULL, Z=NULL, Q=NULL, a=NULL, b=NULL,
    fit='power', ignore_oob_Z=TRUE, plot=TRUE)
estimate_areal_depth=FALSE

#debug query_available_results####
region='all'; site=NULL; year=NULL
region='CT'; site=NULL; year=NULL
region='CT'; site='BUNN'; year=NULL
region='NC'; year=2017; site=NULL
year=2017; region=NULL; site=NULL
query_available_results('NC', 'Eno', '2017')
query_available_results('NC', 'Eno')
query_available_results('FL')
query_available_results('all')
query_available_results('NC', NULL, 2017)
query_available_results(NULL, NULL, 2017)

#debug query_available_data####
region='all'; site=NULL
region='CT'; site=NULL
region='NC'; site='Eno'
region='SE'; site='AbiskoM1'
region=NULL; site=NULL
startdate='2017-03-26'; enddate='2017-04-05'
startdate='2017-11-26'; enddate='2018-04-05'
startdate = enddate = NULL
variable=NULL
variable='DO_mgL'
query_available_data('NC', 'Eno', '2017-01-08', '2017-04-09')
query_available_data('AZ', 'LV', NULL, NULL, 'Depth_m')
query_available_data(NULL, NULL, NULL, NULL, 'Depth_m')
query_available_data('all', NULL, NULL, NULL, NULL)
query_available_data(NULL, NULL, NULL, NULL, 'all')
query_available_data('AZ', 'LV', NULL, NULL, NULL)
query_available_data(NULL, NULL, '2017-04-07', '2017-04-10', NULL)
query_available_data('all', 'LV', '2017-04-07', '2017-04-10', NULL)
query_available_data(NULL, NULL, '2017-04-08', '2017-04-09', 'DO_mgL')

#debug request_results ####
sitecode = 'NC_Eno'
year = '2017'
z = request_results('NC_Eno', '201e')
z = request_results('NC_Eno', '2017')
z = request_results('KS_KANSASR', '2018')
z$data_daily$date[1]

#testing model comparison and API upload infrastructure ####
x = request_data('NC_Eno', '2017-01-26', '2017-02-05')
x = request_data('AZ_OC', '2017-01-01', '2017-12-31')

x$specs

y = prep_metabolism(x, maxhours=3)
fitdata = y
y$specs

o = fit_metabolism(y)

# saveRDS(o, '~/Desktop/test_out.rds')
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
