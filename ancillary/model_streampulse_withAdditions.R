#Basic StreamPULSE data processing pipeline (with additions for plotting,
#outlier detection, time series decomposition, drift correction)
#Updated 8/31/18

#SETUP####

# Install StreamPULSE pipeline tools from GitHub
# This package is currently in development and changes frequently!
# If something doesn't work as expected, first try reinstalling.
library(devtools)
install_github('streampulse/StreamPULSE', dependencies=TRUE)

# Install latest version of streamMetabolizer.
install.packages('streamMetabolizer', dependencies=TRUE,
    repos=c('https://owi.usgs.gov/R','https://cran.rstudio.com'))

#set this to wherever you extracted the "outliers_etc" folder
setwd('~/Desktop/outliers_etc/')

# Attach packages.
library(StreamPULSE)
library(streamMetabolizer)
library(ks)

#PIPELINE####

# Select site and date range for which to acquire StreamPULSE data.
# site_code is a combination of regionID and siteID
# from https://data.streampulse.org/sitelist
site_code = 'NC_UNHC'
start_date = '2017-01-01'
end_date = '2017-01-14'
site_code = 'AZ_LV'
start_date = '2018-05-10'
end_date = '2018-06-06'

# Download data from streampulse.
# Token parameter not needed if downloading public data.
# If you need a token, contact Mike at streampulse.info@gmail.com
streampulse_data = request_data(sitecode=site_code,
    startdate=start_date, enddate=end_date)

#read in rating curve data and sensor offsets
zq = read.csv('ZQ_data.csv')
offsets = read.csv('sensor_offsets.csv')

site = strsplit(site_code, '_')[[1]][2]
Z = zq[zq$site == site, 'level_m']
Q = zq[zq$site == site, 'discharge_cms']
offset = offsets[offsets$site == site, 2] / 100

# Format data for metabolism modeling.
fitdata = prep_metabolism(d=streampulse_data, type='bayes',
    model='streamMetabolizer', interval='15 min')
fitdata = prep_metabolism(d=streampulse_data, type='bayes',
    model='streamMetabolizer', interval='15 min',
    zq_curve=list(sensor_height=offset, Z=Z, Q=Q, fit='power',
        ignore_oob_Z=TRUE, plot=TRUE))

# Fit metabolism model and generate predictions (calls streamMetabolizer
# functions: mm_name, specs, metab, predict_metab).
modelfit = fit_metabolism(fitdata)
# saveRDS(modelfit, '~/Desktop/modelfit.rds') #store output to avoid rerunning everything
# modelfit = readRDS('~/Desktop/modelfit.rds')

#PLOTS####

#load plotting functions from separate file
source('plotting_functions.R')

#extract some useful data
data_daily = modelfit$fit@data_daily
data = modelfit$fit@data
fit = modelfit$fit@fit
#mod_out is one of the objects available for each siteyear in the zip file i sent Alice
mod_out = list('data_daily'=data_daily, 'data'=data, 'fit'=fit)
#predictions is the other thing available for each siteyear
predictions = modelfit$predictions

#make each of the plots from the diagnostic app (some axes and legends will be missing here)
pr = processing_func(predictions, 1, 366)
#you can adjust the plot window by changing the 1 and 366 (DOY bounds)
season_ts_func(pr, suppress_NEP=TRUE, 1, 366)
cumulative_func(pr, 1, 366)
kernel_func(pr)
O2_plot(mod_out, 1, 366)
KvQvER_plot(mod_out)

#OUTLIERS####
source('~/git/streampulse/server_copy/sp/find_outliers.R')
df = fitdata$data
dd = data.frame('solar.time'=df$solar.time, 'DO.obs'=df$DO.obs,
    'DO.sat'=df$DO.sat, 'depth'=df$depth, 'temp.water'=df$temp.water,
    'light'=df$light, 'discharge'=df$discharge)
df = dd
#get list of outliers detected for each input series
out = outlier_detect(fitdata$data) #pretty sure the warnings are chill (line 259)
out = outlier_detect(dd) #pretty sure the warnings are chill (line 259)
z = read.csv('~/Dropbox/streampulse/data/test_outl2.csv')
source('~/git/streampulse/server_copy/sp/find_outliers2.R')
out = outlier_detect(z)

jump_inds = readRDS('~/Desktop/jump_inds.rds')
pos_jumps = readRDS('~/Desktop/pos_jumps.rds')
as.numeric(jump_inds %in% pos_jumps)
runs = rle(as.numeric(jump_inds %in% pos_jumps))
ends = cumsum(runs$lengths)
runs = cbind(values=runs$values, starts=c(1, lag(ends)[-1] + 1),
    stops=ends, lengths=runs$lengths, deparse.level=1)


#TIME SERIES DECOMPOSITION####

#synthesize a series with trend, periodic, and error components
x = 1:100
tr = seq(1, 20, length.out=100)
pe = sin(tr)
er = rnorm(100, 0, 1)
ts = tr + pe + er
ts = ts(ts, deltat=3/100)

plot(pe)
plot(ts)
dec = decompose(ts)
plot(dec) #this doesnt do a good job of extracting the seasonal trend,
#probably because i didn't specify the period precisely. but there's the
#concept anyway. and then you can extract components like so:
plot(dec$seasonal)

#DRIFT CORRECTION####
#(simple slope-based approach)

#synthesize and plot hypothetical interval between cleanings
x = 1:1000
tr = seq(1, 20, length.out=1000)
er = rnorm(1000, 0, 2)
ts = tr + er
ts = ts(ts)
plot(ts, col='black', type='l')

#fit linear model and add to plot
mod = lm(ts ~ x)
abline(mod, col='red')

#get y values on fitted line for each x
y_predictions = predict(mod, data.frame(x=x))
?predict.lm #for help with the function above

#subtract them out and add corrected series to plot
corrected = ts - y_predictions
lines(corrected, col='blue')
