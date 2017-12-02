#setup ####
#install packages from CRAN if necessary
package_list <- c('coda','dplyr','httr','jsonlite','R2jags','tidyr','imputeTS',
    'accelerometry','geoknife')
new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

#install streamMetabolizer from github (stable version)
if (!require("streamMetabolizer")) {
    if (!require("devtools")) install.packages('devtools', repos="http://cran.rstudio.com/")
    library(devtools)
    # install_github('USGS-R/streamMetabolizer', ref='develop') #install from github
    #install specific stable release on github. check the site to find out
    #if the library() text doesnt alert you:
    # https://github.com/USGS-R/streamMetabolizer/releases
    install_github('USGS-R/streamMetabolizer@v0.10.8')
    detach('package:devtools', unload=TRUE)
}

#make sure streamMetabolizer is up to date
# update.packages(oldPkgs=c("streamMetabolizer","unitted"),
#     dependencies=TRUE, repos=c("https://owi.usgs.gov/R",
#         "https://cran.rstudio.com"))

#load all packages
for(i in c(package_list)) library(i, character.only=TRUE)

#experiment ####
rm(list=ls()); cat('\014')
source('~/git/streampulse/model/gapfill_functions.R')
source('~/git/streampulse/model/sp_functions.R')

model_type = "bayes"
model_name = "streamMetabolizer"
site_code = "NC_Eno"
start_date = "2016-01-01"
end_date = "2017-01-01"
streampulse_data = request_data(sitecode=site_code,
    startdate=start_date, enddate=end_date, variables=NULL,
    flags=FALSE, token=NULL)
fitdata = prep_metabolism(d=streampulse_data, type=model_type,
    model=model_name, interval='15 min', fillgaps=TRUE)
modelfit = fit_metabolism(fitdata)
predictions = predict_metabolism(modelfit)


#compare stan and jags (bayes, gpp with er)####

pstan = readRDS('~/git/streampulse/model/temp/pstan.rds')
pjags = readRDS('~/git/streampulse/model/temp/pjags.rds')

pdf(width=7, height=6, compress=FALSE,
    file='~/Dropbox/streampulse/figs/BASE_streamMetab_comparison.pdf')
par(mfrow=c(2,1))

# pstan = predictions
ymin = min(c(pstan$ER,pstan$GPP), na.rm=TRUE)
ymax = max(c(pstan$ER,pstan$GPP), na.rm=TRUE)
ind = which(!is.na(pstan$GPP))
plot(pstan$GPP[ind], type='l', ylim=c(-20,170),#ylim=c(ymin,ymax),
    main='streamMetabolizer', xlab='Date', xaxt='n', ylab='g O2/m^2/d', las=1)
axis(1, seq(1,92,10), pstan$date[ind][seq(1,92,10)])
lines(pstan$ER[ind], col='red')
legend('topright', legend=c('GPP','ER'), lty=1, col=c('black','red'), bty='n')
# saveRDS(pstan, '~/git/streampulse/model/temp/pstan.rds')

# pjags = predictions
# pjags = readRDS('~/git/streampulse/model/temp/pjags.rds')
ymin = min(c(pjags$ER.mean,pjags$GPP.mean), na.rm=TRUE)
ymax = max(c(pjags$ER.mean,pjags$GPP.mean), na.rm=TRUE)
plot(pjags$GPP.mean, type='l', ylim=c(-20,170),#ylim=c(ymin,ymax),
    main='BASE', xlab='Date', xaxt='n', ylab='g O2/m^2/d', las=1)
axis(1, seq(1,92,10), pjags$date[seq(1,92,10)])
lines(pjags$ER.mean, col='red')
legend('topright', legend=c('GPP','ER'), lty=1, col=c('black','red'), bty='n')
# saveRDS(pjags, '~/git/streampulse/model/temp/pjags.rds')
dev.off()

#compare stan and jags (bayes, gpp with gpp, er with er) ####
pdf(width=7, height=6, compress=FALSE,
    file='~/Dropbox/streampulse/figs/BASE_streamMetab_comparison2.pdf')
par(mfrow=c(2,1))

# pstan = predictions
ymin = min(c(pstan$GPP, pjags$GPP.mean), na.rm=TRUE)
ymax = max(c(pstan$GPP, pjags$GPP.mean), na.rm=TRUE)
ind = which(!is.na(pstan$GPP))
plot(pstan$GPP[ind], type='l', ylim=c(ymin,10),#ylim=c(ymin,ymax),
    main='GPP', xlab='Date', xaxt='n', ylab='g O2/m^2/d', las=1, pch=20)
axis(1, seq(1,92,10), pstan$date[ind][seq(1,92,10)])
lines(pjags$GPP.mean, col='red', pch=20)
legend('topright', legend=c('SM','BASE'), lty=1, col=c('black','red'), bty='n')
# saveRDS(pstan, '~/git/streampulse/model/temp/pstan.rds')

# pjags = predictions
# pjags = readRDS('~/git/streampulse/model/temp/pjags.rds')
ymin = min(c(pstan$ER,pjags$ER.mean), na.rm=TRUE)
ymax = max(c(pstan$ER,pjags$ER.mean), na.rm=TRUE)
plot(pstan$ER[ind], type='l', ylim=c(ymin,50),#ylim=c(ymin,ymax),
    main='ER', xlab='Date', xaxt='n', ylab='g O2/m^2/d', las=1, pch=20)
axis(1, seq(1,92,10), pjags$date[seq(1,92,10)])
lines(pjags$ER.mean, col='red', pch=20)
legend('topright', legend=c('SM','BASE'), lty=1, col=c('black','red'), bty='n')
# saveRDS(pjags, '~/git/streampulse/model/temp/pjags.rds')
dev.off()


#compare mle and bayes (TODO) ####

ymin = min(c(predictions$GPP, predictions$ER), na.rm=TRUE)
ymax = max(c(predictions$GPP, predictions$ER), na.rm=TRUE)
plot(predictions$GPP, type='l', ylim=c(ymin, ymax))
lines(predictions$ER, col='red')

#model exploration####
?mm_name()
mm_valid_names('mle')
mm_parse_name(mm_valid_names('mle'))

#get windspeed and air pressure data for streampulse sites ####

library(geoknife)

#load site list
sites = read.csv('/home/mike/Dropbox/streampulse/data/site_data.csv',
    stringsAsFactors=FALSE)
sites = sites[-which(sites$region == 'PR'),]
sites = sites[3,]

#get windspeed and air pressure datasets #(nothing good in the GDP library)
webdatasets = query('webdata')

title(webdatasets)
metdata = webdata(webdatasets[95]) #promising: 77, 95

set_num = grep('speed', abstract(webdatasets), ignore.case=TRUE)
title(webdatasets[set_num])
metdata = webdata(webdatasets[set_num[3]]) #only viable one, but doesnt work for 2016-17?

set_num = grep('pressure', abstract(webdatasets))
title(webdatasets[set_num])
metdata = webdata(webdatasets[set_num[1]])

#format names
sites$name = gsub(' ', '-', sites$name)
stations = as.data.frame(t(sites[,c('lon','lat')]))
colnames(stations) = paste(sites$region, sites$site, sites$name, sep='_')
# stations = data.frame('Eno_River'=c(-79.0968, 36.0715)) #circumvent the above for testing

#set extents
query(metdata, 'variables')
variables(metdata) = 'wind_speed'
query(metdata, 'times')
times(metdata) = as.POSIXct(c("2016-01-01", "2017-01-01"))
stations = simplegeom(stations)

#edit job details
# knife = webprocess()
# knife@wait = TRUE
# knife@processInputs$REQUIRE_FULL_COVERAGE = 'false'
# algs = query(webprocess(), 'algorithms')
# knife@algorithm = algs[4]

#submit job #UI library would do the job, but it's broken
metdata_job = geoknife(stencil=stations, fabric=metdata, wait=TRUE)
metdata_data = result(metdata_job, with.units=TRUE)
check(metdata_job)
head(metdata_data)
# 44.5117, 44.5118, 44.5117, 44.5116, 44.5117
# -71.8379, -71.8378, -71.8377, -71.8378, -71.8379

#instead, lets use noaa
fabric <- webdata(url = 'https://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface/vwnd.sig995.2017.nc',
    variables = "vwnd")
# keep time as c(NA, NA) to get all times
stations = simplegeom(stations)

#edit job details
knife = webprocess(wait=TRUE)

metdata_job = geoknife(stencil=stations, fabric=fabric, knife=knife)
metdata_data = result(metdata_job, with.units=TRUE)
check(metdata_job)
head(metdata_data)

#temp code ####

x = input_data[7800:7850,-(1:2)]
x[20:30,2:3] = NA; x[1:5,4:5] = NA; x[12:36,1] = NA
imputed = sapply(X=x,
    FUN=gap_impute, tol=12, algorithm='interpolation',
    simplify=TRUE)

x3 = gap_impute(x, tol=12)
plot(x3, type='l', col='red')
lines(x, col='black')


