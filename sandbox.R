#setup ####
#install packages from CRAN if necessary
package_list <- c('coda','dplyr','httr','jsonlite','R2jags','tidyr','imputeTS',
    'accelerometry')
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
rm(list=ls()); cat('/014')
source('~/git/streampulse/model/gapfill_functions.R')
source('~/git/streampulse/model/sp_functions.R')

model_type = "mle"
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

#temp code ####

x = input_data[7800:7850,-(1:2)]
x[20:30,2:3] = NA; x[1:5,4:5] = NA; x[12:36,1] = NA
imputed = sapply(X=x,
    FUN=gap_impute, tol=12, algorithm='interpolation',
    simplify=TRUE)

x3 = gap_impute(x, tol=12)
plot(x3, type='l', col='red')
lines(x, col='black')


