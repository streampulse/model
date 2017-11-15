#setup ####
#install packages from CRAN if necessary
package_list <- c('coda','dplyr','httr','jsonlite','R2jags','tidyr')
new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

#install streamMetabolizer from github (development version)
if (!require("streamMetabolizer")) {
    if (!require("devtools")) install.packages('devtools', repos="http://cran.rstudio.com/")
    library(devtools)
    install_github('USGS-R/streamMetabolizer', ref='develop') #install from github
    detach('package:devtools', unload=TRUE)
}

#make sure streamMetabolizer is up to date in any case
update.packages(oldPkgs=c("streamMetabolizer","unitted"),
    dependencies=TRUE, repos=c("https://owi.usgs.gov/R",
        "https://cran.rstudio.com"))









#compare stan and jags ####
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

#bollocks ####
# x = prep_metabolism('NC_Eno', '2016-11-13', '2017-09-10')
# x = prep_metabolism('OC', '2016-01-01', '2017-03-01')
# sitecode='OC'; startdate='2017-06-13'; enddate='2017-09-10'
# variables=''
