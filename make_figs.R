rm(list=ls()); cat('\014')

#setup ####
#install packages from CRAN if necessary
package_list <- c('coda','plyr','dplyr','httr','jsonlite','R2jags',
    'tidyr','imputeTS','accelerometry','geoknife','MetaboPlots')
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

source('~/git/streampulse/model/gapfill_functions.R')
source('~/git/streampulse/model/sp_functions.R')

model_type = "mle"
model_name = "streamMetabolizer"
fillgaps='interpolation'
interval='15 min'
site_code = "AZ_LV"
start_date = "2017-07-07"
end_date = "2017-12-25"
site_code = "AZ_OC"
start_date = "2016-11-13"
end_date = "2017-12-03"
site_code = "NC_Eno"
start_date = "2016-01-01"
end_date = "2017-01-01"
site_code = "NC_Eno"
start_date = "2016-09-30"
end_date = "2016-12-30"
streampulse_data = request_data(sitecode=site_code,
    startdate=start_date, enddate=end_date, variables=NULL,
    flags=FALSE, token=NULL)
# head(streampulse_data)
source('~/git/streampulse/model/gapfill_functions.R')
fitdata = prep_metabolism(d=streampulse_data, type=model_type,
    model=model_name, interval=interval, fillgaps=fillgaps)

modelfit = fit_metabolism(fitdata)

class(fitdata) = "data.frame"
engine = 'stan'; pool_K600='binned'; proc_err = TRUE
modname <- mm_name(type='bayes', pool_K600='binned',
    err_obs_iid=TRUE, err_proc_acor=FALSE, err_proc_iid=TRUE,
    ode_method = 'trapezoid', deficit_src='DO_mod', engine='stan')
modspecs <- specs(modname)
# modspecs$K600_lnQ_nodes_centers <- seq(from=md$logQ_min[n],
#     to = md$logQ_max[n], by = md$Qnodes_steps[n])
modspecs$K600_lnQ_nodes_centers = seq(from=.74, to = 5.25, by = .75)
modelfit <- metab(specs = modspecs, data = fitdata)

# saveRDS(modelfit, paste('~/Dropbox/streampulse/data/models/model_objects/',
    site_code, start_date, end_date, model_type, model_name,
    substr(interval,1,2), fillgaps, sep='_'))

predictions = predict_metabolism(modelfit)


#plot ####
ymin = min(c(predictions$ER,predictions$GPP), na.rm=TRUE)
ymax = max(c(predictions$ER,predictions$GPP), na.rm=TRUE)
ind = which(!is.na(predictions$GPP))
plot(predictions$GPP[ind], type='l', ylim=c(ymin, ymax),
    main='streamMetabolizer', xlab='Date', xaxt='n', ylab='g O2/m^2/d', las=1)
axis(1, seq(1,92,10), predictions$date[ind][seq(1,92,10)])
lines(predictions$ER[ind], col='red')
legend('topright', legend=c('GPP','ER'), lty=1, col=c('black','red'), bty='n')

#phil's package
average_plot(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')
cumulative_plots(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')
# diag_plots(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')
kernel_plot(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')
# p = processing_func(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')
# season_ts_func(p, 'GPP', 'ER')


