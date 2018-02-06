#StreamPULSE model script
#all sites, round 1
#Mike Vlah: vlahm13@gmail.com
#2/5/18

#if you're using RStudio, you can collapse sections with
#Alt+o (Windows/Linux) or Cmd+Opt+o (Mac)

rm(list=ls()); cat('\014') #clear environment and console
setwd('MODIFY_PATH_HERE/streampulse_model_runs_round1') #set working directory

#install/update streamMetabolizer ####

#install latest version if you don't have streamMetabolizer already
install.packages("streamMetabolizer", dependencies=TRUE,
    repos=c("https://owi.usgs.gov/R","https://cran.rstudio.com"))

#update streamMetabolizer if you aready have it
update.packages(oldPkgs=c("streamMetabolizer","unitted"), dependencies=TRUE,
    repos=c("https://owi.usgs.gov/R", "https://cran.rstudio.com"))

#load packages and StreamPULSE functions ####

library(dplyr)
library(httr)
library(jsonlite)
library(streamMetabolizer)
library(tidyr)
library(imputeTS)
library(zoo)
library(accelerometry)
library(geoknife) #for acquiring air pressure if necessary
# library(R2jags) # only required for BASE
# library(coda) # only required for BASE

# This sources our StreamPULSE functions from GitHub (R package coming soon)
sp_functions = GET("https://raw.githubusercontent.com/streampulse/model/master/sp_functions.R")
eval(parse(text=content(sp_functions, as="text", encoding="UTF-8")), envir=.GlobalEnv)

# Source gapfilling functions
gapfill_functions = GET("https://raw.githubusercontent.com/streampulse/model/master/gapfill_functions.R")
eval(parse(text = content(gapfill_functions, as="text", encoding="UTF-8")),
    envir= .GlobalEnv)

# Source BASE functions (not needed for this script)
# BASE_functions = GET("https://raw.githubusercontent.com/streampulse/model/master/BASE_functions.R")
# eval(parse(text = content(BASE_functions, as="text", encoding="UTF-8")),
#     envir= .GlobalEnv)

#model setup (see https://data.streampulse.org/model for details) ####

model_type = 'bayes'
model_name = 'streamMetabolizer'
fillgaps = 'interpolation'
interval = '15 min'

#choose site and dates (commented sites still have issues,
#which are detailed in notes_from_model_runs.csv).
##site_code = "PR_QS"; start_date = "2014-03-10"; end_date = '2015-03-10'#"2017-12-20"
##site_code = "PR_Icacos"; start_date = "2016-06-09"; end_date = '2016-12-15'
#site_code = "CT_BUNN"; start_date = "2015-05-20"; end_date = '2016-05-20'#ed "2016-11-15"
#site_code = "CT_STIL"; start_date = "2015-05-20"; end_date = '2016-05-20'#ed "2016-11-15"
#site_code = "MD_GFVN"; start_date = "2016-02-18"; end_date = '2016-11-19'
# site_code = "MD_GFCP"; start_date = "2016-05-07"; end_date = '2016-11-19'
# site_code = "MD_GFGB"; start_date = "2016-02-19"; end_date = '2016-11-11'
# site_code = "MD_DRKR"; start_date = "2016-02-18"; end_date = '2016-11-19'
# site_code = "MD_POBR"; start_date = "2016-02-19"; end_date = '2016-11-19'
# site_code = "MD_BARN"; start_date = "2016-02-19"; end_date = '2016-11-19'
# site_code = "VT_SLPR"; start_date = "2015-06-03"; end_date = '2016-06-03' #ed 2016-11-11
#site_code = "AZ_LV"; start_date = "2017-08-07"; end_date = "2017-12-25" #sd 2017-07-07 but no o2
#site_code = "AZ_OC"; start_date = "2016-11-15"; end_date = "2017-12-03" #sd 11-13
#site_code = "NC_Eno"; start_date = "2016-07-11"; end_date = "2017-08-30"
##site_code = "NC_UEno"; start_date = "2016-07-12"; end_date = "2017-08-30"
##site_code = "NC_Stony"; start_date = "2016-06-30"; end_date = "2017-08-09"
##site_code = "NC_NHC"; start_date = "2016-09-14"; end_date = "2017-09-13"
##site_code = "NC_UNHC"; start_date = "2016-07-12"; end_date = "2017-08-30"
##site_code = "NC_Mud"; start_date = "2016-07-12"; end_date = "2017-08-30"

#retrieve data from StreamPULSE database
streampulse_data = request_data(sitecode=site_code,
    startdate=start_date, enddate=end_date, variables=NULL,
    flags=TRUE, token=NULL)

#acquire additional data if necessary; format for streamMetabolizer
fitdata = prep_metabolism(d=streampulse_data, type=model_type,
    model=model_name, interval=interval,
    rm_flagged=list('Bad Data', 'Questionable'), fillgaps=fillgaps)

#model ####

#see what you're dealing with (uncomment lines below to save pdf.
#make sure working directory is set at top of script if you do.)
#note that existing plots will overwritten unless you move or rename them.
plotvars = colnames(fitdata)[! colnames(fitdata) %in% c('solar.time')]
# pdf(width=5, height=9, file=paste0('figs/input_',
#     site_code, '_', start_date, '_', end_date, '.pdf'), compress=FALSE)
par(mfrow=c(length(plotvars),1), mar=c(0,0,0,0), oma=c(4,4,.5,.5))
t = as.Date(fitdata$solar.time)
for(i in plotvars){
    plot(fitdata[,i], type='l', xlab='', xaxt='n', xaxs='i', las=2)
    mtext(i, 2, 2.5)
    if(i == plotvars[length(plotvars)]){
        yearstarts = match(unique(substr(t,1,4)), substr(t,1,4))
        monthstarts = match(unique(substr(t,1,7)), substr(t,1,7))
        axis(1, yearstarts, substr(t[yearstarts],1,4), line=1, tick=FALSE)
        axis(1, monthstarts, substr(t[monthstarts],6,7))
    }
}
# dev.off()

#fit model with stock parameters (uncomment next block for more options)
modelfit = fit_metabolism(fitdata)

#fit model with custom parameters (mm_name, specs, metab are from streamMetabolizer)
# class(fitdata) = "data.frame" #just a formality, execute and disregard
# modname = mm_name(type='bayes', pool_K600='binned',
#     err_obs_iid=TRUE, err_proc_acor=FALSE, err_proc_iid=TRUE,
#     ode_method = 'trapezoid', deficit_src='DO_mod', engine='stan')
# modspecs = specs(modname)
#
# #overwrite default ln(Q) node centers
# addis = tapply(log(fitdata$discharge), substr(fitdata$solar.time,1,10), mean)
# modspecs$K600_lnQ_nodes_centers = seq(from=min(addis),
#     to=max(addis), length.out=7)
#
# modelfit = metab(specs=modspecs, data=fitdata)

#check daily k estimates (k>90 is cause for concern)
max(modelfit@fit$daily$K600_daily_mean, na.rm=TRUE) #max daily k
plot(density(modelfit@fit$daily$K600_daily_mean, na.rm=TRUE)) #distribution

#save streamMetabolizer output if desired
saveRDS(modelfit, paste('mod_objects/fit', site_code, start_date, end_date,
    'bayes_binned_obsproc_trapezoid_DO-mod_stan.rds', sep='_'))

#read streamMetabolizer output back in later (uncomment and specify path)
# modelfit = readRDS(paste0('mod_objects/fit_....rds'))

#check estimated daily k-ER correlation
daily_er = modelfit@fit$daily$ER_daily_mean
daily_k = modelfit@fit$daily$K600_daily_mean
plot(daily_k, daily_er) #visual evaluation
daily_k[is.na(daily_k)] = mean(daily_k, na.rm=TRUE)
daily_er[is.na(daily_er)] = mean(daily_er, na.rm=TRUE)
cor(daily_k, daily_er) #Pearson coefficient (> 0.6 is cause for concern)

#explore output ####

#get model predictions
predictions = predict_metabolism(modelfit)

#save predictions if desired
#note that existing prediction objects will overwritten unless
#you move or rename them.
saveRDS(predictions, paste('mod_objects/predictions',
    site_code, start_date, end_date,
    'bayes_binned_obsproc_trapezoid_DO-mod_stan.rds', sep='_'))

#source functions from Phil's forthcoming MetaboPlots package
source('MetaboPlots_prerelease.R')

#plot metabolism, cumulative GPP-ER, GPP x ER kernel density
#uncomment lines below to save plots. note that existing plots will be
#overwritten unless you move or rename them.
p = processing_func(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')
# pdf(width=7, height=7, file=paste0('figs/output_',
#     site_code, '_', start_date, '_', end_date, '.pdf'), compress=FALSE)
diag_plots(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')
# dev.off()

#if the above plot doesn't work, this can get you started
ymin = min(c(predictions$ER,predictions$GPP), na.rm=TRUE)
ymax = max(c(predictions$ER,predictions$GPP), na.rm=TRUE)
plot(predictions$GPP, type='l', ylim=c(ymin, ymax),
    main='', xlab='Date Index', ylab='g O2/m^2/d', las=1)
lines(predictions$ER[ind], col='red')
legend('topright', legend=c('GPP','ER'), lty=1, col=c('black','red'), bty='n')
