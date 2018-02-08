rm(list=ls()); cat('\014') #clear environment and console

#setup ####
#install packages from CRAN if necessary
package_list <- c('coda','plyr','dplyr','httr','jsonlite','R2jags',
    'tidyr','imputeTS','accelerometry','geoknife','MetaboPlots','zoo')
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

#choices ####

model_type = "bayes"
# model_type = "mle"
model_name = "streamMetabolizer"
fillgaps='interpolation'
interval='15 min'
#done, ##problem
##site_code = "PR_QS"; start_date = "2014-03-10"; end_date = '2015-03-10'#ed "2017-12-20"
# site_code = "RI_CorkBrk"; start_date = "2015-06-23"; end_date = '2016-06-22'#ed "2017-01-03", sd=2014-06-23
#site_code = "CT_BUNN"; start_date = "2015-05-20"; end_date = '2016-05-20'#ed "2016-11-15"
#site_code = "CT_HUBB"; start_date = "2015-05-20"; end_date = '2016-05-20'#ed "2016-11-11"
# site_code = "CT_FARM"; start_date = "2015-05-20"; end_date = '2016-05-20'#ed "2016-11-11"
#site_code = "CT_Unio"; start_date = "2015-05-20"; end_date = '2016-05-20'#ed "2016-11-11"
# site_code = "VT_Pass"; start_date = "2015-06-26"; end_date = '2016-06-25' #ed 2016-10-27
#site_code = "VT_POPE"; start_date = "2015-06-03"; end_date = '2016-06-02' #ed 2016-11-11
#site_code = "VT_MOOS"; start_date = "2015-06-03"; end_date = '2016-11-11' #ed 2016-11-11
##site_code = "PR_Icacos"; start_date = "2016-06-09"; end_date = '2016-12-15'
#site_code = "CT_STIL"; start_date = "2015-05-20"; end_date = '2016-05-20'#ed "2016-11-15"
#site_code = "MD_GFVN"; start_date = "2016-02-18"; end_date = '2016-11-19'
#site_code = "MD_GFCP"; start_date = "2016-05-07"; end_date = '2016-11-19'
#site_code = "MD_GFGB"; start_date = "2016-02-19"; end_date = '2016-11-11'
#site_code = "MD_DRKR"; start_date = "2016-02-18"; end_date = '2016-11-19'
# site_code = "MD_POBR"; start_date = "2016-02-19"; end_date = '2016-11-19'
#site_code = "MD_BARN"; start_date = "2016-02-19"; end_date = '2016-11-19'
#site_code = "VT_SLPR"; start_date = "2015-06-03"; end_date = '2016-06-03' #ed 2016-11-11
# site_code = "AZ_LV"; start_date = "2017-08-07"; end_date = "2017-12-25" #sd 2017-07-07 but no o2
#site_code = "AZ_OC"; start_date = "2016-11-15"; end_date = "2017-12-03" #sd 11-13
#site_code = "AZ_SC"; start_date = "2017-02-08"; end_date = "2017-03-28"
# site_code = "AZ_WB"; start_date = "2017-08-04"; end_date = "2017-12-27"
# site_code = "NC_Eno"; start_date = "2016-07-11"; end_date = "2017-08-30"
##site_code = "NC_UEno"; start_date = "2016-07-12"; end_date = "2017-08-30"
##site_code = "NC_Stony"; start_date = "2016-06-30"; end_date = "2017-08-09"
site_code = "NC_NHC"; start_date = "2016-09-14"; end_date = "2017-09-13"
##site_code = "NC_UNHC"; start_date = "2016-07-12"; end_date = "2017-08-30"
##site_code = "NC_Mud"; start_date = "2016-07-12"; end_date = "2017-08-30"
# site_code = 'WI_BEC'; start_date = '2017-01-26'; end_date = '2018-01-25' #2009-10-02; 2018-01-25
# site_code = 'WI_BRW'; start_date = '2017-01-27'; end_date = '2018-01-26' #2014-06-13; 2018-01-26

#run ####
source('~/git/streampulse/model/sp_functions.R')
# source('~/git/streampulse/model/gapfill_functions.R')
streampulse_data = request_data(sitecode=site_code,
    startdate=start_date, enddate=end_date, variables=NULL,
    flags=TRUE, token=NULL)
head(streampulse_data$data)
unique(streampulse_data$data$variable)
sum(streampulse_data$data$flagtype != '')
# z = streampulse_data$data
# z = z$DateTime_UTC[z$variable == 'WaterTemp_C']
# z = z$DateTime_UTC[z$variable == 'DOsat_pct']
# unique(substr(z, 15, 16))

plotvars = unique(streampulse_data$data$variable)
par(mfrow=c(length(plotvars),1), mar=c(0,0,0,0), oma=c(4,4,.5,.5))
t = as.Date(sort(unique(streampulse_data$data$DateTime_UTC)))
for(i in plotvars){
    plot(streampulse_data$data$value[which(streampulse_data$data$variable == i)],
        type='l', xlab='', xaxt='n', xaxs='i', las=2)
    mtext(i, 2, 2.5)
    if(i == plotvars[length(plotvars)]){
        yearstarts = match(unique(substr(t,1,4)), substr(t,1,4))
        monthstarts = match(unique(substr(t,1,7)), substr(t,1,7))
        axis(1, yearstarts, substr(t[yearstarts],1,4), line=1, tick=FALSE)
        axis(1, monthstarts, substr(t[monthstarts],6,7))
    }
}

# sval = streampulse_data$data$value
# svar = streampulse_data$data$variable
# dis_vals = which(svar %in% c('Discharge_m3s', 'USGSDischarge_m3s'))
# sval[dis_vals][sval[dis_vals] <= 0] = NA
# streampulse_data$data$value = sval

# dim(streampulse_data$data)
# source('~/git/streampulse/model/sp_functions.R')
# source('~/git/streampulse/model/gapfill_functions.R')
fitdata = prep_metabolism(d=streampulse_data, type=model_type,
    model=model_name, interval='15 min',
    # model=model_name, interval=interval,
    rm_flagged=list('Bad Data', 'Questionable'), fillgaps=fillgaps)

# plot(fitdata$DO.sat, type='l')
# plot(fitdata$DO.sat, type='l', xlim=c(3625,3645), xaxs='i')
# fitdata$DO.sat[3622:3645] = mean(fitdata$DO.sat)
#
# plot(fitdata$DO.obs, type='l')
# plot(fitdata$DO.obs, type='l', xlim=c((3798-191),(3897-186)), xaxs='i')
# fitdata$DO.obs[(3798-191):(3897-186)] = 11


plotvars = colnames(fitdata)[! colnames(fitdata) %in% c('solar.time')]
# pdf(width=5, height=9,
#     file=paste0('~/Desktop/untracked/sm_figs/input_',
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

modelfit = fit_metabolism(fitdata)

max(modelfit@fit$daily$K600_daily_mean, na.rm=TRUE) #shoudnt be over 90
plot(density(modelfit@fit$daily$K600_daily_mean, na.rm=TRUE))

# class(fitdata) = "data.frame"
# engine = 'stan'; pool_K600='binned'; proc_err = TRUE
# modname <- mm_name(type='bayes', pool_K600='binned',
#     err_obs_iid=TRUE, err_proc_acor=FALSE, err_proc_iid=TRUE,
#     ode_method = 'trapezoid', deficit_src='DO_mod', engine='stan')
# modspecs <- specs(modname)
# addis = tapply(log(fitdata$discharge),
#     substr(fitdata$solar.time,1,10), mean)
# modspecs$K600_lnQ_nodes_centers = seq(from=min(addis),
#     to=max(addis), length.out=7)
# modelfit <- metab(specs = modspecs, data = fitdata)

saveRDS(modelfit, paste('~/Desktop/untracked/sm_out/fit',
    site_code, start_date, end_date,
    'bayes_binned_obsproc_trapezoid_DO-mod_stan.rds', sep='_'))
    # model_type, model_name, substr(interval,1,2), fillgaps, sep='_'))
# modelfit = readRDS(paste0('~/Desktop/untracked/sm_out/fit_AZ_SC_2017-02-08_2017-03-28_bayes_binned_obsproc_trapezoid_DO-mod_stan.rds'))
# modelfit = readRDS(paste0('~/Desktop/untracked/sm_out/fit_AZ_LV_2017-08-07_2017-12-25_bayes_binned_obsproc_trapezoid_DO-mod_stan.rds'))

#check daiy k-er correlation
daily_er = modelfit@fit$daily$ER_daily_mean
daily_k = modelfit@fit$daily$K600_daily_mean
plot(daily_k, daily_er)
daily_k[is.na(daily_k)] = mean(daily_k, na.rm=TRUE)
daily_er[is.na(daily_er)] = mean(daily_er, na.rm=TRUE)
cor(daily_k, daily_er)

predictions = predict_metabolism(modelfit)

saveRDS(predictions, paste('~/Desktop/untracked/sm_out/predictions',
    site_code, start_date, end_date,
    'bayes_binned_obsproc_trapezoid_DO-mod_stan.rds', sep='_'))

#plot ####
# ymin = min(c(predictions$ER,predictions$GPP), na.rm=TRUE)
# ymax = max(c(predictions$ER,predictions$GPP), na.rm=TRUE)
# ind = which(!is.na(predictions$GPP))
# plot(predictions$GPP[ind], type='l', ylim=c(ymin, ymax),
#     main='streamMetabolizer', xlab='Date', xaxt='n', ylab='g O2/m^2/d', las=1)
# axis(1, seq(1,92,10), predictions$date[ind][seq(1,92,10)])
# lines(predictions$ER[ind], col='red')
# legend('topright', legend=c('GPP','ER'), lty=1, col=c('black','red'), bty='n')

#phil's package
# ts = predictions[,c('date','GPP','ER')]
# date_var='date'; gpp_var='GPP'; er_var='ER'
processing_func = function (ts, date_var, gpp_var, er_var) {
    ts[!is.na(ts[, gpp_var]) & ts[, gpp_var] < 0 | !is.na(ts[,
        gpp_var]) & ts[, gpp_var] > 100, gpp_var] <- NA
    ts[!is.na(ts[, er_var]) & ts[, er_var] > 0, "er_var"] <- NA
    if('tbl' %in% class(ts[,date_var])){
        full_dates = as.data.frame(ts[,date_var])
        colnames(full_dates) = 'Date'
        colnames(ts)[which(colnames(ts) == 'date')] = 'Date'
    } else {
        ts_date = ts[,date_var]
        ts$Date <- format(ts_date, "%Y-%m-%d")
        full_dates <- setNames(data.frame(
            seq(from = as.Date(paste(min(format(ts[, date_var], "%Y")),
            "-01-01", sep = "")),
            to = as.Date(paste(max(format(ts[,
            date_var], "%Y")), "-12-31", sep = "")), by = 1)), "Date")
    }
    ts_full <- merge(full_dates, ts, by = c("Date"), all = TRUE)
    ts_full$Year <- as.numeric(format(ts_full[, "Date"], "%Y"))
    ts_full$DOY <- as.numeric(format(ts_full[, "Date"], "%j"))
    ts_full$NPP <- ts_full[, gpp_var] + ts_full[, er_var]
    return(ts_full)
}
p = processing_func(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')

# pdf(width=5, height=9,
#     file=paste('~/Desktop/untracked/sm_figs/cumulative',
#     site_code, start_date, end_date, sep='_'), compress=FALSE)
# cumulative_func(p, 'GPP', 'ER')
# dev.off()
#
# pdf(width=8, height=6,
#     file=paste('~/Desktop/untracked/sm_figs/ts',
#     site_code, start_date, end_date, sep='_'), compress=FALSE)
# season_ts_func(p, 'GPP', 'ER')
# dev.off()
#
# pdf(width=6, height=9,
#     file=paste('~/Desktop/untracked/sm_figs/kernel',
#     site_code, start_date, end_date, sep='_'), compress=FALSE)
# kernel_func(p, 'GPP', 'ER')
# dev.off()

diag_plots = function (ts, date_var, gpp_var, er_var){
    ts_full <- processing_func(ts, date_var, gpp_var, er_var)
    layout(matrix(c(1, 1, 3, 1, 1, 4, 2, 2, 5), 3, 3, byrow = TRUE),
        widths = c(1, 1, 2))
    par(cex = 0.6)
    par(mar = c(3, 4, 0.1, 0.1), oma = c(3, 0.5, 0.5, 0.5))
    par(tcl = -0.25)
    par(mgp = c(2, 0.6, 0))
    kernel_func(ts_full, gpp_var, er_var)
    par(mar = c(0, 4, 2, 0.1))
    season_ts_func(ts_full, gpp_var, er_var)
    par(mar = c(0, 4, 0, 0.1))
    cumulative_func(ts_full, gpp_var, er_var)
}
pdf(width=7, height=7,
    file=paste0('~/Desktop/untracked/sm_figs/output2_',
    site_code, '_', start_date, '_', end_date, '.pdf'), compress=FALSE)
diag_plots(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')
dev.off()

# average_plot(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')
# cumulative_plots(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')
# kernel_plot(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')
