rm(list=ls()); cat('\014')

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

model_type = "bayes"
# model_type = "mle"
model_name = "streamMetabolizer"
fillgaps='interpolation'
interval='15 min'
#done, ##problem
##site_code = "PR_QS"; start_date = "2014-03-10"; end_date = '2015-03-10'#"2017-12-20"
#site_code = "AZ_LV"; start_date = "2017-08-07"; end_date = "2017-12-25" #sd 2017-07-07 but no o2
site_code = "AZ_OC"; start_date = "2016-11-15"; end_date = "2017-12-03" #sd 11-13
#site_code = "NC_Eno"; start_date = "2016-07-11"; end_date = "2017-08-30"
##site_code = "NC_UEno"; start_date = "2016-07-12"; end_date = "2017-08-30"
##site_code = "NC_Stony"; start_date = "2016-06-30"; end_date = "2017-08-09"
##site_code = "NC_NHC"; start_date = "2016-09-14"; end_date = "2017-09-13"
##site_code = "NC_UNHC"; start_date = "2016-07-12"; end_date = "2017-08-30"
##site_code = "NC_Mud"; start_date = "2016-07-12"; end_date = "2017-08-30"
streampulse_data = request_data(sitecode=site_code,
    startdate=start_date, enddate=end_date, variables=NULL,
    flags=TRUE, token=NULL)
# head(streampulse_data$data)
# dim(streampulse_data$data)
source('~/git/streampulse/model/sp_functions.R')
source('~/git/streampulse/model/gapfill_functions.R')
fitdata = prep_metabolism(d=streampulse_data, type=model_type,
    model=model_name, interval=interval,
    rm_flagged='none', fillgaps=fillgaps)
    # rm_flagged=list('Bad Data', 'Questionable'), fillgaps=fillgaps)

plot(fitdata$DO.sat, type='l')
plot(fitdata$DO.sat, type='l', xlim=c(3625,3645), xaxs='i')
fitdata$DO.sat[3622:3645] = mean(fitdata$DO.sat)

plot(fitdata$DO.obs, type='l')
plot(fitdata$DO.obs, type='l', xlim=c((3798-191),(3897-186)), xaxs='i')
fitdata$DO.obs[(3798-191):(3897-186)] = 11

fitdata$depth[fitdata$depth <= 0] = 0.00001
plot(fitdata$depth, type='l')

modelfit = fit_metabolism(fitdata)

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
modelfit = readRDS(paste0('~/Desktop/untracked/sm_out/fit_NC_Eno_2016-07-11_2017-08-30_bayes_binned_obsproc_trapezoid_DO-mod_stan.rds'))

predictions = predict_metabolism(modelfit)

#check daiy k-er correlation
daily_er = modelfit@fit$daily$ER_daily_mean
daily_k = modelfit@fit$daily$K600_daily_mean
plot(daily_k, daily_er)
daily_k[is.na(daily_k)] = mean(daily_k, na.rm=TRUE)
daily_er[is.na(daily_er)] = mean(daily_er, na.rm=TRUE)
cor(daily_k, daily_er)

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
    file=paste0('~/Desktop/untracked/sm_figs/plots_',
    site_code, '_', start_date, '_', end_date, '.pdf'), compress=FALSE)
diag_plots(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')
dev.off()

# average_plot(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')
# cumulative_plots(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')
# kernel_plot(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')
