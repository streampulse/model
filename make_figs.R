rm(list=ls()); cat('\014') #clear environment and console

#setup ####
# library(devtools) #devtools is necessary for installing packages from GitHub.
# install_github('streampulse/StreamPULSE', dependencies=TRUE)
# install.packages('streamMetabolizer', dependencies=TRUE,
# repos=c('https://owi.usgs.gov/R','https://cran.rstudio.com'))

library(StreamPULSE)
library(streamMetabolizer)
library(ks)
library(RColorBrewer)
processing_func = function (ts) {
    gpp = ts$GPP; gppup = ts$GPP.upper; gpplo = ts$GPP.lower
    er = ts$ER; erup = ts$ER.upper; erlo = ts$ER.lower
    ts[!is.na(gpp) & gpp < 0 | !is.na(gpp) & gpp > 100, 'GPP'] = NA
    ts[!is.na(er) & er > 0, 'ER'] = NA
    # if('tbl' %in% class(ts$date)){

    # full_dates = as.data.frame(ts$date)
    # colnames(full_dates) = 'Date'

    # colnames(ts)[which(colnames(ts) == 'date')] = 'Date'
    ts_full = as.data.frame(ts)
    # } else {
    #     ts_date = ts$date
    #     ts$Date = format(ts_date, "%Y-%m-%d")
    #     full_dates = setNames(data.frame(
    #         seq(from = as.Date(paste(min(format(ts$date, "%Y")),
    #             "-01-01", sep = "")),
    #             to = as.Date(paste(max(format(ts[,
    #                 date_var], "%Y")), "-12-31", sep = "")), by = 1)), "Date")
    # }

    # ts_full = merge(full_dates, ts, by = c("Date"), all = TRUE)
    ts_full$Year = as.numeric(format(ts_full$date, "%Y"))
    ts_full$DOY = as.numeric(format(ts_full$date, "%j"))
    ts_full$NPP = ts_full$GPP + ts_full$ER
    return(ts_full)
}
season_ts_func = function (ts_full){

    ts_full = ts_full[-c(1,nrow(ts_full)),-c(1,8,9,10,11)]

    avg_trajectory = aggregate(ts_full, by=list(ts_full$DOY),
        FUN=mean, na.rm=TRUE)

    gpp = avg_trajectory$GPP
    gppup = avg_trajectory$GPP.upper; gpplo = avg_trajectory$GPP.lower
    er = avg_trajectory$ER
    erup = avg_trajectory$ER.upper; erlo = avg_trajectory$ER.lower
    doy = avg_trajectory$DOY

    sd_trajectory = aggregate(ts_full, by=list(ts_full$DOY),
        FUN=sd, na.rm=TRUE)

    llim = min(c(gpplo, erlo), na.rm=TRUE)
    ulim = max(c(gppup, erup), na.rm=TRUE)
    # plot(avg_trajectory[, gpp_var],
    plot(doy, avg_trajectory$GPP,
        type="l", col="red", xlab='', ylab=expression(paste("gO"[2] *
                " m"^"-2" * " d"^"-1")), ylim=c(llim, ulim), lwd=2,
        xaxt='n', xlim=c(1, 366))
    polygon(x=c(doy, rev(doy)),
        y=c(gpplo, rev(gppup)), col=adjustcolor('red',alpha.f=0.15),
        border=NA)
    # t = avg_trajectory$Date
    # yearstarts = match(unique(substr(t,1,4)), substr(t,1,4))
    # monthstarts = match(unique(substr(t,1,7)), substr(t,1,7))
    # axis(1, yearstarts, rep('', length(yearstarts)), lwd.ticks=2, tck=-0.05)
    # axis(1, yearstarts, substr(t[yearstarts],1,4), line=1, tick=FALSE)
    # month_abbs = month.abb[as.numeric(substr(t[monthstarts],6,7))]
    # axis(1, monthstarts[-1], month_abbs[-1])
    lines(doy, avg_trajectory$ER, col="steelblue", lwd=2)
    polygon(x=c(doy, rev(doy)),
        y=c(erlo, rev(erup)), col=adjustcolor('steelblue',alpha.f=0.15),
        border=NA)
    lines(avg_trajectory$DOY, avg_trajectory$NPP, col="darkorchid3", lwd=2)
    abline(h=0)
    legend("topleft", inset=c(0.1, -0.15), ncol=3, xpd=TRUE,
        c("GPP", "NEP", "ER"), bty="n", lty=c(1, 1, 1),
        lwd=2, col=c("red", "darkorchid3", "steelblue"))
    month_labs = month.abb
    month_labs[seq(2, 12, 2)] = ''
    axis(1, seq(1, 365, length.out=12), month_labs)
}
cumulative_func = function (ts_full){
    ts_full = ts_full[-c(1,nrow(ts_full)),-c(1,8,9,10)]
    na_rm = na.omit(ts_full)
    na_rm$csum_gpp = ave(na_rm$GPP, na_rm$Year, FUN=cumsum)
    na_rm$csum_er = ave(na_rm$ER, na_rm$Year, FUN=cumsum)
    na_rm$csum_npp = ave(na_rm$NPP, na_rm$Year, FUN=cumsum)
    lim = max(abs(na_rm[, c("csum_gpp", "csum_er")]))
    pal = rev(brewer.pal(7, "Spectral"))
    cols = setNames(data.frame(unique(na_rm$Year), pal[1:length(unique(na_rm[,
        "Year"]))]), c("Year", "color"))
    csum_merge = merge(na_rm, cols, by="Year", type="left")
    plot(csum_merge$DOY, csum_merge$csum_gpp, pch=20,
        cex=1.5, col=paste(csum_merge$color), type='p',
        ylim=c(0, lim), xaxt="n", xlim=c(0, 366), ylab="Cumulative GPP")
    legend("topleft", paste(c(cols$Year)), lwd=c(1, 1),
        col=paste(cols$color), cex=0.9)
    plot(csum_merge$DOY, csum_merge$csum_er, pch=20,
        cex=1.5, col=paste(csum_merge$color), type='p',
        ylim=c(-lim, 0), xaxt="n", xlim=c(0, 366), ylab="Cumulative ER")
    plot(csum_merge$DOY, csum_merge$csum_npp, pch=20,
        cex=1.5, col=paste(csum_merge$color), ylim=c(-lim, lim),
        ylab="Cumulative NPP", xlim=c(0, 366),
        xlab='', type='p', xaxt='n')
    month_labs = month.abb
    month_labs[seq(2, 12, 2)] = ''
    axis(1, seq(1, 365, length.out=12), month_labs)
    abline(h=0, col="grey60", lty=2)
}
kernel_func = function (ts_full, main){
    ts_full = ts_full[-c(1,nrow(ts_full)),-c(1,8,9,10)]

    kernel = kde(na.omit(ts_full[, c('GPP', 'ER')]))
    k_lims = max(abs(c(min(ts_full$ER, na.rm=TRUE),
        max(ts_full$GPP, na.rm=TRUE))))
    plot(kernel, xlab=expression(paste("GPP (gO"[2] * " m"^"-2" *
            " d"^"-1" * ")")),
        ylab=expression(paste("ER (gO"[2] * " m"^"-2" * " d"^"-1" * ")")),
        ylim=c(-k_lims, 0), xlim=c(0, k_lims), display="filled.contour2",
        col=c(NA, "gray80", "gray60", "gray40"))
    mtext(main, 3, line=-2)
    abline(0, -1)
    legend("topright", c("75%", "50%", "25%"), bty="n", lty=c(1,
        1, 1), lwd=2, col=c("gray80", "gray60", "gray40"))
}
# ts=predictions[,c('date','GPP','ER','GPP.lower','GPP.upper','ER.lower',
    # 'ER.upper')]; date_var='date'; gpp_var='GPP'
# er_var='ER'; main='PR_Icacos'
diag_plots = function (ts, main){
    ts_full = processing_func(ts)
    layout(matrix(c(1, 1, 3, 1, 1, 4, 2, 2, 5), 3, 3, byrow=TRUE),
        widths=c(1, 1, 2))
    par(cex=0.6)
    par(mar=c(3, 4, 0.1, 0.1), oma=c(3, 0.5, 0.5, 0.5))
    par(tcl=-0.25)
    par(mgp=c(2, 0.6, 0))
    kernel_func(ts_full, main)
    par(mar=c(0, 4, 2, 0.1))
    season_ts_func(ts_full)
    par(mar=c(0, 4, 0, 0.1))
    cumulative_func(ts_full)
}
# ts_full = processing_func(predictions)
# diag_plots(predictions, 'PR_Icacos')


#choices ####

model_type = "bayes"
# model_type = "mle"
model_name = "streamMetabolizer"
fillgaps='interpolation'
interval='15 min'
#done, ##problem
site_code = "SE_AbiskoM1"; start_date = "2016-06-28"; end_date = '2016-10-25'
site_code = "SE_AbiskoM2"; start_date = "2016-06-14"; end_date = '2016-09-22'
site_code = "SE_AbiskoM6"; start_date = "2016-06-10"; end_date = '2016-10-23'
site_code = "SE_AbiskoM9"; start_date = "2016-06-15"; end_date = '2016-09-08'
site_code = "SE_AbiskoM10"; start_date = "2016-06-15"; end_date = '2016-08-16'
site_code = "SE_AbiskoM16"; start_date = "2016-05-27"; end_date = '2016-10-30'
site_code = "SE_AbiskoM17"; start_date = "2016-06-15"; end_date = '2016-09-08'
site_code = "SE_M18"; start_date = "2016-06-08"; end_date = '2016-10-22'
site_code = "SE_M6"; start_date = "2016-06-08"; end_date = '2016-10-22'
site_code = "PR_QS"; start_date = "2015-01-01"; end_date = '2015-11-30'#ed "2017-12-20"
# site_code = "PR_RioIcacosTrip"
# site_code = "PR_Prieta"
site_code = "RI_CorkBrk"; start_date = "2016-01-01"; end_date = '2016-12-31'#ed "2017-01-03", sd=2014-06-23
#site_code = "CT_BUNN"; start_date = "2015-05-20"; end_date = '2016-05-20'#ed "2016-11-15"
#site_code = "CT_HUBB"; start_date = "2015-05-20"; end_date = '2016-05-20'#ed "2016-11-11"
#site_code = "CT_FARM"; start_date = "2015-05-20"; end_date = '2016-05-20'#ed "2016-11-11"
#site_code = "CT_Unio"; start_date = "2015-05-20"; end_date = '2016-05-20'#ed "2016-11-11"
#site_code = "VT_Pass"; start_date = "2015-06-26"; end_date = '2016-06-25' #ed 2016-10-27
#site_code = "VT_POPE"; start_date = "2015-06-03"; end_date = '2016-06-02' #ed 2016-11-11
#site_code = "VT_MOOS"; start_date = "2015-06-03"; end_date = '2016-11-11' #ed 2016-11-11
##site_code = "PR_Icacos"; start_date = "2016-06-09"; end_date = '2016-12-15'
#site_code = "CT_STIL"; start_date = "2015-05-20"; end_date = '2016-05-20'#ed "2016-11-15"
#site_code = "MD_GFVN"; start_date = "2016-02-18"; end_date = '2016-11-19'
#site_code = "MD_GFCP"; start_date = "2016-05-07"; end_date = '2016-11-19'
#site_code = "MD_GFGB"; start_date = "2016-02-19"; end_date = '2016-11-11'
#site_code = "MD_DRKR"; start_date = "2016-02-18"; end_date = '2016-11-19'
# site_code = "MD_POBR"; start_date = "2016-02-19"; end_date = '2016-11-19'
site_code = "MD_BARN"; start_date = "2016-02-19"; end_date = '2016-11-19'
#site_code = "VT_SLPR"; start_date = "2015-06-03"; end_date = '2016-06-03' #ed 2016-11-11
# site_code = "AZ_LV"; start_date = "2017-08-07"; end_date = "2017-12-25" #sd 2017-07-07 but no o2
#site_code = "AZ_OC"; start_date = "2016-11-15"; end_date = "2017-12-03" #sd 11-13
#site_code = "AZ_SC"; start_date = "2017-02-08"; end_date = "2017-03-28"
# site_code = "AZ_WB"; start_date = "2017-08-04"; end_date = "2017-12-27"
site_code = "NC_Eno"; start_date = "2017-01-01"; end_date = "2017-08-30" #2016-07-11;
# site_code = "NC_UEno"; start_date = "2017-01-01"; end_date = "2017-12-31" #2016-07-12; 2018-02-09
# site_code = "NC_Stony"; start_date = "2016-12-01"; end_date = "2017-06-01" #2016-06-30; 2017-08-09
# site_code = "NC_NHC"; start_date = "2016-09-14"; end_date = "2017-09-13"
site_code = "NC_Mud"; start_date = "2017-01-01"; end_date = "2017-12-31" #2016-07-12; 2018-01-25
site_code = "NC_Eno"; start_date = "2017-01-01"; end_date = "2017-10-20"#2016-07-11
# site_code = "NC_UEno"; start_date = "2017-01-01"; end_date = "2017-12-31" #2016-07-12; 2018-02-09
# site_code = "NC_Stony"; start_date = "2016-12-01"; end_date = "2017-06-01" #2016-06-30; 2017-08-09
# site_code = "NC_NHC"; start_date = "2017-01-01"; end_date = "2018-01-01" #2016-09-14; 2018-01-25
# site_code = "NC_UNHC"; start_date = "2017-01-01"; end_date = "2017-12-11" #2016-07-12
# site_code = 'WI_BEC'; start_date = '2017-01-26'; end_date = '2018-01-25' #2009-10-02; 2018-01-25
# site_code = 'WI_BRW'; start_date = '2017-01-27'; end_date = '2018-01-26' #2014-06-13; 2018-01-26
site_code = 'FL'; start_date = '2017-01-27'; end_date = '2018-01-26' #2014-06-13; 2018-01-26


#run ####
# source('~/git/streampulse/model/sp_functions.R')
# source('~/git/streampulse/model/gapfill_functions.R')

streampulse_data = request_data(sitecode=site_code,
    startdate=start_date, enddate=end_date, variables=NULL,
    flags=TRUE)
    # flags=TRUE, token='67f2d1a026b6c9e3446e') #gerard
    # flags=TRUE, token='7e4cd63a38de5d4a4715') #maria
    flags=TRUE, token='cfb849bbcbe2aa5859d0') #miguel
head(streampulse_data$data)
unique(streampulse_data$data$variable)
sum(streampulse_data$data$flagtype != '')
nrow(streampulse_data$data)
# streampulse_data$data[streampulse_data$data$variable=='Light_PAR',
#     c('variable','flagtype')]
# z = streampulse_data$data
# z = z$DateTime_UTC[z$variable == 'WaterTemp_C']
# z = z$DateTime_UTC[z$variable == 'DOsat_pct']
# unique(substr(z, 15, 16))

plotvars = unique(streampulse_data$data$variable)
par(mfrow=c(length(plotvars),1), mar=c(0,0,0,0), oma=c(4,4,.5,.5))
t = as.Date(sort(unique(streampulse_data$data$DateTime_UTC)))
for(i in plotvars){
    pvar = streampulse_data$data$value[which(streampulse_data$data$variable == i)]
    plot(pvar,
        type='l', xlab='', xaxt='n', xaxs='i', las=2)
    lin = ifelse(which(plotvars == i) %% 2 == 0, 2.5, 3)
    mtext(i, 2, cex=0.6, line=lin)
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


zq = read.csv('~/Dropbox/streampulse/data/rating_curves/ZQ_data.csv')
site = strsplit(site_code, '_')[[1]][2]
Z = zq[zq$site == site, 'level_m']
Q = zq[zq$site == site, 'discharge_cms']
# dim(streampulse_data$data)
# source('~/git/streampulse/model/sp_functions.R')
# source('~/git/streampulse/model/gapfill_functions.R')

# source('~/git/streampulse/model/sp_functions.R')
fitdata = prep_metabolism(d=streampulse_data, type='bayes',
    model='streamMetabolizer', interval='15 min',
    rm_flagged=list('Bad Data', 'Questionable'), fillgaps=fillgaps,
    # zq_curve=list(sensor_height=NULL, fit='power',
    # zq_curve=list(sensor_height=NULL, a=.316, b=9.529, fit='power',
        # plot=TRUE),
    # zq_curve=list(sensor_height=NULL, Z=NULL, Q=NULL, a=2.5662, b=2.7534, #nhc
    # zq_curve=list(sensor_height=NULL, Z=NULL, Q=NULL, a=1.7257, b=3.0586, #ueno
    # zq_curve=list(sensor_height=NULL, Z=Z, Q=Q, a=NULL, b=NULL,
        # fit='linear', plot=TRUE),
        # fit='power', plot=TRUE),
    estimate_areal_depth=TRUE)

# apply(fitdata[,-1], 2, function(x) sum(is.infinite(x)))

# plot(fitdata$DO.sat, type='l')
# plot(fitdata$DO.sat, type='l', xlim=c(3625,3645), xaxs='i')
# fitdata$DO.sat[3622:3645] = mean(fitdata$DO.sat)
#
# plot(fitdata$DO.obs, type='l')
# plot(fitdata$DO.obs, type='l', xlim=c((3798-191),(3897-186)), xaxs='i')
# fitdata$DO.obs[(3798-191):(3897-186)] = 11


plotvars = colnames(fitdata)[! colnames(fitdata) %in% c('solar.time')]
pdf(width=5, height=9,
    file=paste0('~/Desktop/untracked/sm_figs/input_',
    site_code, '_', start_date, '_', end_date, '.pdf'), compress=FALSE)
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
dev.off()

# source('~/git/streampulse/model/sp_functions.R')
modelfit = fit_metabolism(fitdata)


modelfit = readRDS(paste0('~/Desktop/untracked/mod_fit_objects/fit_FL_WS1500_2017-01-01_2017-12-31_bayes_binned_obsproc_trapezoid_DO-mod_stan.rds'))
max(modelfit@fit$daily$K600_daily_mean, na.rm=TRUE) #shoudnt be over 90
plot(density(modelfit@fit$daily$K600_daily_mean, na.rm=TRUE))
modelfit@fit$daily$K600_daily_mean
saveRDS(modelfit@fit, '~/Desktop/untracked/WS1500_fit.rds')
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

#check daiy k-er correlation
par(mfrow=c(1,1))
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
# predictions = readRDS('~/Desktop/untracked/region_zips/PR_streamMetabolizer_output/predictions/predictions_PR_Icacos_2016-10-01_2016-12-01_bayes_binned_obsproc_trapezoid_DO-mod_stan.rds')

# predictions = readRDS(paste0('~/Desktop/untracked/streampulse_model_runs_',
#     'round1/mod_objects/predictions_',
#     site_code, '_', start_date, '_', end_date, '_bayes_binned_obsproc_trapezoid',
#     '_DO-mod_stan.rds'))

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
# p = processing_func(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')

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
# ts_full=p; gpp_var='GPP'; er_var='ER'

pdf(width=7, height=7,
    file=paste0('~/Desktop/untracked/sm_figs/output2_',
    site_code, '_', start_date, '_', end_date, '.pdf'), compress=FALSE)
diag_plots(predictions[,c('date','GPP','ER')], 'PR_Icacos')
dev.off()

# average_plot(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')
# cumulative_plots(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')
# kernel_plot(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')

#replot all ####
fdir = '~/Desktop/untracked/streampulse_model_runs_round1/mod_objects/'
fdir = '~/Desktop/untracked/sm_out/'
f = list.files(fdir)
for(i in f){
    fullpath = paste0(fdir, i)
    predictions = readRDS(fullpath)
    # p = processing_func(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')

    yrs = substr(predictions$date,1,4)
    years_represented = unique(yrs)
    if(length(years_represented > 1)){
        mode_year = names(which.max(table(yrs)))
        predictions = predictions[substr(predictions$date, 1, 4) == mode_year,]
    }

    splt = strsplit(i, '_')[[1]]
    site_code = paste0(splt[2:3], collapse='_')
    start_date = splt[4]
    end_date = splt[5]
    pdf(width=6, height=5.5,
        file=paste0('~/Desktop/untracked/all_round1_figs/output_',
        site_code, '_', start_date, '_', end_date, '.pdf'), compress=FALSE)
    diag_plots(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER',
        main=site_code)
    dev.off()
}
