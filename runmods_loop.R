# setup ####

rm(list=ls()); cat('\014')
# install.packages('devtools')
# library(devtools) #devtools is necessary for installing packages from GitHub.
# install_github('streampulse/StreamPULSE', dependencies=TRUE)
# install.packages('streamMetabolizer', dependencies=TRUE,
#                  repos=c('https://owi.usgs.gov/R','https://cran.rstudio.com'))
# install.packages('ks')
# install.packages('RColorBrewer')

library(StreamPULSE)
library(streamMetabolizer)
library(ks)
library(RColorBrewer)

#expects sm_figs and sm_out directories at this location
# setwd('C:/Users/vlahm/Desktop/untracked')
# setwd('D:/untracked')
setwd('~/Desktop/untracked')

#metaboplots funcs
processing_func = function (ts) {
    gpp = ts$GPP; gppup = ts$GPP.upper; gpplo = ts$GPP.lower
    er = ts$ER; erup = ts$ER.upper; erlo = ts$ER.lower
    # ts[!is.na(gpp) & gpp < 0 | !is.na(gpp) & gpp > 100, 'GPP'] = NA
    # ts[!is.na(er) & er > 0, 'ER'] = NA
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
season_ts_func = function (ts_full, suppress_NEP=FALSE){

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
        y=c(gpplo, rev(gppup)), col=adjustcolor('red', alpha.f=0.3),
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
        y=c(erlo, rev(erup)), col=adjustcolor('steelblue', alpha.f=0.3),
        border=NA)
    abline(h=0)
    if(suppress_NEP){
        # plot(1,1, col=adjustcolor('red',alpha.f=0.2))
        legend("topleft", inset=c(0, -0.1), ncol=2, xpd=TRUE,
            legend=c("GPP", "ER"), bty="n", lty=1,
            lwd=2, col=c("red", "steelblue"),
            x.intersp=c(.5,.5))#, text.width=.05)
        legend('topright', inset=c(0.1, -0.1), ncol=2, xpd=TRUE,
            bty="n", lty=1, legend=c('','95CI'),
            col=c(adjustcolor('red', alpha.f=0.3),
                adjustcolor('steelblue', alpha.f=0.3)),
            x.intersp=c(-.1,.5), text.width=.05, lwd=3)

    } else {
        lines(avg_trajectory$DOY, avg_trajectory$NPP, col="darkorchid3", lwd=2)
        # plot(1,1, col=adjustcolor('red',alpha.f=0.2))
        legend("topleft", ncol=3, xpd=TRUE,
        # legend("topleft", ncol=5, xpd=TRUE,
            c("GPP", "NEP", "ER"), bty="n", lty=1,
            # c("GPP", "NEP", "ER", '', '95CI'), bty="n", lty=1,
            lwd=2, col=c("red", "darkorchid3", "steelblue"),
                # adjustcolor('red', alpha.f=0.3),
                # adjustcolor('steelblue', alpha.f=0.3)),
            inset=c(0, -0.1))
        # x.intersp=c(.5,.5,.5,.3,1.3), text.width=.05)
    }
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
diag_plots = function (ts, main, suppress_NEP=FALSE){
    ts_full = processing_func(ts)
    layout(matrix(c(1, 1, 3, 1, 1, 4, 2, 2, 5), 3, 3, byrow=TRUE),
        widths=c(1, 1, 2))
    par(cex=0.6)
    par(mar=c(3, 4, 0.1, 0.1), oma=c(3, 0.5, 0.5, 0.5))
    par(tcl=-0.25)
    par(mgp=c(2, 0.6, 0))
    kernel_func(ts_full, main)
    par(mar=c(0, 4, 2, 0.1))
    season_ts_func(ts_full, suppress_NEP)
    par(mar=c(0, 4, 0, 0.1))
    cumulative_func(ts_full)
}

# choose sites and dates ####
# site_deets = data.frame(
#     site_code = c('SE_M18','SE_M6'),#c('PR_QS', 'PR_Icacos'),
#     start_date = c('2016-06-08', '2016-06-08'),#c('2017-03-01', '2017-03-01'),
#     end_date = c('2016-10-12', '2016-10-12'),#c('2017-12-31', '2017-12-31'),
#     stringsAsFactors=FALSE)

# write.csv(site_deets, 'site_deets.csv', row.names=FALSE)
site_deets = read.csv('site_deets.csv', stringsAsFactors=FALSE)
# site_deets = site_deets[substr(site_deets$site_code, 1, 2) != 'SE',]
# site_deets = site_deets[site_deets$site_code != 'MD_DRKR',]
site_deets = site_deets[site_deets$skip != 'x',]
site_deets = site_deets[site_deets$site_code == 'MD_DRKR',]
site_deets = site_deets[site_deets$site_code == 'FL_SF700',]
site_deets = site_deets[site_deets$site_code == 'PR_QS',]
md_drkr
md_gfvn
pr_qs
dim(site_deets)
site_deets = site_deets[29:57,]
site_deets = site_deets[49:53,]#NC minus Eno
site_deets = site_deets[48:54,]#NC minus Eno
a = readRDS('~/git/streampulse/server_copy/sp/shiny/data/modOut_NC_Eno_2017-01-01_2017-12-31_bayes_binned_obsproc_trapezoid_DO-mod_stan.rds')
b = readRDS('~/git/streampulse/server_copy/sp/shiny/data/predictions_NC_Eno_2017-01-01_2017-12-31_bayes_binned_obsproc_trapezoid_DO-mod_stan.rds')
a$data_daily$date[1]
a$data_daily$date[length(a$data_daily$date)]
site_deets = site_deets[1,]

# run ####
results = matrix('', ncol=5, nrow=nrow(site_deets))
colnames(results) = c('Region', 'Site', 'Result', 'Kmax', 'K_ER_cor')

zq = read.csv('~/Dropbox/streampulse/data/rating_curves/ZQ_data.csv')
offsets = read.csv('~/Dropbox/streampulse/data/rating_curves/sensor_offsets.csv')

for(i in 1:nrow(site_deets)){

    site_code = site_deets$site_code[i]
    token = site_deets$tok[i]
    start_date = site_deets$start_date[i]
    end_date = site_deets$end_date[i]

    results[i,1] = strsplit(site_code, '_')[[1]][1]
    results[i,2] = strsplit(site_code, '_')[[1]][2]

    #request
    streampulse_data = try(request_data(sitecode=site_code,startdate=start_date,
                                    enddate=end_date, token=token))

    if(class(streampulse_data) == 'try-error'){
        results[i,3] = 'request error'
        next
    }

    #get rating curve data if needed (comment if not)
    site = strsplit(site_code, '_')[[1]][2]
    Z = zq[zq$site == site, 'level_m']
    Q = zq[zq$site == site, 'discharge_cms']
    offset = offsets[offsets$site == site, 2] / 100

    #prep
    fitdata = try(prep_metabolism(d=streampulse_data, type='bayes',
        model='streamMetabolizer',
        zq_curve=list(sensor_height=offset, Z=Z, Q=Q, fit='power',
            ignore_oob_Z=TRUE, plot=TRUE)))

    if(class(fitdata) == 'try-error'){
        results[i,3] = 'prep error'
        next
    }

    #plot input
    plotvars = colnames(fitdata)[! colnames(fitdata) %in% c('solar.time')]
    pdf(width=5, height=9,
        file=paste0('sm_figs/input_',
                    site_code, '_', start_date, '_', end_date, '.pdf'), compress=FALSE)
    par(mfrow=c(length(plotvars),1), mar=c(0,0,0,0), oma=c(4,4,.5,.5))
    t = as.Date(fitdata$solar.time)
    for(j in plotvars){
        plot(fitdata[,j], type='l', xlab='', xaxt='n', xaxs='i', las=2)
        mtext(j, 2, 2.5)
        if(j == plotvars[length(plotvars)]){
            yearstarts = match(unique(substr(t,1,4)), substr(t,1,4))
            monthstarts = match(unique(substr(t,1,7)), substr(t,1,7))
            axis(1, yearstarts, substr(t[yearstarts],1,4), line=1, tick=FALSE)
            axis(1, monthstarts, substr(t[monthstarts],6,7))
        }
    }
    dev.off()

    #fit
    modelfit = try(fit_metabolism(fitdata))
    # modelfit = readRDS('sm_out/fit_NC_Eno_2017-01-01_2017-02-01_bayes_binned_obsproc_trapezoid_DO-mod_stan.rds')

    if(class(modelfit) == 'try-error'){
        results[i,3] = 'fit error'
        next
    }

    #save fit object
    saveRDS(modelfit, paste('sm_out/fit',
        site_code, start_date, end_date,
        'bayes_binned_obsproc_trapezoid_DO-mod_stan.rds', sep='_'))

    #get diagnostic stats
    results[i,4] = as.character(round(max(modelfit@fit$daily$K600_daily_mean,
        na.rm=TRUE), 2))

    daily_er = modelfit@fit$daily$ER_daily_mean
    daily_k = modelfit@fit$daily$K600_daily_mean
    daily_k[is.na(daily_k)] = mean(daily_k, na.rm=TRUE)
    daily_er[is.na(daily_er)] = mean(daily_er, na.rm=TRUE)
    results[i,5] = as.character(round(cor(daily_k, daily_er), 2))

    #predict
    predictions = predict_metabolism(modelfit)

    #save prediction object
    saveRDS(predictions, paste('sm_out/predictions',
        site_code, start_date, end_date,
        'bayes_binned_obsproc_trapezoid_DO-mod_stan.rds', sep='_'))
    # predictions = readRDS('~/Desktop/untracked/PR/predictions_PR_Icacos_2017-03-01_2017-12-31_bayes_binned_obsproc_trapezoid_DO-mod_stan.rds')

    #plot output
    pdf(width=7, height=7,
        file=paste0('sm_figs/output_',
            site_code, '_', start_date, '_', end_date, '.pdf'), compress=FALSE)
    diag_plots(predictions, site_code, suppress_NEP=TRUE)
    dev.off()

    results[i,3] = 'Run Finished'

}

write.csv(results, 'results1.csv', row.names=FALSE)
