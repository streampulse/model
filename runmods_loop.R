# setup ####

rm(list=ls()); cat('\014')
# library(devtools) #devtools is necessary for installing packages from GitHub.
# install_github('streampulse/StreamPULSE', dependencies=TRUE)
# install.packages('streamMetabolizer', dependencies=TRUE,
                 # repos=c('https://owi.usgs.gov/R','https://cran.rstudio.com'))

library(StreamPULSE)
library(streamMetabolizer)
library(ks)
library(RColorBrewer)

setwd('C:/Users/vlahm/Desktop/untracked')

#metaboplots funcs
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
season_ts_func = function (ts_full, gpp_var, er_var){
    avg_trajectory <- aggregate(ts_full, by = list(ts_full[,
        "DOY"]), FUN = mean, na.rm = TRUE)
    sd_trajectory <- aggregate(ts_full, by = list(ts_full[,
        "DOY"]), FUN = sd, na.rm = TRUE)
    llim <- min(c(avg_trajectory[, er_var],
        avg_trajectory[, gpp_var]), na.rm = T)
    ulim <- max(c(avg_trajectory[, er_var],
        avg_trajectory[, gpp_var]), na.rm = T)
    # plot(avg_trajectory[, gpp_var],
    plot(avg_trajectory$DOY, avg_trajectory[, gpp_var],
        type = "l", col = "red", xlab = '', ylab = expression(paste("gO"[2] *
                " m"^"-2" * " d"^"-1")), ylim = c(llim, ulim), lwd = 2,
        xaxt='n', xlim=c(1, 366))
    month_labs = month.abb
    month_labs[seq(2, 12, 2)] = ''
    axis(1, seq(1, 365, length.out=12), month_labs)
    # t = avg_trajectory$Date
    # yearstarts = match(unique(substr(t,1,4)), substr(t,1,4))
    # monthstarts = match(unique(substr(t,1,7)), substr(t,1,7))
    # axis(1, yearstarts, rep('', length(yearstarts)), lwd.ticks=2, tck=-0.05)
    # axis(1, yearstarts, substr(t[yearstarts],1,4), line=1, tick=FALSE)
    # month_abbs = month.abb[as.numeric(substr(t[monthstarts],6,7))]
    # axis(1, monthstarts[-1], month_abbs[-1])
    lines(avg_trajectory[, "DOY"], avg_trajectory[, er_var],
        # lines(avg_trajectory[, er_var],
        col = "steelblue", lwd = 2)
    lines(avg_trajectory[, "DOY"], avg_trajectory[, "NPP"],
        # lines(avg_trajectory$NPP,
        col = "grey60", lwd = 2)
    abline(h = 0)
    legend("topleft", inset = c(0.1, -0.15), ncol = 3, xpd = TRUE,
        c("GPP", "NEP", "ER"), bty = "n", lty = c(1, 1, 1),
        lwd = 2, col = c("red", "grey60", "steelblue"))
}
cumulative_func = function (ts_full, gpp_var, er_var){
    na_rm <- na.omit(ts_full[, c("Year", "DOY", gpp_var, er_var,
        "NPP")])
    na_rm$csum_gpp <- ave(na_rm[, gpp_var], na_rm[, "Year"],
        FUN = cumsum)
    na_rm$csum_er <- ave(na_rm[, er_var], na_rm[, "Year"], FUN = cumsum)
    na_rm$csum_npp <- ave(na_rm[, "NPP"], na_rm[, "Year"], FUN = cumsum)
    lim <- max(abs(na_rm[, c("csum_gpp", "csum_er")]))
    pal <- rev(brewer.pal(7, "Spectral"))
    cols <- setNames(data.frame(unique(na_rm[, "Year"]), pal[1:length(unique(na_rm[,
        "Year"]))]), c("Year", "color"))
    csum_merge <- merge(na_rm, cols, by = "Year", type = "left")
    plot(csum_merge$DOY, csum_merge[, "csum_gpp"], pch = 20,
        cex = 1.5, col = paste(csum_merge[, "color"]), ylim = c(0,
            lim), xaxt = "n", xlim=c(0, 366),  ylab = "Cumulative GPP",
        type='p')
    legend("topleft", paste(c(cols[, "Year"])), lwd = c(1, 1),
        col = paste(cols[, "color"]), cex = 0.9)
    plot(csum_merge$DOY, csum_merge[, "csum_er"], pch = 20,
        cex = 1.5, col = paste(csum_merge[, "color"]), ylim = c(-lim,
            0), xaxt = "n", xlim=c(0, 366), ylab = "Cumulative ER",
        type='p')
    plot(csum_merge$DOY, csum_merge[, "csum_npp"], pch = 20,
        cex = 1.5, col = paste(csum_merge[, "color"]), ylim = c(-lim,
            lim), ylab = "Cumulative NPP", xlim=c(0, 366),
        xlab = '', type='p', xaxt='n')
    month_labs = month.abb
    month_labs[seq(2, 12, 2)] = ''
    axis(1, seq(1, 365, length.out=12), month_labs)
    abline(h = 0, col = "grey60", lty = 2)
}
kernel_func = function (ts_full, gpp_var, er_var){
    kernel <- kde(na.omit(ts_full[, c(gpp_var, er_var)]))
    k_lims <- max(abs(c(min(ts_full[, er_var], na.rm = TRUE),
        max(ts_full[, gpp_var], na.rm = TRUE))))
    plot(kernel, xlab = expression(paste("GPP (gO"[2] * " m"^"-2" *
            " d"^"-1" * ")")), ylab = expression(paste("ER (gO"[2] *
                    " m"^"-2" * " d"^"-1" * ")")), ylim = c(-k_lims, 0),
        xlim = c(0, k_lims), display = "filled.contour2", col = c(NA,
            "gray80", "gray60", "gray40"))
    abline(0, -1)
    legend("topright", c("75%", "50%", "25%"), bty = "n", lty = c(1,
        1, 1), lwd = 2, col = c("gray80", "gray60", "gray40"))
}
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

# choose sites and dates ####
site_deets = data.frame(
    site_code=c('FL_SF700','FL_SF2500','FL_SF2800','FL_NR1000','FL_ICHE2700',
                'FL_WS1500','PR_QS'),
    start_date = c('2017-01-01','2017-01-01','2017-01-01','2017-01-01',
                   '2017-01-01','2017-01-01','2015-01-01'),
    end_date = c('2017-12-31','2017-12-31','2017-12-31','2017-12-31',
                 '2017-12-31','2017-12-31','2015-11-31'),
    stringsAsFactors=FALSE)

site_deets = site_deets[1,]
site_deets[,3] = '2017-02-01'
site_deets[,1] = 'NC_Eno'

# run ####
results = matrix('', ncol=5, nrow=nrow(site_deets))
colnames(results) = c('Region', 'Site', 'Result', 'Kmax', 'K_ER_cor')

for(i in 1:nrow(site_deets)){

    site_code = site_deets$site_code[i]
    start_date = site_deets$start_date[i]
    end_date = site_deets$end_date[i]

    results[i,1] = strsplit(site_code, '_')[[1]][1]
    results[i,2] = strsplit(site_code, '_')[[1]][2]

    #request
    streampulse_data = try(request_data(sitecode=site_code,startdate=start_date,
                                    enddate=end_date))

    if(class(streampulse_data) == 'try-error'){
        results[i,3] = 'request error'
        next
    }

    #prep
    fitdata = try(prep_metabolism(d=streampulse_data,type='bayes',
        model='streamMetabolizer'))

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
    #predictions = readRDS('sm_out/predictions_NC_Eno_2017-01-01_2017-02-01_bayes_binned_obsproc_trapezoid_DO-mod_stan.rds')

    #plot output
    pdf(width=7, height=7,
        file=paste0('sm_figs/output2_',
            site_code, '_', start_date, '_', end_date, '.pdf'), compress=FALSE)
    diag_plots(predictions[,c('date','GPP','ER')], 'date', 'GPP', 'ER')
    dev.off()

    results[i,3] = 'Run Finished'

}

write.csv(results, 'results.csv', row.names=FALSE)
