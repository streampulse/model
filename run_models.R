#StreamPULSE model script
#Mike Vlah: vlahm13@gmail.com
#last update 3/14/2018

#if you're using RStudio, you can collapse sections with
#Alt+o (Windows/Linux) or Cmd+Opt+o (Mac)

rm(list=ls()); cat('\014') #clear environment and console

#install/update and load streamMetabolizer and StreamPULSE packages ####

#install latest version if you don't have streamMetabolizer already
install.packages('streamMetabolizer', dependencies=TRUE,
    repos=c('https://owi.usgs.gov/R','https://cran.rstudio.com'))

#update streamMetabolizer if you aready have it
update.packages(oldPkgs=c('streamMetabolizer','unitted'), dependencies=TRUE,
    repos=c('https://owi.usgs.gov/R', 'https://cran.rstudio.com'))

#install StreamPULSE package
library(devtools) #devtools is necessary for installing packages from GitHub.
install_github('streampulse/StreamPULSE', dependencies=TRUE)

#load packages
library(streamMetabolizer)
library(StreamPULSE)

#model setup (see documentation for details) ####

model_type = 'bayes'
model_name = 'streamMetabolizer'

#choose site and dates (see data.streampulse.org/sitelist for site codes).
#see data.streampulse.org/visualize for available dates.
site_code = 'SE_M6'; start_date = '2016-06-08'; end_date = '2016-10-12'
site_code = 'SE_M18'; start_date = '2016-06-08'; end_date = '2016-10-12'

#retrieve data from StreamPULSE database
streampulse_data = request_data(sitecode=site_code,
    startdate=start_date, enddate=end_date)

#acquire additional data if necessary; format for streamMetabolizer
fitdata = prep_metabolism(d=streampulse_data, type=model_type,
    model=model_name)

#model ####

#see what you're dealing with
plotvars = colnames(fitdata)[! colnames(fitdata) %in% c('solar.time')]
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

#fit model with stock parameters (uncomment next block for more options)
modelfit = fit_metabolism(fitdata)

#fit model with custom parameters (mm_name, specs, metab are from streamMetabolizer)
# class(fitdata) = 'data.frame' #just a formality, execute and disregard
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

#plot ####

#load additional packages for plotting
library(ks)
library(RColorBrewer)

#source MetaboPlots functions (package coming soon?)
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

#plot
diag_plots(predictions, site_code, suppress_NEP=TRUE)
