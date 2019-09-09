library(stringr)

setwd('~/git/streampulse/model/ancillary/K_sensitivity_analysis/')
sites = dir('modout')

K_vals = c(2, 8, 16, 32, 64)
# K_vals = c(1, 8, 16, 32, 64)

#setup ####

annual = array(NA, dim=c(3, length(sites), length(K_vals), 1),
    dimnames=list(c('GPP', 'ER', 'K600'), sites, K_vals, NULL))
seasonal = array(NA, dim=c(3, length(sites), length(K_vals), 4),
    dimnames=list(c('GPP', 'ER', 'K600'), sites, K_vals, NULL))
daily = array(NA, dim=c(3, length(sites), length(K_vals), 365),
    dimnames=list(c('GPP', 'ER', 'K600'), sites, K_vals, NULL))

# annual = array(rnorm(3*3*5*1, 0, 1), dim=c(3, length(sites), length(K_vals), 1),
#     dimnames=list(c('GPP', 'ER', 'K600'), sites, K_vals, NULL))
# seasonal = array(rnorm(3*3*5*4, 0, 4), dim=c(3, length(sites), length(K_vals), 4),
#     dimnames=list(c('GPP', 'ER', 'K600'), sites, K_vals, NULL))
# daily = array(rnorm(3*3*5*365, 0, 365), dim=c(3, length(sites), length(K_vals), 365),
#     dimnames=list(c('GPP', 'ER', 'K600'), sites, K_vals, NULL))

for(i in 1:length(sites)){
    site = sites[i]
    fns = na.omit(str_match(dir(paste0('modout/', site)), '.*?_daily.csv'))
    for(fn in fns){
        d = read.csv(paste0('modout/', site, '/', fn),
            colClasses=c('date'='Date'))
        padded_d = rbind(d, matrix(NA, nrow=365 - nrow(d), ncol=ncol(d),
            dimnames=list(NULL, colnames(d))))
        fixed_K = str_match(fn, '([0-9]+)_daily.csv$')[,2]

        err = FALSE
        tryCatch(daily['GPP', site, fixed_K, ] <- padded_d$GPP_mean,
            error=function(e) err <<- TRUE)
        if(err) next

        daily['ER', site, fixed_K, ] = padded_d$ER_mean
        daily['K600', site, fixed_K, ] = padded_d$K600_daily_mean
        seasonal['GPP', site, fixed_K, ] = tapply(padded_d$GPP_mean,
            cut(1:nrow(padded_d), breaks=4), mean, na.rm=TRUE)
        seasonal['ER', site, fixed_K, ] = tapply(padded_d$ER_mean,
            cut(1:nrow(padded_d), breaks=4), mean, na.rm=TRUE)
        seasonal['K600', site, fixed_K, ] = tapply(padded_d$K600_daily_mean,
            cut(1:nrow(padded_d), breaks=4), mean, na.rm=TRUE)
        annual['GPP', site, fixed_K, ] = mean(padded_d$GPP_mean, na.rm=TRUE)
        annual['ER', site, fixed_K, ] = mean(padded_d$ER_mean, na.rm=TRUE)
        annual['K600', site, fixed_K, ] = mean(padded_d$K600_daily_mean, na.rm=TRUE)
    }
}

#annual plot ####
png(width=7, height=6, units='in', filename='figs/annual.png', type='cairo',
    res=300)
par(mfrow=c(1,1), mar=c(4, 4, 3, 4), oma=c(1,1,1,1))
plot(annual['GPP', sites[1], , ], xaxt='n', xlab='K600', ylab='gO2m^-2d^-1',
    ylim=c(range(annual[c('GPP', 'ER'), , , ])), type='l', col='red',
    main='Annual averages')
axis(1, 1:5, K_vals)
lines(annual['ER', sites[1], , ], col='blue')
for(i in 2:length(sites)){
    lines(annual['GPP', sites[i], , ], col='red', lty=i)
    lines(annual['ER', sites[i], , ], col='blue', lty=i)
}
legend(x='topright', legend=c('GPP', 'ER'), col=c('red', 'blue'),
    inset=c(-0.15,0), xpd=TRUE, lty=1, bty='n', cex=0.6)
legend(x='bottomright', legend=sapply(strsplit(sites, '_'), function(x) x[2]),
    col='black', inset=c(-0.15,0), xpd=TRUE, lty=1:length(sites), bty='n',
    cex=0.6)
dev.off()

#quarterly plots ####
png(width=7, height=7, units='in', filename='figs/quarterly.png', type='cairo',
    res=300)
par(mfrow=c(2,2), mar=c(3, 3, 3, 3), oma=c(3, 3, 3, 3))
for(s in 1:4){
    plot(seasonal['GPP', sites[1], , s], xaxt='n', xlab='K600', ylab='gO2m^-2d^-1',
        ylim=c(range(seasonal[c('GPP', 'ER'), , , ], na.rm=TRUE)), type='l',
        col='red', main=paste('Quarter', s))
    axis(1, 1:5, K_vals)
    lines(seasonal['ER', sites[1], , s], col='blue')
    for(i in 2:length(sites)){
        lines(seasonal['GPP', sites[i], , s], col='red', lty=i)
        lines(seasonal['ER', sites[i], , s], col='blue', lty=i)
    }
    if(s == 2){
        legend(x='bottomright', legend=c('GPP', 'ER'), col=c('red', 'blue'),
            inset=c(-0.4,0), xpd=NA, lty=1, bty='n', cex=1)
    }
    if(s == 4){
        legend(x='topright', legend=sapply(strsplit(sites, '_'), function(x) x[2]),
            col='black', inset=c(-0.4,0), xpd=NA, lty=1:length(sites), bty='n',
            cex=1)
    }
}
mtext('Quarterly averages', 3, outer=TRUE, font=2, cex=1.5)
mtext('gO2m^-2d^-1', 2, outer=TRUE, font=2, cex=1.5)
mtext('K600', 1, outer=TRUE, font=2, cex=1.5)
dev.off()

#daily dist plot ####
png(width=7, height=6, units='in', filename='figs/daily.png', type='cairo',
    res=300)
par(mfrow=c(3,1), mar=c(1, 4, 3, 5), oma=c(4,2,2,2))
for(s in 1:length(sites)){
    boxplot(t(daily['GPP', sites[i], , ]), border='red', at=seq(0.9, 4.9, 1),
        pars=list(boxcol='transparent', medlty='blank', medpch=16,
            whisklty=c(1,1), medcex=0.7, outcex=0, staplelty='blank'),
        xaxt='n', main=sites[s], bty='l',
        ylim=range(daily[c('GPP', 'ER'), , , ], na.rm=TRUE))
    boxplot(t(daily['ER', sites[i], , ]), border='blue', at=seq(1.1, 5.1, 1),
        pars=list(boxcol='transparent', medlty='blank', medpch=16,
            whisklty=c(1,1), medcex=0.7, outcex=0, staplelty='blank'),
        add=TRUE, xaxt='n', bty='l',
        ylim=range(daily[c('GPP', 'ER'), , , ], na.rm=TRUE))
    if(s == 2){
        legend(x='right', legend=c('GPP', 'ER'), col=c('red', 'blue'),
            inset=c(-0.15,0), xpd=NA, lty=1, bty='n', cex=1)
    }
    if(s == 3){
        axis(1, at=1:5, labels=K_vals)
    }
}
mtext('Daily distributions', 3, outer=TRUE, font=2, cex=1.5)
mtext('gO2m^-2d^-1', 2, outer=TRUE, font=2, cex=1)
mtext('K600', 1, outer=TRUE, font=2, cex=1, line=2)
dev.off()
