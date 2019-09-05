library(stringr)

setwd('~/git/streampulse/model/ancillary/K_sensitivity_analysis/')
sites = dir('modout')

annual = array(NA, dim=c(3, length(sites), 1),
    dimnames=list(c('GPP', 'ER', 'K600'), sites, NULL))
seasonal = array(dim=c(3, length(sites), 4),
    dimnames=list(c('GPP', 'ER', 'K600'), sites, NULL))
daily = array(NA, dim=c(3, length(sites), 365),
    dimnames=list(c('GPP', 'ER', 'K600'), sites, NULL))

for(i in 1:length(sites)){
    site = sites[i]
    fn = na.omit(str_match(dir(paste0('modout/', site)), '.*?_daily.csv'))
    fn = na.omit(str_match(dir(paste0('modout/', site)), '.*?_daily.csv'))
    d = read.csv(paste0('modout/', site, '/', fn),
        colClasses=c('date'='Date'))
    padded_d = rbind(d, matrix(NA, nrow=365 - nrow(d), ncol=ncol(d),
        dimnames=list(NULL, colnames(d))))

    err = FALSE
    tryCatch(daily['GPP', site, ] <- padded_d$GPP_mean,
        error=function(e) err <<- TRUE)
    if(err) next

    daily['ER', site, ] = padded_d$ER_mean
    daily['K600', site, ] = padded_d$K600_daily_mean
    seasonal['GPP', site, ] = tapply(padded_d$GPP_mean,
        cut(1:nrow(padded_d), breaks=4), mean, na.rm=TRUE)
    seasonal['ER', site, ] = tapply(padded_d$ER_mean,
        cut(1:nrow(padded_d), breaks=4), mean, na.rm=TRUE)
    seasonal['K600', site, ] = tapply(padded_d$K600_daily_mean,
        cut(1:nrow(padded_d), breaks=4), mean, na.rm=TRUE)
    annual['GPP', site, ] = mean(padded_d$GPP_mean, na.rm=TRUE)
    annual['ER', site, ] = mean(padded_d$ER_mean, na.rm=TRUE)
    annual['K600', site, ] = mean(padded_d$K600_daily_mean, na.rm=TRUE)
}

d = read.csv('modout/NC_NHC/NC_NHC_2017-01-01_2017-12-31_daily.csv',
    stringsAsFactors=FALSE)
d$GPP_mean
