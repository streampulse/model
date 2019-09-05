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
    d = read.csv(paste0('modout/', site, '/', fn),
        colClasses=c('date'='Date'))
    rbind(d, matrix(NA, nrow=365 - nrow(d), ncol=ncol(d),
        dimnames=list(c()) #here
    daily['GPP', site, ] = padded_d
    seasonal['GPP', site, ] = tapply(d$GPP_mean, cut(d$date, #then here
}
d = read.csv('modout/NC_NHC/NC_NHC_2017-01-01_2017-12-31_daily.csv',
    stringsAsFactors=FALSE)
d$GPP_mean
