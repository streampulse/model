# rm(list=ls()); cat('\014')

# install.packages('devtools')
# library(devtools)
# install_github('streampulse/StreamPULSE', ref='master', dependencies=TRUE)
# install.packages('streamMetabolizer', dependencies=TRUE,
#                  repos=c('https://owi.usgs.gov/R','https://cran.rstudio.com'))

library(StreamPULSE)
library(streamMetabolizer)
library(plyr)
library(dplyr)

setwd('~/Desktop/untracked')
outfolder = 'sp20190808'

site_deets = read.csv('~/git/streampulse/model/site_deets.csv',
    stringsAsFactors=FALSE)
# site_deets = site_deets[substr(site_deets$site_code, 1, 2) == 'RI',]

results = matrix('', ncol=4, nrow=nrow(site_deets))
colnames(results) = c('Region', 'Site', 'Year', 'Result')#, 'Kmax', 'K_ER_cor')

zq = read.csv('~/git/streampulse/model/ZQ_data.csv')
offsets = read.csv('~/git/streampulse/model/sensor_offsets.csv')

good_days = read.csv('~/Downloads/SP_data_RhatK600_sub1.1.csv',
    stringsAsFactors=FALSE)

runs = unique(select(good_days, State_Site, Year))
runs$State_Site = sub(' ', '_', runs$State_Site)
rownames(runs) = 1:nrow(runs)
runs = runs[1:11,]

for(i in 1:nrow(runs)){

    site_code = runs$State_Site[i]
    # token = site_deets$tok[i]
    start_date = paste0(runs$Year[i], '-01-01')
    end_date = paste0(runs$Year[i], '-12-31')

    df = filter(good_days, State_Site == sub('_', ' ', site_code),
            Year == as.numeric(substr(start_date, 1, 4))) %>%
        select(solar.time, DO.obs, DO.sat, depth, temp.water, light, discharge) %>%
        mutate(solar.time=as.POSIXct(solar.time, tz='UTC'))


    int = site_deets$int[i]
    write(paste(i, site_code, substr(start_date, 1, 4), Sys.time()),
        '~/Desktop/untracked/model_run_log.txt', append=TRUE)

    results[i,1] = strsplit(site_code, '_')[[1]][1]
    results[i,2] = strsplit(site_code, '_')[[1]][2]
    results[i,3] = substr(start_date, 1, 4)


    #set specs
    bayes_name_new = mm_name(type='bayes', pool_K600="binned", err_obs_iid=TRUE,
        err_proc_iid = TRUE, ode_method = "trapezoid", deficit_src='DO_mod',
        engine='stan')
    bayes_specs_new = specs(bayes_name_new)

    n_KQ_nodes = 7
    Q_by_day = tapply(log(df$discharge), substr(df$solar.time, 1, 10), mean)
    bayes_specs_new$K600_lnQ_nodes_centers = seq(from=min(Q_by_day, na.rm=TRUE),
        to=max(Q_by_day, na.rm=TRUE), length.out=n_KQ_nodes)

    #fit
    fit_err = FALSE
    tryCatch({
        dat_metab = metab(bayes_specs_new, data=df)
        dat_fit = get_fit(dat_metab)
    }, error=function(e) {fit_err <<- TRUE})
    if(fit_err){
        results[i,4] = 'fit error'
        next
    }

    #save deets
    dirs = list.dirs(outfolder, recursive=FALSE)
    write_dir = paste(outfolder, site_code, sep='/')
    if(! site_code %in% dirs){
        dir.create(write_dir)
    }
    # write_dir = paste('sp20190717', site_code, substr(end_date, 1, 4), sep='/')

    fn_prefix = paste0(write_dir, '/', site_code, '_', start_date, '_',
        end_date, '_')
    write.csv(dat_fit$daily, paste0(fn_prefix, "daily.csv"), row.names=FALSE)
    write.csv(dat_fit$overall, paste0(fn_prefix, "overall.csv"), row.names=FALSE)
    write.csv(dat_fit$KQ_overall, paste0(fn_prefix, "KQ_overall.csv"), row.names=FALSE)
    specs_out = data.frame(unlist(get_specs(dat_metab)))
    write.csv(specs_out, paste0(fn_prefix, 'specs.csv'))
    daily_out = get_data_daily(dat_metab)
    write.csv(daily_out, paste0(fn_prefix, 'datadaily.csv'), row.names=FALSE)
    data_out = get_data(dat_metab)
    write.csv(data_out, paste0(fn_prefix, 'mod_and_obs_DO.csv'), row.names=FALSE)

    results[i,4] = 'Run Finished'

}

write.csv(results, row.names=FALSE,
    file=paste0('~/Desktop/untracked/', outfolder, '/results', run_id, '.csv'))
# #determine "high" DOsat amplitude by which to filter
# z = model_fit
# k_rhats = z$fit@fit$daily[,c('date','K600_daily_Rhat')]
# mo = get_data(z$fit)
# k_rhats = k_rhats[k_rhats$K600_daily_Rhat <= 1.2,]
# k_rhats = na.omit(k_rhats)
# mo = mo[mo$date %in% k_rhats$date,]
# rngs = tapply(mo$DO.sat, mo$date, range)
# dosat = data.frame(matrix(unlist(rngs), ncol=2, byrow=TRUE,
#     dimnames=list(NULL, c('min', 'max'))))
# dosat$range = dosat$max - dosat$min
