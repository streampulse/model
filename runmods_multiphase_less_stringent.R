# rm(list=ls()); cat('\014')

# install.packages('devtools')
# library(devtools)
# install_github('streampulse/StreamPULSE', ref='master', dependencies=TRUE)
# install.packages('streamMetabolizer', dependencies=TRUE,
#                  repos=c('https://owi.usgs.gov/R','https://cran.rstudio.com'))

library(StreamPULSE)
library(streamMetabolizer)
library(plyr)

setwd('~/Desktop/untracked')

site_deets = read.csv('~/git/streampulse/model/site_deets2.csv',
    stringsAsFactors=FALSE)
site_deets = site_deets[1:4,]
# site_deets = site_deets[5:8,]
# site_deets = site_deets[9:12,]
# site_deets = site_deets[13:15,]
# site_deets = site_deets[16:18,]
# site_deets = site_deets[substr(site_deets$site_code, 1, 2) == 'RI',]
# run_id = 2
# run_set = 1:15 + 15 * (run_id - 1)
# site_deets = site_deets[run_set,]

results = matrix('', ncol=4, nrow=nrow(site_deets))
colnames(results) = c('Region', 'Site', 'Year', 'Result')#, 'Kmax', 'K_ER_cor')

zq = read.csv('~/git/streampulse/model/ZQ_data.csv')
offsets = read.csv('~/git/streampulse/model/sensor_offsets.csv')

for(i in 1:nrow(site_deets)){

    site_code = site_deets$site_code[i]
    token = site_deets$tok[i]
    start_date = site_deets$start_date[i]
    end_date = site_deets$end_date[i]
    int = site_deets$int[i]
    write(paste(i, site_code, substr(start_date, 1, 4), Sys.time()),
        '~/Desktop/untracked/model_run_log.txt', append=TRUE)

    # #establish K600_lnQ_nodes_meanlog prior
    # av_veloc = site_deets$velocity_ms[i]
    # av_slope = site_deets$slope_prop[i]
    # av_depth = site_deets$depth_m[i]
    # av_disch = site_deets$discharge_m3s[i]
    # logK = log(5.62) + 0.504*log(av_veloc * av_slope) -
    #     0.575*log(av_depth) - 0.0892*log(av_disch)

    results[i,1] = strsplit(site_code, '_')[[1]][1]
    results[i,2] = strsplit(site_code, '_')[[1]][2]
    results[i,3] = substr(start_date, 1, 4)

    #request
    if(site_code == 'NC_ColeMill'){
        streampulse_data = try(request_data(sitecode=site_code, startdate=start_date,
            enddate=end_date, token=token,
            variables=c('DO_mgL','DOsat_pct','satDO_mgL','WaterPres_kPa',
                'Depth_m','WaterTemp_C','Light_PAR','AirPres_kPa')))
    } else if(substr(site_code, 1, 2) == 'WI'){
        streampulse_data = try(request_data(sitecode=site_code, startdate=start_date,
            enddate=end_date, token=token,
            variables=c('DO_mgL','DOsat_pct','satDO_mgL','WaterPres_kPa',
                'Level_m','WaterTemp_C','Light_PAR','AirPres_kPa','Discharge_m3s')))
    } else {
        streampulse_data = try(request_data(sitecode=site_code,startdate=start_date,
            enddate=end_date, token=token))
    }

    if(class(streampulse_data) == 'try-error'){
        results[i,4] = 'request error'
        next
    }

    if(site_code %in% c('NC_UEno','NC_Stony','NC_NHC','NC_Mud','NC_UNHC')){

        #use rating curve data for some sites
        site = strsplit(site_code, '_')[[1]][2]
        Z = zq[zq$site == site, 'level_m']
        Q = zq[zq$site == site, 'discharge_cms']
        offset = offsets[offsets$site == site, 2] / 100

        if(site_code == 'NC_Stony'){
            fitdata = try(prep_metabolism(d=streampulse_data, type='bayes',
                model='streamMetabolizer', interval=int,
                zq_curve=list(sensor_height=offset, Z=Z, Q=Q, fit='power',
                    ignore_oob_Z=TRUE, plot=TRUE), estimate_areal_depth=TRUE,
                estimate_PAR=TRUE, retrieve_air_pres=TRUE))
        } else {
            fitdata = try(prep_metabolism(d=streampulse_data, type='bayes',
                model='streamMetabolizer', interval=int,
                zq_curve=list(sensor_height=offset, Z=Z, Q=Q, fit='power',
                    ignore_oob_Z=TRUE, plot=TRUE), estimate_areal_depth=FALSE,
                estimate_PAR=TRUE, retrieve_air_pres=TRUE))
        }
    } else if(substr(site_code, 1, 2) == 'WI' &&
            ! 'Level_m' %in% streampulse_data$data$variable){
        fitdata = try(prep_metabolism(d=streampulse_data, type='bayes',
            model='streamMetabolizer', interval=int,
            estimate_areal_depth=TRUE, estimate_PAR=TRUE,
            retrieve_air_pres=TRUE))
    } else {
        fitdata = try(prep_metabolism(d=streampulse_data, type='bayes',
            model='streamMetabolizer', interval=int,
            estimate_areal_depth=FALSE, estimate_PAR=TRUE,
            retrieve_air_pres=TRUE))
    }

    if(class(fitdata) == 'try-error'){
        results[i,4] = 'prep error'
        next
    }

    #remove storm days
    df = fitdata$data
    df$Date = as.Date(df$solar.time)
    df_list = split(df, df$Date)
    df_diff = ldply(lapply(df_list,
        function(x){
            range_diff = max(x$discharge) - min(x$discharge)
            return(range_diff)
        }), data.frame)
    colnames(df_diff) = c("Date", "Diff")
    df_diff = df_diff[which(df_diff$Diff <
            (0.5 * median(df$discharge, na.rm=TRUE))),]
    df = df[which(as.character(df$Date) %in% df_diff$Date),]

    #exclude low DO observations, NA rows, impossible water temperatures
    df = df[which(df$DO.obs > 0),]
    df = na.omit(df)
    fitdata$data$temp.water[fitdata$data$temp.water > 40] = NA
    df$Date = NULL

    if(! nrow(df)){
        results[i,4] = 'no modelable days'
        next
    }

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
    dirs = list.dirs('sp20190812', recursive=FALSE)
    write_dir = paste('sp20190812', site_code, sep='/')
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
    file=paste0('~/Desktop/untracked/sp20190812/results', run_id, '.csv'))
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
