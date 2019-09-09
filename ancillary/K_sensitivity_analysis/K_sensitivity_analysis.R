library(StreamPULSE)
library(streamMetabolizer)
library(plyr)

setwd('~/git/streampulse/model/ancillary/K_sensitivity_analysis/')

site_deets = read.csv('~/git/streampulse/model/site_deets.csv',
    stringsAsFactors=FALSE)
site_deets = site_deets[substr(site_deets$site_code, 1, 2) == 'NC',]
site_deets = site_deets[site_deets$start_date == '2017-01-01', ]
site_deets = site_deets[! site_deets$site_code %in% c('NC_Stony', 'NC_ColeMill',
    'NC_UEno', 'NC_UNHC'),]
site_deets$end_date
site_deets$start_date[site_deets$site_code %in% c('NC_NHC')] = '2017-01-01'
site_deets$end_date[site_deets$site_code %in% c('NC_NHC')] = '2017-12-31'
site_deets$start_date[site_deets$site_code %in% c('NC_Eno')] = '2016-01-01'
site_deets$end_date[site_deets$site_code %in% c('NC_Eno')] = '2016-12-31'

run_id = 1

results = matrix('', ncol=4, nrow=nrow(site_deets))
colnames(results) = c('Region', 'Site', 'Year', 'Result')#, 'Kmax', 'K_ER_cor')

zq = read.csv('~/git/streampulse/model/ZQ_data.csv')
offsets = read.csv('~/git/streampulse/model/sensor_offsets.csv')

# i=1; K=2
for(i in 1:nrow(site_deets)){
    # for(K in c(1, 2, 4, 8, 16, 32, 64, 128)){
    for(K in run_id){

        site_code = site_deets$site_code[i]
        token = site_deets$tok[i]
        start_date = site_deets$start_date[i]
        end_date = site_deets$end_date[i]
        int = site_deets$int[i]
        write(paste(i, site_code, substr(start_date, 1, 4), Sys.time()),
            '~/Desktop/untracked/model_run_log.txt', append=TRUE)

        results[i,1] = strsplit(site_code, '_')[[1]][1]
        results[i,2] = strsplit(site_code, '_')[[1]][2]
        results[i,3] = substr(start_date, 1, 4)

        #request
        if(site_code == 'NC_ColeMill'){
            streampulse_data = try(request_data(sitecode=site_code, startdate=start_date,
                enddate=end_date, token=token,
                variables=c('DO_mgL','DOsat_pct','satDO_mgL','WaterPres_kPa',
                    'Depth_m','WaterTemp_C','Light_PAR','AirPres_kPa')))
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

        #exclude NA rows, impossible water temperatures, omit discharge
        df = fitdata$data
        df = na.omit(df)
        fitdata$data$temp.water[fitdata$data$temp.water > 40] = NA
        df$Date = NULL
        df$discharge = NULL

        if(! nrow(df)){
            results[i,4] = 'no modelable days'
            next
        }

        #set specs
        bayes_name_new = mm_name(type='bayes', pool_K600="none", err_obs_iid=TRUE,
            err_proc_iid = TRUE, ode_method = "trapezoid", deficit_src='DO_mod',
            engine='stan')
        bayes_specs_new = specs(bayes_name_new)
        bayes_specs_new$K600_daily_meanlog = exp(K)
        bayes_specs_new$K600_daily_sdlog = 0.001

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
        dirs = list.dirs('.', recursive=FALSE)
        write_dir = paste('modout', site_code, sep='/')
        if(! site_code %in% dirs){
            dir.create(write_dir, recursive=TRUE)
        }

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
}

write.csv(results, row.names=FALSE, file=paste0('results', run_id, '.csv'))
