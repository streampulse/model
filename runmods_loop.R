# setup ####

# rm(list=ls()); cat('\014')

# install.packages('devtools')
# library(devtools)
# install_github('streampulse/StreamPULSE', ref='master', dependencies=TRUE)
# install.packages('streamMetabolizer', dependencies=TRUE,
#                  repos=c('https://owi.usgs.gov/R','https://cran.rstudio.com'))
# install.packages('ks')
# install.packages('RColorBrewer')

library(tidyverse)
library(StreamPULSE)
library(streamMetabolizer)
library(logging)

dir.create('logs', showWarnings = FALSE)
logging::basicConfig()
logging::addHandler(logging::writeToFile,
                    logger = 'sp',
                    file = 'logs/0_master.log')

# library(ks)
# library(RColorBrewer)

#expects sm_figs and sm_out directories at this location
# setwd('C:/Users/vlahm/Desktop/untracked')
# setwd('K:/untracked')
setwd('~/Desktop/untracked')

site_deets = read.csv('~/git/streampulse/model/site_deets_2020run.csv',
                      stringsAsFactors=FALSE)

Mode <- function(x, na.rm = TRUE){

    if(na.rm){
        x <- na.omit(x)
    }

    ux <- unique(x)
    mode_out <- ux[which.max(tabulate(match(x, ux)))]
    return(mode_out)

}

# generate new sitedata file (comment this before running multiple sessions)####

# site_deets_orig = read.csv('~/git/streampulse/model/site_deets_postAZ.csv',
#     stringsAsFactors=FALSE)
# site_deets_updatedyears = site_deets_orig %>%
#     as_tibble() %>%
#     group_by(site_code) %>%
#     summarize(across(everything(), first),
#               .groups = 'drop') %>%
#     mutate(start_date = '2020-01-01',
#            end_date = '2020-12-31',
#            skip = '',
#            run_status = '')
# site_deets_updatedyears <- mutate(site_deets_updatedyears,
#        start_date = '2019-01-01',
#        end_date = '2019-12-31') %>%
#     bind_rows(site_deets_updatedyears)
#
# site_deets_new = read_csv('~/Downloads/all_basic_site_data.csv') %>%
#     mutate(site_code = paste(regionID, siteID, sep='_')) %>%
#     filter(! is.na(firstRecord) & ! is.na(lastRecord),
#            embargoDaysRemaining == 0,
#            ! site_code %in% site_deets_orig$site_code,
#            ! grepl('(nwis\\-|\\-down|\\-up)', site_code)) %>%
#     select(site_code, firstRecord, lastRecord) %>%
#     mutate(year = purrr::map2(substr(firstRecord, 1, 4),
#                               substr(lastRecord, 1, 4),
#                               seq,
#                               by=1)) %>%
#     unnest(cols = year) %>%
#     mutate(start_date = paste0(year, '-01-01'),
#            end_date = paste0(year, '-12-31'),
#            velocity_ms = NA, depth_m = NA, discharge_m3s = NA, slope_prop = NA,
#            int = '15 min',
#            skip = '',
#            run_status = '',
#            tok = site_deets_orig$tok[1]) %>%
#     select(-firstRecord, -lastRecord, -year)
#
# #for now, also remove Alice's sites
# site_deets_new = site_deets_new[! grepl('NC_', site_deets_new$site_code), ]
#
# site_deets = bind_rows(site_deets_updatedyears, site_deets_new)
# site_deets$int[site_deets$site_code == 'FL_NR1000'] = '60 min'
#
# write_csv(site_deets, '~/git/streampulse/model/site_deets_2020run.csv')

# filter ####

# site_deets = site_deets[site_deets$skip != 'x',]
# site_deets = site_deets[site_deets$run_status == 'r',]
# site_deets = site_deets[substr(site_deets$site_code, 1, 2) == 'NC',]

# run ####
results = matrix('', ncol=4, nrow=nrow(site_deets))
colnames(results) = c('Region', 'Site', 'Year', 'Result')#, 'Kmax', 'K_ER_cor')

zq = read.csv('~/git/streampulse/model/ZQ_data.csv')
offsets = read.csv('~/git/streampulse/model/sensor_offsets.csv')

#KEEP SITE DEETS THE SAME AND RUN 2:4
# for(i in 1:nrow(site_deets)){
runstart = 169
runend = 209
for(i in runstart:runend){

    site_code = site_deets$site_code[i]
    token = site_deets$tok[i]
    start_date = site_deets$start_date[i]
    end_date = site_deets$end_date[i]
    int = site_deets$int[i] #there's a check for this below

    #establish K600_lnQ_nodes_meanlog prior (no parameter to access this within
    #fit metab at the moment
    av_veloc = site_deets$velocity_ms[i]
    av_slope = site_deets$slope_prop[i]
    av_depth = site_deets$depth_m[i]
    av_disch = site_deets$discharge_m3s[i]
    logK = log(5.62) + 0.504*log(av_veloc * av_slope) -
        0.575*log(av_depth) - 0.0892*log(av_disch)

    results[i,1] = strsplit(site_code, '_')[[1]][1]
    results[i,2] = strsplit(site_code, '_')[[1]][2]
    results[i,3] = substr(start_date, 1, 4)

    #request

    if(site_code == 'NC_ColeMill'){
        streampulse_data = try(request_data(sitecode=site_code,startdate=start_date,
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

    intervals = streampulse_data$data %>%
        group_by(variable) %>%
        summarize(interval = Mode(diff(DateTime_UTC)),
                  .groups = 'drop')

    print(intervals)

    interval = intervals %>%
        summarize(interval = Mode(interval)) %>%
        pull(interval)

    if(interval != 15){
        logging::logwarn(paste('check interval:', site_code, start_date, end_date,
                               paste(intervals$interval, collapse = ','),
                               sep = ' | '),
                         logger = 'sp.module')
        int = paste0(interval, 'min')
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

    #plot input
    # plotvars = colnames(fitdata)[! colnames(fitdata) %in% c('solar.time')]
    # pdf(width=5, height=9,
    #     file=paste0('sm_figs/input_',
    #                 site_code, '_', start_date, '_', end_date, '.pdf'), compress=FALSE)
    # par(mfrow=c(length(plotvars),1), mar=c(0,0,0,0), oma=c(4,4,.5,.5))
    # t = as.Date(fitdata$solar.time)
    # for(j in plotvars){
    #     plot(fitdata[,j], type='l', xlab='', xaxt='n', xaxs='i', las=2)
    #     mtext(j, 2, 2.5)
    #     if(j == plotvars[length(plotvars)]){
    #         yearstarts = match(unique(substr(t,1,4)), substr(t,1,4))
    #         monthstarts = match(unique(substr(t,1,7)), substr(t,1,7))
    #         axis(1, yearstarts, substr(t[yearstarts],1,4), line=1, tick=FALSE)
    #         axis(1, monthstarts, substr(t[monthstarts],6,7))
    #     }
    # }
    # dev.off()

    #fit
    fitdata$data$temp.water[fitdata$data$temp.water > 40] = NA
    modelfit = try(fit_metabolism(fitdata, skip_prompts=TRUE))
    # modelfit = readRDS('sm_out/fit_NC_Eno_2017-01-01_2017-02-01_bayes_binned_obsproc_trapezoid_DO-mod_stan.rds')

    if(class(modelfit) == 'try-error'){
        results[i,4] = 'fit error'
        next
    }

    # #save fit object
    # saveRDS(modelfit, paste('sm_out/fit',
    #     site_code, start_date, end_date,
    #     'bayes_binned_obsproc_trapezoid_DO-mod_stan.rds', sep='_'))

    #get diagnostic stats
    # results[i,4] = as.character(round(max(modelfit@fit$daily$K600_daily_mean,
    #     na.rm=TRUE), 2))
    #
    # daily_er = modelfit@fit$daily$ER_daily_mean
    # daily_k = modelfit@fit$daily$K600_daily_mean
    # daily_k[is.na(daily_k)] = mean(daily_k, na.rm=TRUE)
    # daily_er[is.na(daily_er)] = mean(daily_er, na.rm=TRUE)
    # results[i,5] = as.character(round(cor(daily_k, daily_er), 2))

    #predict
    # predictions = predict_metabolism(modelfit)

    #save prediction object
    # saveRDS(predictions, paste('sm_out/predictions',
    #     site_code, start_date, end_date,
    #     'bayes_binned_obsproc_trapezoid_DO-mod_stan.rds', sep='_'))
    # predictions = readRDS('~/Desktop/untracked/PR/predictions_PR_Icacos_2017-03-01_2017-12-31_bayes_binned_obsproc_trapezoid_DO-mod_stan.rds')

    #plot output
    # pdf(width=7, height=7,
    #     file=paste0('sm_figs/output_',
    #         site_code, '_', start_date, '_', end_date, '.pdf'), compress=FALSE)
    # diag_plots(predictions, site_code, suppress_NEP=TRUE)
    # dev.off()

    results[i,4] = 'Run Finished'

}

write.csv(results, paste0('results_', site_deets$site_code[runstart], '-',
                          site_deets$site_code[runend], '.csv'), row.names=FALSE)
