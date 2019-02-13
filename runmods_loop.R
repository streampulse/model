# setup ####

rm(list=ls()); cat('\014')
# install.packages('devtools')
# library(devtools)
# install_github('streampulse/StreamPULSE', ref='master', dependencies=TRUE)
# install.packages('streamMetabolizer', dependencies=TRUE,
#                  repos=c('https://owi.usgs.gov/R','https://cran.rstudio.com'))
# install.packages('ks')
# install.packages('RColorBrewer')

library(StreamPULSE)
library(streamMetabolizer)
# library(ks)
# library(RColorBrewer)

#expects sm_figs and sm_out directories at this location
# setwd('C:/Users/vlahm/Desktop/untracked')
setwd('K:/untracked')
# setwd('~/Desktop/untracked')

# choose sites and dates ####
# site_deets = data.frame(
#     site_code = c('SE_M18','SE_M6'),#c('PR_QS', 'PR_Icacos'),
#     start_date = c('2016-06-08', '2016-06-08'),#c('2017-03-01', '2017-03-01'),
#     end_date = c('2016-10-12', '2016-10-12'),#c('2017-12-31', '2017-12-31'),
#     stringsAsFactors=FALSE)

# write.csv(site_deets, 'site_deets.csv', row.names=FALSE)

#filter ####
site_deets = read.csv('site_deets.csv', stringsAsFactors=FALSE)
site_deets$int = '15 min'
site_deets$int[site_deets$site_code == 'RI_CorkBrk'] = '30 min'
site_deets$int[(site_deets$site_code == 'FL_NR1000' &
                  site_deets$start_date == '2016-01-01')] = '60 min'

# site_deets = site_deets[substr(site_deets$site_code, 1, 2) != 'SE',]
site_deets$skip[15:nrow(site_deets)] = ''
site_deets = site_deets[site_deets$skip != 'x',]
site_deets = site_deets[1:31,]
site_deets = site_deets[32:62,]

# run ####
results = matrix('', ncol=4, nrow=nrow(site_deets))
colnames(results) = c('Region', 'Site', 'Year', 'Result')#, 'Kmax', 'K_ER_cor')

zq = read.csv('C:/Users/vlahm/Desktop/model/ZQ_data.csv')
offsets = read.csv('C:/Users/vlahm/Desktop/model/sensor_offsets.csv')

for(i in 1:nrow(site_deets)){

    site_code = site_deets$site_code[i]
    token = site_deets$tok[i]
    start_date = site_deets$start_date[i]
    end_date = site_deets$end_date[i]

    results[i,1] = strsplit(site_code, '_')[[1]][1]
    results[i,2] = strsplit(site_code, '_')[[1]][2]
    results[i,3] = substr(start_date, 1, 4)

    #request
    streampulse_data = try(request_data(sitecode=site_code,startdate=start_date,
                                    enddate=end_date, token=token))

    if(class(streampulse_data) == 'try-error'){
        results[i,3] = 'request error'
        next
    }

    if(site_code %in% c('NC_UEno','NC_Stony','NC_NHC','NC_Mud','NC_UNHC')){

        #use rating curve data for some sites
        site = strsplit(site_code, '_')[[1]][2]
        Z = zq[zq$site == site, 'level_m']
        Q = zq[zq$site == site, 'discharge_cms']
        offset = offsets[offsets$site == site, 2] / 100

        fitdata = try(prep_metabolism(d=streampulse_data, type='bayes',
            model='streamMetabolizer', interval=site_deets$int[i],
            zq_curve=list(sensor_height=offset, Z=Z, Q=Q, fit='power',
                ignore_oob_Z=TRUE, plot=TRUE)))
    } else {
        fitdata = try(prep_metabolism(d=streampulse_data, type='bayes',
            model='streamMetabolizer', interval=site_deets$int[i]))
    }

    if(class(fitdata) == 'try-error'){
        results[i,3] = 'prep error'
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

    results[i,3] = 'Run Finished'

}

write.csv(results, 'results1.csv', row.names=FALSE)
