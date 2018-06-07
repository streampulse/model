library(streamMetabolizer)
library(StreamPULSE)

setwd('~/Desktop/untracked/single_year_fits/')
fits = list.files(pattern='fit_')
for(f in fits){
    print(f)
    filename_tail = substr(f, 5, nchar(f))
    fit = readRDS(f)
    # predictions = predict_metabolism(fit)
    data_daily = fit@data_daily
    data = fit@data
    fit = fit@fit
    mod_out = list('data_daily'=data_daily, 'data'=data, 'fit'=fit)
    saveRDS(mod_out,
        paste0('~/git/streampulse/server_copy/sp/shiny/data/modOut_',
            filename_tail))
    # saveRDS(predictions,
    #     paste0('~/git/streampulse/server_copy/sp/shiny/data/predictions_',
    #         filename_tail))
}
