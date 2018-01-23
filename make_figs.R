rm(list=ls()); cat('\014')
source('~/git/streampulse/model/gapfill_functions.R')
source('~/git/streampulse/model/sp_functions.R')

model_type = "mle"
model_name = "streamMetabolizer"
site_code = "AZ_LV"
start_date = "2017-07-07"
end_date = "2017-12-25"
streampulse_data = request_data(sitecode=site_code,
    startdate=start_date, enddate=end_date, variables=NULL,
    flags=FALSE, token=NULL)
head(streampulse_data)
source('~/git/streampulse/model/gapfill_functions.R')
fitdata = prep_metabolism(d=streampulse_data, type=model_type,
    model=model_name, interval='15 min', fillgaps='interpolation')

# modelfit = fit_metabolism(fitdata)

class(fitdata) = "data.frame"
engine = 'stan'; pool_K600='binned'; proc_err = TRUE
modname <- mm_name(type='bayes', pool_K600='binned',
    err_obs_iid=TRUE, err_proc_acor=FALSE, err_proc_iid=TRUE,
    ode_method = 'trapezoid', deficit_src='DO_mod', engine='stan')
modspecs <- specs(modname)
# modspecs$K600_lnQ_nodes_centers <- seq(from=md$logQ_min[n],
#     to = md$logQ_max[n], by = md$Qnodes_steps[n])
modspecs$K600_lnQ_nodes_centers = seq(from=.74, to = 5.25, by = .75)
modelfit <- metab(specs = modspecs, data = fitdata)

predictions = predict_metabolism(modelfit)

#plot
ymin = min(c(predictions$ER,predictions$GPP), na.rm=TRUE)
ymax = max(c(predictions$ER,predictions$GPP), na.rm=TRUE)
ind = which(!is.na(predictions$GPP))
plot(predictions$GPP[ind], type='l', ylim=c(ymin, ymax),
    main='streamMetabolizer', xlab='Date', xaxt='n', ylab='g O2/m^2/d', las=1)
axis(1, seq(1,92,10), predictions$date[ind][seq(1,92,10)])
lines(predictions$ER[ind], col='red')
legend('topright', legend=c('GPP','ER'), lty=1, col=c('black','red'), bty='n')
