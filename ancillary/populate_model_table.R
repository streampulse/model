library(RMariaDB)
library(DBI)
# library(stringr)

setwd('/home/mike/git/streampulse/server_copy/sp/shiny/data/')
pw = readLines('/home/mike/Dropbox/stuff_2/credentials/spdb.txt')

con = dbConnect(MariaDB(), dbname='sp', username='root', password=pw)

mods = list.files(pattern='modOut')

time_now = Sys.time()
attr(time_now, 'tzone') = 'UTC'

model_deets = data.frame()
for(i in 1:length(mods)){
    spec_str = substr(mods[i], 8, nchar(mods[i])-4)
    spec_vec = strsplit(spec_str, '_')[[1]]

    preds = readRDS(paste0('predictions_', spec_str, '.rds'))
    mod_out = readRDS(mods[i])

    modyear = spec_vec[3]
    # modyear = ifelse(substr(spec_vec[3],1,4) == substr(spec_vec[4],1,4),
    #     substr(spec_vec[3],1,4), 0)

    rmse = sqrt(mean((mod_out$data$DO.mod - mod_out$data$DO.obs)^2, na.rm=TRUE))

    gpp_upperCI = abs(preds$GPP.upper - preds$GPP)
    gpp_lowerCI = abs(preds$GPP.lower - preds$GPP)
    gpp_95ci = mean(c(gpp_upperCI, gpp_lowerCI), na.rm=TRUE)
    er_upperCI = abs(preds$ER.upper - preds$ER)
    er_lowerCI = abs(preds$ER.lower - preds$ER)
    er_95ci = mean(c(er_upperCI, er_lowerCI), na.rm=TRUE)

    prop_pos_er = sum(preds$ER > 0, na.rm=TRUE) / length(preds$ER)
    prop_neg_gpp = sum(preds$GPP < 0, na.rm=TRUE) / length(preds$GPP)

    pearson = cor(mod_out$fit$daily$ER_mean, mod_out$fit$daily$K600_daily_mean,
        use='na.or.complete')

    # coverage = as.numeric(as.Date(spec_vec[4]) - as.Date(spec_vec[3]))
    coverage = as.numeric(as.Date(preds$date[nrow(preds)]) -
            as.Date(preds$date[1]))

    rc = ifelse(spec_vec[1] == 'NC' && spec_vec[2] != 'Eno', TRUE, FALSE)

    kmax = max(mod_out$fit$daily$K600_daily_mean, na.rm=TRUE)

    model_deets = rbind(model_deets, data.frame(region=spec_vec[1],
        site=spec_vec[2], start_date=as.Date(preds$date[1]),
        end_date=as.Date(preds$date[nrow(preds)]),
        requested_variables='all',
        year=modyear, run_finished=time_now, model='streamMetabolizer',
        method='bayes', engine='stan',
        rm_flagged="Bad Data,Questionable", used_rating_curve=rc,
        pool='binned', proc_err=TRUE, obs_err=TRUE, proc_acor=FALSE,
        ode_method='trapezoid', deficit_src='DO-mod', interv='15 min',
        fillgaps='interpolation', estimate_areal_depth=TRUE, O2_GOF=rmse,
        GPP_95CI=gpp_95ci, ER_95CI=er_95ci, prop_pos_ER=prop_pos_er,
        prop_neg_GPP=prop_neg_gpp, ER_K600_cor=pearson, coverage=coverage,
        kmax=kmax))
}

dbWriteTable(con, 'model', model_deets, append=TRUE)

dbDisconnect(con)
