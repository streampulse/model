library(StreamPULSE)
library(dplyr)

#prepare for debug gapfill
site_code = 'AZ_OC'; start_date = '2018-06-01'; end_date = '2018-07-01'
sp_data = request_data(sitecode=site_code,
    startdate=start_date, enddate=end_date)
sp_data_prepped = prep_metabolism(d=sp_data, type='bayes',
    model='streamMetabolizer')
sp_data_prepped$data$DO.obs[10:20] = NA
sp_data_prepped$data$DO.obs[100:200] = NA
sp_data_prepped$data[1000:1200,] = NA

df=data.frame(sp_data_prepped$data)
colnames(df) = c('DateTime_UTC', 'DO_mgL', 'DOsat_pct', 'Depth_m', 'WaterTemp_C',
    'Light_lux', 'Discharge_m3s')
df = df[!is.na(df$DateTime_UTC),]
df = mutate(df, DateTime_UTC=as.character(DateTime_UTC))

maxspan_days=5
knn=3

sort_records = function(pl){

    df = pl$df
    df = mutate(df, DateTime_UTC=as.POSIXct(DateTime_UTC, tz='UTC'))
    df = df[order(df$DateTime_UTC),]

    pl = list(df=df, messages=list())
    return(pl)
}

homogenize_intervals = function(pl){

    messages = pl$messages
    df = pl$df

    int_runs = rle(diff(as.numeric(df$DateTime_UTC)))
    uniq_ints = na.omit(unique(int_runs$values))

    if(length(uniq_ints) > 1){
        messages = append(messages,
            paste('Multiple sample intervals detected:',
                paste(uniq_ints / 60, collapse=', '), 'mins.'))
    }

    # get min interval that accounts for >1% of values
    int_count = tapply(int_runs$lengths, int_runs$values, sum)
    int_prop = int_count / sum(int_count)
    min_int = min(as.numeric(names(int_prop[int_prop > 0.01]))) / 60
    min_int = as.difftime(min_int, unit='mins')

    #conform dataframe to this interval
    min_int_df = data.frame(DateTime_UTC=seq.POSIXt(df$DateTime_UTC[1],
        df$DateTime_UTC[nrow(df)], by=min_int))
    df = right_join(df, min_int_df, by='DateTime_UTC')

    pl = list(df=df, messages=messages)
    return(pl)
}

s = pl$df$DO

lin_interp_series = function(s, max_fill_samp){

    gaps_bool = is.na(s)
    na_runs = rle(gaps_bool)
    na_runs = data.frame(lengths=na_runs$lengths, values=na_runs$values,
        ends=cumsum(na_runs$lengths), starts=c(1, ends[-length(ends)] + 1))
    long_na_runs = filter(na_runs, values == TRUE, lengths > max_fill_samp)
    long_na_inds = unlist(mapply(seq, long_na_runs$starts, long_na_runs$ends))

    #linearly interpolate all gaps
    s[gaps_bool] = approx(s, xout=which(gaps_bool))$y

    #replace NAs where gaps are too big to safely interpolate
    s[long_na_inds] = NA

    return(s)
}

linear_interp = function(pl, max_fill_hrs=3){

    messages = pl$messages
    df = pl$df

    samp_int = as.numeric(difftime(df$DateTime_UTC[2],  df$DateTime_UTC[1],
        units='hours'))

    skip_ind = which(colnames(df) %in% c('DateTime_UTC', 'Depth_m',
        'Discharge_m3s'))
    if('Depth_m' %in% colnames(df)){
        df$Depth_m = lin_interp_series(df$Depth_m, max_fill_samp=
    df[, -dt_ind] = apply(df_no_dt[, -dt_ind], 2, lin_interp_series,
        max_fill_samp=max_fill_hrs / samp_int)

    ina <- is.na(x)
    csum <- cumsum(!ina)
    # which gaps are too long
    x[ina] <- approx(x, xout = which(ina))$y
    if(!is.infinite(tol)){ # if not infinite tol, remove threshold values
        wg <- a.n(names(which(table(csum) > tol)))
        x[which(csum%in%wg)[-1]] <- NA
    }
    x

}

# xf = df
pl = list(df=df, messages=list())
pl = sort_records(pl)
pl = homogenize_intervals(pl)
pl = linear_interp(pl)
