
# x = input_data[,8]; tol=12; algorithm='interpolation'
# x=select(daypoints,-date,-time); tol=0; samp=samp; algorithm='mean'
series_impute = function(x, tol, samp, algorithm, ...){
    # records locations of NA runs longer than tol, imputes all gaps,
    # then replaces NAs for long runs.
    # samp is the number of samples per day
    # algorithm and ... are passed to imputeTS::na.seasplit

    na_locations = which(is.na(x))

    if(length(na_locations)){ #if NAs in x, get indices of long runs
        runs = rle2(diff(na_locations), indices=TRUE, return.list=FALSE)
        runs = runs[runs[,'values'] == 1 & runs[,'lengths'] >= tol-1, ,
            drop=FALSE]
        long_na_runs = mapply(seq, runs[,'starts'], runs[,'stops'] + 1,
            SIMPLIFY=FALSE)
        long_na_run_indices = na_locations[unlist(long_na_runs)]
    }

    #impute
    if(length(x) > samp * 2){ #don't leverage periodicity for <2 days
        x = ts(x, start=2, frequency=samp) #add periodicity info
    } else {
        warning(paste('Less than 2 days of data, so interpolating without',
            'periodicity information.'))
    }
    imputed = as.numeric(suppressWarnings(na.seadec(x,
        algorithm=algorithm, ...)))

    #restore NAs where there were long runs
    if(length(na_locations)){
        imputed[long_na_run_indices] = NA
    }

    return(imputed)
}

#for sums of squared differences between standardized values from NA days and
#all other days
sum_sq_diff = function(x, narow){
    sum((narow - x)^2, na.rm=TRUE)
}

# find similar k days
# df=select(daily_averages, -date); k=3; minobs=2
top_k = function(df, k, minobs){
    # df is dataframe of daily values with a row for each date and a colum for each variable
    # k is number of similar days to find
    # minobs is the min. observations required for a day

    df = df[-c(1,nrow(df)),] #first and last obs are artifacts. removing them.
    cc = complete.cases(df)
    narows = which(!cc)
    fullrows = which(cc)
    D = matrix(NA, nrow(df), k) #holds nearest neighbors for each row with NAs

    enough_full_rows = TRUE
    if(length(fullrows) < 30){
        message(paste('Not filling gaps via k nearest neighbors approach;',
            'fewer than 30 days without NAs for comparison.'))
        enough_full_rows = FALSE
    }

    if(length(narows) & enough_full_rows){

        df = apply(df, 2, scale) #standardize variables for comparison of rows

        #for each row with missing data, find k similar rows with none missing
        for (i in narows){

            #skip days with less than minobs. cant say much about similarity
            if(sum(!is.na(df[i,])) < minobs){ next }

            #get sums of squared differences for each full day
            daily_deltas = apply(df[fullrows,], MARGIN=1, FUN=sum_sq_diff,
                narow=df[i,])
            D[i,] = fullrows[order(daily_deltas)[1:k]] #k nearest days
        }
    }

    #restore bogus days at beginning and end of D (for compatibility)
    D = rbind(rep(NA, ncol(D)), D, rep(NA, ncol(D)))
    return(D) #matrix with rows for NA days and columns for knn days
}

# data prep function for fill_gaps
# adds snap points, gets average data
prep_missing = function(df, nearest_neighbors, daily_averages, mm, samp){
    # df is the data frame
    # daily_averages is the data frame of days
    # mm is the missing days
    # nearest_neighbors is the matching neighbors
    #
    ### MISSING DATA
    missing = filter(df, date %in% daily_averages$date[mm])
    # if missing first and last obs, extend timeseries to include neighbor days
    #  unless first/last day
    if(any(!complete.cases(missing)[c(1,nrow(missing))])){
        if(mm[1]!=1) mm = c(mm[1]-1, mm) # if not first day, extend obs range
        if(tail(mm,1)!=nrow(daily_averages)) mm = c(mm, tail(mm,1)+1) # if not last day, extend obs range
        missing = filter(df, date%in%daily_averages$date[mm]) # missing one step interpolation
        # if any data missing still (i.e., first and last day), add avg of first and last obs
        #  this should catch most NAs for filling missing
        if(any(!complete.cases(missing)[c(1,nrow(missing))])){
            missing[c(1,nrow(missing)),-c(1,2)] =
                apply(missing[c(1,nrow(missing)),-c(1,2)], 2, mean, na.rm=T)
        }
    }
    ndays = length(mm)
    ### SIMILAR DATA
    # grab similar days - pairs of missing day and similar day
    ss = data.frame(date=daily_averages$date[mm],
        match=daily_averages$date[t(nearest_neighbors[mm,])])
    ss = ss[complete.cases(ss),]
    similar = left_join(ss, df, by=c("match"="date")) %>%
        select(-match) %>% group_by(date, time) %>% summarize_all(mean) %>%
        ungroup()
    # make sure that the dates in similar and missing line up
    missing = right_join(missing, select(similar,date,time),
        by=c("date","time"))
    ### DAILY SNAP POINTS
    # add snap points at beginning/end each new day to rescale and
    # match the daily trends
    newdaypoints = which(missing$time %in% missing$time[c(1,nrow(missing))])
    daypoints = missing[newdaypoints,]
    if(any(is.na(daypoints))){
        dayfill = series_impute(select(daypoints,-date,-time), tol=0,
            samp=samp, algorithm='mean')
        missing[newdaypoints,] = data.frame(date=daypoints$date,
            time=daypoints$time, dayfill)
    }
    list(missing=select(missing,-date,-time),
        similar=select(similar,-date,-time), index=select(similar,date,time))
}

# df=input_data; lim=0; samp=samples_per_day; g=1
fill_missing = function(df, nearest_neighbors, daily_averages,
    date_index, maxspan_days, samp, lim=0){

    # df is the input data frame, all numeric data
    # nearest_neighbors are the similar days for each day
    # daily_averages is the daily data
    # date_index are the identifiers
    # maxspan_days is the maximum number of days to gap fill
    # lim is the minimum number of days to fill
    #     will not fill gaps that are less than this
    #     used for testing (b/c the test data have pre-existing gaps)

    # the days that need filling in daily_averages
    filld = which(complete.cases(nearest_neighbors))

    if(length(filld)){ #skip the rest if no data are missing

        # groups for blocks of missing data
        group = cumsum(c(TRUE, diff(filld) > 1))

        for(g in unique(group)){

            # grab missing days
            mm = filld[group == g]
            if(length(mm) >= lim && length(mm) <= maxspan_days){
                pp = prep_missing(df, nearest_neighbors, daily_averages, mm,
                    samp)
                dy = (pp$missing - pp$similar)
                dyhat = linear_fill(dy)
                filled = pp$similar + dyhat
                df[which(paste(df$date,df$time) %in%
                        paste(pp$index$date,pp$index$time)),-c(1,2)] = filled
            }
        }
    }
    data.frame(date_index, select(df,-date,-time), stringsAsFactors=FALSE)
}

# df=dd; maxspan_days=5; knn=3; sint=desired_int; algorithm=fillgaps
gap_fill = function(df, maxspan_days=5, knn=3, sint, algorithm, ...){
    # df is data frame, requires one column as POSIXct date time and the
    # other columns as numeric. the order of columns does not matter.
    # int is the interval between samples

    # check if all but one column is numeric
    if( !(length(which(sapply(df, function(x) inherits(x, "numeric")))) ==
            ncol(df)-1) ){
        stop("ERROR: All but one column in df must be numeric.")
    }
    # check if a posix column exists
    wposix = which(sapply(df, function(x) inherits(x, "POSIXct")))
    if( !(length(wposix)==1) ){
        stop("ERROR: Need at least one column in df with POSIXct datetime.")
    }

    # find POSIXct column, that is the one that we need to break into
    # date and time
    dtcol = colnames(df)[wposix]

    # kind of goofy to do this by date and time, but that's because I
    # translated the code from Python
    input_data = df %>% mutate(date=as.Date(df[,dtcol]),
        time=strftime(df[,dtcol], format="%H:%M:%S")) %>%
        select(-one_of(dtcol)) %>% select(date, time, everything())
    date_index = df %>% select(one_of(dtcol)) # index data

    #get daily sampling frequency so imputation can leverage periodicity
    samples_per_day = 24 * 60 / as.double(sint, units='mins')

    #put knn gapfiller here, before in-line gapfiller

    # impute in-line gaps (runs of NAs within a column) in df
    imputed = as.data.frame(sapply(X=input_data[,-(1:2)],
        FUN=series_impute, tol=12, samp=samples_per_day,
        algorithm=algorithm, ...,
        simplify=FALSE)) #gaps >= tol will not be filled
    input_data = data.frame(select(input_data, date, time), imputed)

    # get averages for days with full sample coverage; otherwise NA (via mean())
    # nearly_complete_day = samples_per_day * 0.95 #could also use partial days
    daily_averages = input_data %>% select(-time) %>% group_by(date) %>%
        summarize_all(funs((n() == samples_per_day) * mean(.)))

    # find k nearest neighbors for each day index
    nearest_neighbors = top_k(select(daily_averages, -date), k=knn, minobs=3)

    filled = fill_missing(input_data, nearest_neighbors, daily_averages,
        date_index, maxspan_days, samp=samples_per_day)

    return(filled)
}
