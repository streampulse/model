
# x = input_data[,8]; tol=12; algorithm='interpolation'
series_impute = function(x, tol, algorithm, ...){
    # records locations of NA runs longer than tol, imputes all gaps,
    # then replaces NAs for long runs.
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
    if(length(x) > samples_per_day * 2){ #don't leverage periodicity for <2 days
        x = ts(x, start=2, frequency=samples_per_day) #add periodicity info
    } else {
        warning(paste('Less than 2 days of data, so interpolating without',
            'periodicity information.'))
    }
    imputed = as.numeric(suppressWarnings(na.seasplit(x,
        algorithm=algorithm, ...)))

    #restore NAs where there were long runs
    if(length(na_locations)){
        imputed[long_na_run_indices] = NA
    }

    return(imputed)
}

# find similar k days
# df=select(daily_averages, -date); k=3; minobs=2
top_k = function(df, k, minobs=3){
    # df is dataframe of daily values with a row for each date and a colum for each variable
    # k is number of similar days to find
    # minobs is the min. observations required for a day
    D = matrix(NA, nrow(df), k) # empty dataframe for matches for each day
    # rows with missing data, exclude first and last
    cc = complete.cases(df)
    narows = which(!cc)

    if(length(narows)){ #skip the rest if there are no missing data
        narows = narows[1:(length(narows)-1)]
        # for each row with missing data, find similar rows without missing data
        for (i in narows){
            # skip days with less than n observations. hard to say anything about similarity
            if(sum(!is.na(df[i,])) < minobs){ next } # alternative: skip if all NA
            # difference between the observations and the rest of the dataframe
            daily_deltas = t(as.numeric(df[i,]) - t(df))
            # initialize NA row indexes
            delt = rep(NA, nrow(df))
            # candidate rows are only rows without NA
            for(j in which(cc)){
                tt = daily_deltas[j,] # the candidate row
                # skip comparisons with all NA rows (or self) or where limited comparison data
                if(j == i || all(is.na(tt))){ next }
                wn = which(is.na(tt)) # which are NA
                # euclidean distance, not normalized by the length of the vectors
                # (because we don't care about absolute distance)
                delt[j] = tt[-wn] %*% tt[-wn]
            }
            # get top k matches
            D[i,] = order(delt)[1:k]
        }
    }
    # returns a matrix with rows for each day and columns with k nearest neighbor days
    D
}

# data prep function for fill_gaps
# adds snap points, gets average data
prep_missing = function(df, nearest_neighbors, daily_averages, mm){
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
            missing[c(1,nrow(missing)),-c(1,2)] = apply(missing[c(1,nrow(missing)),-c(1,2)], 2, mean, na.rm=T)
        }
    }
    ndays = length(mm)
    ### SIMILAR DATA
    # grab similar days - pairs of missing day and similar day
    ss = data.frame(date=daily_averages$date[mm], match=daily_averages$date[t(nearest_neighbors[mm,])])
    ss = ss[complete.cases(ss),]
    similar = left_join(ss, df, by=c("match"="date")) %>%
        select(-match) %>% group_by(date, time) %>% summarize_all(mean) %>% ungroup()
    # make sure that the dates in similar and missing line up
    missing = right_join(missing, select(similar,date,time), by=c("date","time"))
    ### DAILY SNAP POINTS
    # add snap points at beginning/end each new day to rescale and match the daily trends
    newdaypoints = which(missing$time %in% missing$time[c(1,nrow(missing))])
    daypoints = missing[newdaypoints,]
    if(any(is.na(daypoints))){
        dayfill = linear_fill(select(daypoints,-date,-time))
        missing[newdaypoints,] = data.frame(date=daypoints$date, time=daypoints$time, dayfill)
    }
    list(missing=select(missing,-date,-time), similar=select(similar,-date,-time), index=select(similar,date,time))
}

# df=input_data; lim=0
fill_missing = function(df, nearest_neighbors, daily_averages,
    date_index, maxspan_days, lim=0){
    # df is the input data frame, all numeric data
    # nearest_neighbors are the similar days for each day
    # daily_averages is the daily data
    # date_index are the identifiers
    # maxspan_days is the maximum number of days to gap fill
    # lim is the minimum number of days to fill
    #     will not fill gaps that are less than this
    #     used for testing (b/c the test data have pre-existing gaps)
    #
    # the days that need filling in daily_averages
    filld = which(complete.cases(nearest_neighbors))

    if(length(filld)){ #skip the rest if no data are missing
        # groups for blocks of missing data
        group = cumsum(c(T, diff(filld) > 1))
        for(g in unique(group)){
            # grab missing days
            mm = filld[group==g]
            if( length(mm)>=lim && length(mm)<=maxspan_days ){
                pp = prep_missing(df, nearest_neighbors, daily_averages, mm)
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

# df=dd; maxspan_days=5; knn=3
gap_fill = function(df, maxspan_days=5, knn=3){
    # df is data frame, requires one column as POSIXct date time and the
    # other columns as numeric. the order of columns does not matter.

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
    samples_per_day = 24 * 60 / as.double(desired_int, units='mins')

    # impute in-line gaps (runs of NAs within a column) in df
    imputed = as.data.frame(sapply(X=input_data[,-(1:2)],
            FUN=series_impute, tol=12, algorithm='interpolation',
            simplify=FALSE)) #gaps >= tol will not be filled
    input_data = data.frame(select(input_data, date, time), imputed)

    # get averages for days with full sample coverage; otherwise NA (via mean())
    # nearly_complete_day = samples_per_day * 0.95 #could also use partial days
    daily_averages = input_data %>% select(-time) %>% group_by(date) %>%
        summarize_all(funs((n() == samples_per_day) * mean(.)))

    # find k nearest neighbors for each day index
    nearest_neighbors = top_k(select(daily_averages, -date), k=knn, minobs=2)

    filled = fill_missing(input_data, nearest_neighbors, daily_averages,
        date_index, maxspan_days)

    filled
}
