library(dplyr)
# as numeric shortcut
a.n <- function(x) as.numeric(x)

# linearly fill NAs
na_fill <- function(x, tol=Inf){
    # x is a vector of data
    # tol is max number of steps missing (if greater, it retains NA)
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

# linearly fill dataframe - works but could be optimized to not "apply"
linear_fill <- function(df, tol=Inf){
    # df is dataframe
    apply(df, 2, na_fill, tol=tol)
}

# find similar k days
# df=select(daily_averages, -date); k=3; minobs=2
top_k <- function(df, k, minobs=3){
    # df is dataframe of daily values with a row for each date and a colum for each variable
    # k is number of similar days to find
    # minobs is the min. observations required for a day
    D <- matrix(NA, nrow(df), k) # empty dataframe for matches for each day
    # rows with missing data, exclude first and last
    cc <- complete.cases(df)
    narows <- which(!cc)
    narows <- narows[1:(length(narows)-1)]
    # for each row with missing data, find similar rows without missing data
    for (i in narows){
        # skip days with less than n observations... hard to say anything about similarity
        if(sum(!is.na(df[i,])) < minobs){ next } # alternative: skip if all na
        # difference between the observations and the rest of the dataframe
        daily_deltas <- t(a.n(df[i,]) - t(df))
        # initialize NA row indexes
        delt <- rep(NA, nrow(df))
        # candidate rows are only rows without nan
        for(j in which(cc)){
            tt <- daily_deltas[j,] # the candidate row
            # skip comparisons with all NA rows (or self) or where limited comparison data
            if(j == i || all(is.na(tt))){ next }
            wn <- which(is.na(tt)) # which are NA
            # euclidean distance, not normalized by the length of the vectors
            #  (b/c we don't care about absolute distance)
            delt[j] <- tt[-wn] %*% tt[-wn]
        }
        # get top k matches
        D[i,] <- order(delt)[1:k]
    }
    # returns a matrix with rows for each day and columns with k nearest neighbor days
    D
}

# data prep function for fill_gaps
# adds snap points, gets average data
prep_missing <- function(df, nearest_neighbors, daily_averages, mm){
    # df is the data frame
    # daily_averages is the data frame of days
    # mm is the missing days
    # nearest_neighbors is the matching neighbors
    #
    ### MISSING DATA
    missing <- filter(df, date %in% daily_averages$date[mm])
    # if missing first and last obs, extend timeseries to include neighbor days
    #  unless first/last day
    if(any(!complete.cases(missing)[c(1,nrow(missing))])){
        if(mm[1]!=1) mm = c(mm[1]-1, mm) # if not first day, extend obs range
        if(tail(mm,1)!=nrow(daily_averages)) mm = c(mm, tail(mm,1)+1) # if not last day, extend obs range
        missing <- filter(df, date%in%daily_averages$date[mm]) # missing one step interpolation
        # if any data missing still (i.e., first and last day), add avg of first and last obs
        #  this should catch most NAs for filling missing
        if(any(!complete.cases(missing)[c(1,nrow(missing))])){
            missing[c(1,nrow(missing)),-c(1,2)] <- apply(missing[c(1,nrow(missing)),-c(1,2)], 2, mean, na.rm=T)
        }
    }
    ndays <- length(mm)
    ### SIMILAR DATA
    # grab similar days - pairs of missing day and similar day
    ss <- data.frame(date=daily_averages$date[mm], match=daily_averages$date[t(nearest_neighbors[mm,])])
    ss <- ss[complete.cases(ss),]
    similar <- left_join(ss, df, by=c("match"="date")) %>%
        select(-match) %>% group_by(date, time) %>% summarize_all(mean) %>% ungroup()
    # make sure that the dates in similar and missing line up
    missing <- right_join(missing, select(similar,date,time), by=c("date","time"))
    ### DAILY SNAP POINTS
    # add snap points at beginning/end each new day to rescale and match the daily trends
    newdaypoints <- which(missing$time %in% missing$time[c(1,nrow(missing))])
    daypoints <- missing[newdaypoints,]
    if(any(is.na(daypoints))){
        dayfill <- linear_fill(select(daypoints,-date,-time))
        missing[newdaypoints,] <- data.frame(date=daypoints$date, time=daypoints$time, dayfill)
    }
    list(missing=select(missing,-date,-time), similar=select(similar,-date,-time), index=select(similar,date,time))
}

# df=input_data; lim=0
fill_missing <- function(df, nearest_neighbors, daily_averages,
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
    filld <- which(complete.cases(nearest_neighbors))
    # groups for blocks of missing data
    group <- cumsum(c(T, diff(filld) > 1))
    for(g in unique(group)){
        # grab missing days
        mm <- filld[group==g]
        if( length(mm)>=lim && length(mm)<=maxspan_days ){
            pp <- prep_missing(df, nearest_neighbors, daily_averages, mm)
            dy <- (pp$missing - pp$similar)
            dyhat <- linear_fill(dy)
            filled <- pp$similar + dyhat
            df[which(paste(df$date,df$time) %in%
                    paste(pp$index$date,pp$index$time)),-c(1,2)] <- filled
        }
    }
    data.frame(date_index, select(df,-date,-time), stringsAsFactors=FALSE)
}

# df=dd; maxspan_days=5; knn=3
gap_fill <- function(df, maxspan_days=5, knn=3){
    # df is data frame, requires one column as POSIXct date time and the other columns as numeric
    #  - the order of columns does not matter
    # currently expects 15 min interval
    #  - in future can be expanded to check time interval and automatically fix

    # check if all but one column is numeric
    if( !(length(which(sapply(df, function(x) inherits(x, "numeric")))) ==
            ncol(df)-1) ){
        stop("ERROR: All but one column in df must be numeric.")
    }
    # check if a posix column exists
    wposix <- which(sapply(df, function(x) inherits(x, "POSIXct")))
    if( !(length(wposix)==1) ){
        stop("ERROR: Need at least one column in df with POSIXct datetime.")
    }

    # find POSIXct column, that is the one that we need to break into
    # date and time
    dtcol <- colnames(df)[wposix]

    # kind of goofy to do this by date and time, but that's because I
    # translated the code from Python
    input_data <- df %>% mutate(date=as.Date(df[,dtcol]),
        time=strftime(df[,dtcol], format="%H:%M:%S")) %>%
        select(-one_of(dtcol)) %>% select(date, time, everything())
    # index data
    date_index <- df %>% select(one_of(dtcol))

    # linearly interpolate df
    linearfill_data <- input_data %>% select(-date, -time) %>%
        linear_fill(tol=12)
    input_data <- data.frame(select(input_data, date, time), linearfill_data)

    # get daily averages if full day observations, otherwise NA
    daily_averages <- input_data %>% select(-time) %>% group_by(date) %>%
        summarize_all(funs((n()==96)*mean(.)))

    # find k nearest neighbors for each day index
    nearest_neighbors <- top_k(select(daily_averages, -date), k=knn, minobs=2)

    filled <- fill_missing(input_data, nearest_neighbors, daily_averages,
        date_index, maxspan_days)

    filled
}
