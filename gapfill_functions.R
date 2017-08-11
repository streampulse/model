library(dplyr)
# as numeric shortcut
a.n <- function(x) as.numeric(x)

# gap fill
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
top_k <- function(df, k, minobs=3){
    # df is dataframe
    # k is number of similar days to find
    # minobs is the min. observations required for a day
    D <- matrix(NA, nrow(df), k) # empty dataframe for matches for each day
    # rows with missing data, exclude first and last
    cc <- complete.cases(df)
    narows <- which(!cc)
    narows <- narows[1:(length(narows)-1)]
    # for each row with missing data,
    #  find similar rows without missing data
    for (i in narows){
        # skip days with less than n observations... hard to say anything about similarity
        if(sum(!is.na(df[i,])) < minobs){ next } # alternative: skip if all na
        # difference between the observations and the rest of the dataframe
        xdelt <- t(a.n(df[i,]) - t(df))
        # initialize NA row indexes
        delt <- rep(NA, nrow(df))
        # candidate rows are only rows without nan
        for(j in which(cc)){
            tt <- xdelt[j,] # the candidate row
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
    D
}


# data prep function for fill_gaps
# adds snap points, gets average data
prep_missing <- function(df, ns, xd, mm){
    # df is the data frame
    # xd is the data frame of days
    # mm is the missing days
    # ns is the matching neighbors
    #
    ### MISSING DATA
    missing <- filter(df, date%in%xd$date[mm])
    # if missing first and last obs, extend timeseries to include neighbor days
    #  unless first/last day
    if(any(!complete.cases(missing)[c(1,nrow(missing))])){
        if(mm[1]!=1) mm = c(mm[1]-1, mm) # if not first day, extend obs range
        if(tail(mm,1)!=nrow(xd)) mm = c(mm, tail(mm,1)+1) # if not last day, extend obs range
        missing <- filter(df, date%in%xd$date[mm]) # missing one step interpolation
        # if any data missing still (i.e., first and last day), add avg of first and last obs
        #  this should catch most NAs for filling missing
        if(any(!complete.cases(missing)[c(1,nrow(missing))])){
            missing[c(1,nrow(missing)),-c(1,2)] <- apply(missing[c(1,nrow(missing)),-c(1,2)], 2, mean, na.rm=T)
        }
    }
    ndays <- length(mm)
    ### SIMILAR DATA
    # grab similar days - tuples of missing day and similar day
    ss <- data.frame(date=xd$date[mm], match=xd$date[t(ns[mm,])])
    ss <- ss[complete.cases(ss),]
    similar <- left_join(ss, df, by=c("match"="date")) %>%
        select(-match) %>% group_by(date, time) %>% summarize_all(mean) %>% ungroup()
    # make sure that the dates in similar and missing line up
    missing <- right_join(missing, select(similar,date,time), by=c("date","time"))
    ### DAILY SNAP POINTS
    # add snap points at each new day to rescale
    newdaypoints <- which(missing$time %in% missing$time[c(1,nrow(missing))])
    daypoints <- missing[newdaypoints,]
    if(any(is.na(daypoints))){
        dayfill <- linear_fill(select(daypoints,-date,-time))
        missing[newdaypoints,] <- data.frame(date=daypoints$date, time=daypoints$time, dayfill)
    }
    list(missing=select(missing,-date,-time), similar=select(similar,-date,-time), index=select(similar,date,time))
}

fill_missing <- function(df, ns, xd, ids, maxspan, lim=0){
    # df is the data frame
    # ns are the similar days for each day
    # xd is the daily data
    # ids are the identifiers
    # maxspan is the maximum number of days to gap fill
    # lim is the minimum number of days to fill
    #     will not fill gaps that are less than this
    #     used for testing (b/c the test data have pre-existing gaps)
    #
    # the days that need filling in xd
    filld <- which(complete.cases(ns))
    # groups for blocks of missing data
    group <- cumsum(c(T,diff(filld)>1))
    for(g in unique(group)){
        # grab missing days
        mm <- filld[group==g]
        if(length(mm)>=lim && length(mm)<maxspan){
            pp <- prep_missing(df, ns, xd, mm)
            dy <- (pp$missing - pp$similar)
            dyhat <- linear_fill(dy)
            filled <- pp$similar + dyhat
            df[which(paste(df$date,df$time)%in%paste(pp$index$date,pp$index$time)),-c(1,2)] <- filled
        }
    }
    data.frame(ids,select(df,-date,-time), stringsAsFactors=FALSE)
}

gap_fill <- function(df, maxspan_days=5, knn=3){
    df <- select(df, -DateTime_UTC)
    # require 15 min interval

    # kind of goofy to do this by date and time, but that's because I translated the code from Python
    dds <- df %>% mutate(date=as.Date(solar.time), time=strftime(solar.time, format="%H:%M:%S")) %>%
        select(-c(region, site, solar.time)) %>% select(date, time, everything())
    # index data
    ids <- df %>% select(c(region, site, solar.time))

    # linearly interpolate df
    ddl <- dds %>% select(-date, -time) %>% linear_fill(tol=12)
    dds <- data.frame(select(dds,date,time),ddl)

    # get daily averages if full day observations, otherwise NA
    xd <- dds %>% select(-time) %>% group_by(date) %>%
        summarize_all(funs((n()==96)*mean(.)))

    ns <- top_k(select(xd,-date), k=knn, minobs=2)

    filled <- fill_missing(dds, ns, xd, ids, maxspan=maxspan_days)

    filled
}

# test
# df <- read.csv("CT_Farmington.csv")
# df_filled <- gap_fill(df)
