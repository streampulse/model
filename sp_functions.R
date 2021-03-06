
request_data = function(sitecode, startdate=NULL, enddate=NULL, variables=NULL,
    flags=TRUE, token=NULL){
    # Download data from the streampulse platform

    # sitecode is a site name
    # startdate and enddate are YYYY-MM-DD strings, e.g., '1983-12-09'
    # variables is a vector of c('variable_one', ..., 'variable_n')
    # flags is logical, include flag data or not

    #basic checks
    if(length(sitecode)>1){
        stop("Please only enter one site to model.", call.=FALSE)
    }
    if(!is.null(startdate) & !is.null(enddate)){
        if(as.Date(enddate) < as.Date(startdate)){
            stop("Start date is after end date.", call.=FALSE)
        }
    }
    if(!all(is.null(variables)) & !all(is.character(variables))){
        stop('Argument to "variables" must be a character vector.', call.=FALSE)
    }
    if(all(is.null(variables))){
        variables = c('DO_mgL','DOsat_pct','satDO_mgL','WaterPres_kPa',
            'Depth_m','WaterTemp_C','Light_PAR','AirPres_kPa','Discharge_m3s')
        # 'Level_m', #level is currently not used, but could be soon
        cat(paste0('Requesting all variables potentially useful for ',
            'metabolism modeling.\n'))
    } else {
        variables = unique(variables) #just makin' sure
        message(paste0('You may omit the "variables" parameter to ',
            'automatically retrieve\n\tall variables necessary for metabolism ',
            'modeling.'))
    }

    #assemble url based on user input
    u = paste0("http://data.streampulse.org/api?sitecode=",sitecode)
    if(!is.null(startdate)) u = paste0(u,"&startdate=",startdate)
    if(!is.null(enddate)) u = paste0(u,"&enddate=",enddate)
    u = paste0(u,"&variables=",paste0(variables, collapse=","))
    if(flags) u = paste0(u,"&flags=true")
    cat(paste0('\nAPI call: ',u,'\n\n'))

    #retrieve json; read into r object; format date
    if(is.null(token)){
        r = httr::GET(u)
    }else{
        r = httr::GET(u, httr::add_headers(Token = token))
    }
    json = httr::content(r, as="text", encoding="UTF-8")
    d = jsonlite::fromJSON(json)
    #d = RJSONIO::fromJSON(json) # supposed to take care of NaN
    d$data$DateTime_UTC = as.POSIXct(d$data$DateTime_UTC,tz="UTC")

    #rearrange columns if flag data present
    if(flags){
        notflagcols = colnames(d$data)[which(!colnames(d$data) %in%
                c('flagtype','flagcomment'))]
        d$data = cbind(d$data[,notflagcols], d$data[,c('flagtype','flagcomment')])

        #remove unneeded flag table returned by API
        d$flags = NULL
    }

    #format and print variables acquired
    retrieved_vars = unique(d$data$variable)
    n_print_batches = ceiling(length(retrieved_vars) / 5)
    cat('Retrieved the following variables:\n')
    for(b in 1:n_print_batches){

        s = (5 * (b - 1)) + 1 #first index to grab

        if(b == n_print_batches){
            e = length(retrieved_vars)
        } else {
            e = s + 4
        }

        cat('\t', paste0(retrieved_vars[s:e], collapse=', '), '\n')
    }

    #do the same for those not acquired (if necessary)
    missing_vars = variables[! variables %in% retrieved_vars]
    if(length(missing_vars)){

        n_print_batches = ceiling(length(missing_vars) / 5)
        cat('Could not find:\n')
        for(b in 1:n_print_batches){

            s = (5 * (b - 1)) + 1 #first index to grab

            if(b == n_print_batches){
                e = length(missing_vars)
            } else {
                e = s + 4
            }

            cat('\t', paste0(missing_vars[s:e], collapse=', '), '\n')
        }
    }

    return(d)
}

FindandCollect_airpres <- function(lat, long, start_datetime, end_datetime) {
    tf = tempfile()
    download.file("ftp://ftp.ncdc.noaa.gov/pub/data/noaa/isd-history.txt",tf,mode="wb")
    noaa.sites <- read.fwf(tf, skip = 22, header = F, widths = c(7,6,30, 5, 3, 6, 8, 9, 8, 9, 8), comment.char = "", col.names = c("USAF", "WBAN", "STATION NAME", "CTRY", "ST", "CALL", "LAT", "LON", "ELEV(M)", "BEGIN", "END"), flush = TRUE)
    noaa.sites <- na.omit(noaa.sites)
    noaa.sites <- noaa.sites %>%
        mutate(LAT = as.numeric(as.character(LAT))) %>%
        mutate(LON = as.numeric(as.character(LON))) %>%
        filter(LAT < (lat + 5) & LAT > (lat - 5) & LON < (long + 5) & LON > (long - 5))
    pt1 <- cbind(rep(long, length.out = length(noaa.sites$LAT)), rep(lat, length.out = length(noaa.sites$LAT)))
    pt2 <- cbind(noaa.sites$LON, noaa.sites$LAT)
    dist <- diag(distm(pt1, pt2, fun = distHaversine))/1000
    noaa.sites$dist <- dist
    tmp <- which((as.numeric(substr(noaa.sites$END,1,4)) >= as.numeric(substr(end_datetime, 1, 4))) & as.numeric(substr(noaa.sites$BEGIN,1,4)) <= as.numeric(substr(start_datetime, 1, 4)))
    noaa.sites <- noaa.sites[tmp,]
    noaa.sites <- noaa.sites[with(noaa.sites, order(dist)),]

    yrs <- seq(as.numeric(substr(start_datetime, 1, 4)),as.numeric(substr(end_datetime, 1, 4)), by = 1)
    for (i in 1:length(noaa.sites$dist)) {
        k <- i
        available <- vector(mode = 'logical', length = length(yrs))
        USAF <- as.character(noaa.sites$USAF[i])
        if(nchar(as.character(noaa.sites$WBAN[i])) == 5){
            WBAN <- as.character(noaa.sites$WBAN[i])
        } else {
            WBAN <- paste0(0,as.character(noaa.sites$WBAN[i]))
        }
        y <- as.data.frame(matrix(NA, nrow = 1, ncol = 12))
        for(j in 1:length(yrs)){
            tf = tempfile()
            download.file(paste0("ftp://ftp.ncdc.noaa.gov/pub/data/noaa/isd-lite/",yrs[j],"/",USAF,"-",WBAN,"-",yrs[j],".gz"),tf,mode="wb")
            x = read.table(tf)
            x[x==-9999] = NA
            if(length(which(!is.na(x$V7))) >= 0.9 * length(x$V7)) {
                available[j] <- TRUE
                y <- rbind(x,y)
            }else {
                break
            }
        }
        if(length(yrs) == length(which(available))){
            break
        }
    }
    y <- y[!is.na(y$V1),]
    colnames(y) = c("y","m","d","h","air_temp","dewtemp","air_kPa","winddir","sindspeed","skycover","precip1h","precip6h")
    y$air_kPa = y$air_kPa/100
    y$air_temp = y$air_temp/10
    y$DateTime_UTC = parse_datetime(paste0(y$y,"-",sprintf("%02d",y$m),"-",sprintf("%02d",y$d)," ",sprintf("%02d",y$h),":00:00 0"), "%F %T %Z")
    y <- y[with(y, order(DateTime_UTC)),]
    y = as_tibble(y) %>% select(DateTime_UTC,air_temp,air_kPa)
    ss = tibble(DateTime_UTC=seq(y$DateTime_UTC[1], y$DateTime_UTC[nrow(y)], by=900))
    xx = left_join(ss, y, by = "DateTime_UTC")
    xx = mutate(xx, air_temp=na.approx(air_temp), air_kPa=na.approx(air_kPa))
    daterng = c(start_datetime, end_datetime)
    xtmp = xx %>% filter(DateTime_UTC>=daterng[1] & DateTime_UTC<=daterng[2])
    select(xtmp, DateTime_UTC, air_kPa, air_temp)
    # print(noaa.sites[k,])
    return(select(xtmp, DateTime_UTC, air_kPa, air_temp))
}

retrieve_air_pressure = function(md, dd){

    # md2 <<- md
    # dd2 <<- dd
    # stop('a')
    # md = md2
    # dd = dd2

    lat = md$lat
    long = md$lon
    start_datetime = dd$DateTime_UTC[1]
    end_datetime = dd$DateTime_UTC[nrow(dd)]

    cat('Collecting air pressure data from NCDC (NOAA).\n')
    capture.output(df <- suppressWarnings(FindandCollect_airpres(lat, long,
        start_datetime, end_datetime)))
    cat('\n')

    df_out = df %>% mutate(AirPres_kPa = air_kPa) %>%
        select(DateTime_UTC, AirPres_kPa) %>% as.data.frame()

    return(df_out)
}

estimate_discharge = function(Z=NULL, Q=NULL, a=NULL, b=NULL,
    sh=NULL, dd=NULL, fit, plot){

    # dd2 <<- dd
    # stop('a')
    # dd = dd2
    # if(is.numeric(sh)){ #then need to calculate depth. based on:
        #https://web.archive.org/web/20170617070623/http://www.onsetcomp.com/files/support/tech-notes/onsetBCAguide.pdf

    if(is.null(plot)) plot = TRUE
    if(is.null(fit)) fit = 'power'
    if(plot){
        defpar = par()
        if(is.null(Z)){
            par(mfrow=c(1,1), mar=c(4,4,1,1), oma=c(0,0,0,0))
        } else {
            par(mfrow=c(2,1), mar=c(4,4,1,1), oma=c(0,0,0,0))
        }
    }

    if(! 'Depth_m' %in% colnames(dd)){

        cat(paste0('No depth data detected. Estimating level (AKA stage) from ',
            'water pressure.\n'))
        if(any(! c('WaterPres_kPa','AirPres_kPa') %in% colnames(dd))){
            stop(paste0('Air and/or water pressure not detected.',
                '\n\tNot enough information to proceed.'),
                call.=FALSE)
        }

        #remove air pressure so all pressure is from hydraulic head
        hyd_pres = dd$WaterPres_kPa - dd$AirPres_kPa

        if(! 'WaterTemp_C' %in% colnames(dd)){
            warning(paste0('Water temperature not detected.',
                '\n\tAssuming density of water is 1 g/cm^3.'),
                call.=FALSE)

            fl_dens = 1000 #g/m^3

        } else {

            #compute fluid density
            wat_temp = dd$WaterTemp_C
            T1 = 16.945176 * wat_temp
            T1b = 16.879850e-03 * wat_temp
            T2 = 7.9870401e-03 * wat_temp^2
            T3 = 46.170461e-06 * wat_temp^3
            T4 = 105.56302e-09 * wat_temp^4
            T5 = 280.54253e-12 * wat_temp^5
            fl_dens = (999.83952 + T1 - T2 - T3 + T4 - T5) / (1 + T1b)
        }

        #convert fluid density to lb/ft^3
        fl_dens = 0.0624279606 * fl_dens

        #convert hydraulic pressure to fluid depth
        ft_to_m = 0.3048
        kPa_to_psi = 0.1450377
        psi_to_psf = 144.0
        fl_depth = ft_to_m * (kPa_to_psi * psi_to_psf * hyd_pres) / fl_dens

    } else { #else we have depth already
        fl_depth = dd$Depth_m
    }

    #correct for sensor height above bed if sensor_height supplied
    if(!is.null(sh)){

        depth = fl_depth + sh #sh needs to be in meters

        message(paste0('Computing depth as level (AKA stage) plus ',
            'sensor height.\n\tMake sure other parameters supplied to ',
            'zq_curve are also based on depth,\n\trather than level,',
            ' or else omit sensor_height argument.'))

        dep_or_lvl = 'Depth'

        # cat('Quantiles of computed depth (m):\n')

    } else {

        depth = fl_depth

        message(paste0('Without sensor_height argument, ZQ rating curve will ',
            'be based on level (AKA stage),\n\trather than depth. Make sure ',
            'other parameters supplied to zq_curve are\n\talso based on level,',
            ' or else include sensor_height argument.'))

        dep_or_lvl = 'Level'

        # cat('Quantiles of computed level (m):\n')
    }
    # print(round(quantile(depth, na.rm=TRUE), 4))

    #generate rating curve if Z and Q sample data supplied
    if(!is.null(Z)){

        Q = Q[order(Z)]
        Z = Z[order(Z)] #just making sure there's no funny business

        # Q2 <<- Q
        # Z2 <<- Z
        # a2 <<- a
        # b2 <<- b
        # stop()
        # Q = Q2
        # Z = Z2
        # a = a2
        # b = b2

        #try to fit power model
        if(fit == 'power'){
            # mod = tryCatch(nls(Q ~ (a * Z^b), start=list(a=0.1, b=1)),
            #     error=function(e){
            #         stop(paste0('Failed to fit rating curve.\n\tThis is worth ',
            #             'mentioning to Mike: vlahm13@gmail.com.\n\t',
            #             'Note that you can fit your own curve and then supply\n\t',
            #             'a and b (of Q=aZ^b) directly.'), call.=FALSE)
            #     })
            mod = try(nls(Q ~ (a * Z^b), start=list(a=0.1, b=1)), silent=TRUE)
            if(class(mod) == 'try-error'){
                mod = try(nls(Q ~ (a * Z^b), start=list(a=1, b=1)), silent=TRUE)
                if(class(mod) == 'try-error'){
                    stop(paste0('Failed to fit rating curve.\n\tThis is worth ',
                        'mentioning to Mike: vlahm13@gmail.com.\n\t',
                        'Note that you can fit your own curve and then supply',
                        '\n\ta and b (of Q=aZ^b) directly.'), call.=FALSE)
                }
            }
        } else { #try to fit exponential
            if(fit == 'exponential'){
                mod = tryCatch(nls(Q ~ (a * exp(b * Z)),
                    start=list(a=0.01, b=1)),
                    error=function(e){
                        stop(paste0('Failed to fit rating curve. Try using ',
                            'fit="power" instead.\n\tAlso ',
                            'note that you can fit your own curve and then ',
                            'supply\n\ta and b (of Q=ae^(bZ)) directly.'),
                            call.=FALSE)
                    })
            } else { #fit linear
                mod = nls(Q ~ a * Z + b, start=list(a=1, b=0))
            }
        }

        if(plot == TRUE){
            plot(Z, Q, xlab='', ylab='Q samp. data (cms)',
                las=1, main=paste0('Rating curve fit (', fit, ')'))
            mtext('Z sample data (m)', 1, 2)
            plotseq = round(seq(Z[1], Z[length(Z)], diff(range(Z))/50), 2)
            lines(plotseq, predict(mod, list(Z=plotseq)))
        }

        params = summary(mod)$parameters
        a = params[1,1]
        b = params[2,1]

        #display info about curve
        eqns = list('power'='Q=aZ^b', 'exponential'='Q=ae^(bZ)',
            'linear'='Q=aZ+b')
        cat(paste0('Rating curve summary\n\tFit: ', fit, '\n\tEquation: ',
            eqns[fit], '\n\tParameters: a=', round(a, 3), ', b=',
            round(b, 3), '\n'))
        maxZ = round(max(Z, na.rm=TRUE), 2)
        maxD = round(max(depth, na.rm=TRUE), 2)
        if(maxD > maxZ){
            warning(paste0('Max observed ', tolower(dep_or_lvl), ' = ',
                maxD, '. Max observed input Z = ', maxZ, '.\n\tDischarge ',
                'estimates for ', tolower(dep_or_lvl), ' > ', maxZ,
                ' may be untrustworthy.'), call.=FALSE)
        }

    } #else a and b have been supplied directly

    #estimate discharge using a and b params from rating curve
    if(fit == 'power'){
        discharge = a * depth^b
    } else {
        if(fit == 'exponential'){
            discharge = a * exp(b * depth)
        } else {
            discharge = a * depth + b
        }
    }

    if(plot){
        plot(depth, discharge, xlab='',
            xlim=c(max(min(depth, na.rm=TRUE), 0), max(depth, na.rm=TRUE)),
            ylab='Est. Q (cms)', main='Rating curve prediction', las=1)
        mtext(paste(dep_or_lvl, 'series data (m)'), 1, 2)
    }
    suppressWarnings(par(defpar))

    return(discharge)
}

# rm_flagged=list('Bad Data', 'Questionable')
# d=streampulse_data; type='bayes'; model='streamMetabolizer'; interval='15 min'; fillgaps='interpolation'
# zq_curve=list(sensor_height=NULL, Z=NULL, Q=NULL, a=NULL, b=NULL,
#     fit='power', plot=TRUE)
# estimate_areal_depth=TRUE
prep_metabolism = function(d, model="streamMetabolizer", type="bayes",
    interval='15 min', rm_flagged=list('Bad Data', 'Questionable'),
    fillgaps='interpolation',
    zq_curve=list(sensor_height=NULL, Z=NULL, Q=NULL, a=NULL, b=NULL,
        fit='power', plot=TRUE),
    estimate_areal_depth=TRUE, ...){
    # zq_curve=list(Z=NULL, Q=NULL, a=NULL, b=NULL), ...){

    # type is one of "bayes" or "mle"
    # model is one of "streamMetabolizer" or "BASE"
    # interval is the desired gap between successive observations. should be a
        # multiple of your sampling interval.
    # rm_flagged is a list containing some combination of 'Bad Data',
        # 'Questionable', or 'Interesting'. Values corresponding to flags of
        # the specified type(s) will be replaced with NA, then interpolated
        # if fillgaps is not 'none'. rm_flagged can also be set to 'none' to
        # retain all flag information.
    # fillgaps must be one of the imputation methods available to
        # imputeTS::na.seasplit or 'none'
    # ... passes additional arguments to na.seasplit
    # get_windspeed and get_airpressure query NOAA's ESRL-PSD.
        # units are m/s and pascals, respectively

    # checks
    if(model=="BASE"){
        stop('BASE not yet supported', call.=FALSE) ####
        type="bayes" #can't use mle mode with BASE
    }
    if(!grepl('\\d+ (min|hour)', interval, perl=TRUE)){ #correct interval format
        stop(paste('Interval must be of the form "length [space] unit"\n\twhere',
            'length is numeric and unit is either "min" or "hour".'), call.=FALSE)
    }
    if(!fillgaps %in% c('interpolation','locf','mean','random','kalman','ma',
        'none')){
        stop(paste0("fillgaps must be one of 'interpolation', 'locf', 'mean',",
            "\n\t'random', 'kalman', 'ma', or 'none'"), call.=FALSE)
    }
    if(any(! rm_flagged %in% list('Bad Data','Questionable','Interesting')) &
        any(rm_flagged != 'none')){
        stop(paste0("rm_flagged must either be 'none' or a list containing any",
            " of:\n\t'Bad Data', 'Questionable', 'Interesting'."), call.=FALSE)
    }
    if(any(rm_flagged != 'none') & ! 'flagtype' %in% colnames(d$data)){
        stop(paste0('No flag data available.\n\t',
            'Call request_data again with flags=TRUE.'), call.=FALSE)
    }
    if(! 'list' %in% class(zq_curve)){
        stop('Argument "zq_curve" must be a list.', call.=FALSE)
    }

    ab_supplied = zq_supplied = FALSE
    using_zq_curve = !all(unlist(lapply(zq_curve[c('sensor_height',
        'a','b','Z','Q')], is.null)))
    if(using_zq_curve){

        #unpack arguments supplied to zq_curve
        sensor_height = Z = Q = a = b = fit = plot = NULL
        if(!is.null(zq_curve$sensor_height)){
            sensor_height = zq_curve$sensor_height
        }
        if(!is.null(zq_curve$Z)) Z = zq_curve$Z
        if(!is.null(zq_curve$Q)) Q = zq_curve$Q
        if(!is.null(zq_curve$a)) a = zq_curve$a
        if(!is.null(zq_curve$b)) b = zq_curve$b
        if(!is.null(zq_curve$fit)){
            fit = zq_curve$fit
            if(! fit %in% c('power', 'exponential', 'linear')){
                stop(paste0('Argument to "fit" must be one of: "power", ',
                    '"exponential", "linear".'), call.=FALSE)
            }
        }
        if(!is.null(zq_curve$plot)) plot = zq_curve$plot

        # message(paste0('NOTE: You have specified arguments to zq_curve.\n\t',
        #     'These are only needed if time-series data for discharge cannot',
        #     '\n\tbe found.'))
        if(is.numeric(a) & is.numeric(b)){
            ab_supplied = TRUE
        }
        if(length(Z) > 1 & length(Q) > 1){
            zq_supplied = TRUE
        }
        if(ab_supplied & zq_supplied){
            warning(paste0('Parameters (a, b) and data (Z, Q) supplied for ',
                'rating curve.\n\tOnly one set needed, so ignoring data.'),
                call.=FALSE)
            zq_supplied = FALSE
        } else {
            if(!ab_supplied & !zq_supplied){
                stop(paste0('Argument zq_curve must include either Z and Q as ',
                    'vectors of data\n\tor a and b as parameters of a rating ',
                    'curve.'), call.=FALSE)
            }
        }
    }

    #### Format data for models
    cat(paste("Formatting data for ",model,".\n", sep=""))
    dd = d$data

    #check for consistent sample interval (including cases where there are gaps
    #between samples and where the underying sample pattern changes)
    varz = unique(dd$variable)
    ints_by_var = data.frame(var=varz, int=rep(NA, length(varz)))
    for(i in 1:length(varz)){

        #get lengths and values for successive repetitions of the same
        #sample interval (using run length encoding)
        dt_by_var = sort(unique(dd$DateTime_UTC[dd$variable == varz[i]]))
        run_lengths = rle(diff(as.numeric(dt_by_var)))
        if(length(run_lengths$lengths) != 1){

            # if gaps or interval change, get mode interval
            uniqv = unique(run_lengths$values)
            input_int = uniqv[which.max(tabulate(match(run_lengths$values,
                uniqv)))] / 60 #this gets mode

            if(any(uniqv %% min(uniqv) != 0)){ #if underlying pattern changes
                warning(paste0('Sample interval is not consistent for ', varz[i],
                    '\n\tGaps will be introduced!\n\t',
                    'Using the most common interval: ',
                    as.character(input_int), ' mins.'), call.=FALSE)
            } else {
                message(paste0(length(run_lengths$lengths)-1,
                    ' sample gap(s) detected in ', varz[i], '.'))
            }

            #store the (most common) sample interval for each variable
            ints_by_var[i,2] = as.difftime(input_int, unit='mins')

        } else {

            # if consistent, just grab the diff between the first two times
            ints_by_var[i,2] = difftime(dt_by_var[2],  dt_by_var[1],
                units='mins')
        }
    }

    #will later coerce all vars to the longest sample interval
    if(length(unique(ints_by_var$int)) > 1){
        input_int = max(ints_by_var$int)
        message(paste0('Multiple sample intervals detected across variables (',
            paste(unique(ints_by_var$int), collapse=' min, '),
            ' min).\n\tUsing ', input_int, ' so as not to introduce gaps.'))
    } else {
        input_int = ints_by_var$int[1] #all the same
    }

    #remove (replace with NA) flagged data if desired
    if(any(rm_flagged != 'none')){
        flags_to_remove = unique(unlist(rm_flagged))
        dd$value[dd$flagtype %in% flags_to_remove] = NA
    }
    if('flagtype' %in% colnames(dd)){
        dd = subset(dd, select=-c(flagtype, flagcomment))
    }

    # Use USGS level and discharge if missing local versions
    #level currenty not used, but could be used in the absence of discharge and depth
    if("USGSLevel_m" %in% dd$variable && !"Level_m" %in% dd$variable){
        dd$variable[dd$variable=="USGSLevel_m"] = "Level_m"
    }
    if("USGSDischarge_m3s" %in% dd$variable && !"Discharge_m3s" %in% dd$variable){
        dd$variable[dd$variable=="USGSDischarge_m3s"] = "Discharge_m3s"
    }


    # dd2 <<- dd; d2 <<- d
    # dd = dd2; d = d2

    vd = unique(dd$variable) # variables
    dd = tidyr::spread(dd, variable, value) # spread out data
    md = d$sites # metadata

    #deal with "level"/"depth" naming issue. These measure the same thing.
    #UPDATE: level=stage=gage_height=vertical distance from sensor to surface
    #depth = vertical distance from bed to surface.
    #can be depth-at-gage OR averaged across width OR
    #averaged across area defined by width and O2 turnover distance.
    #this last metric is what the model actually needs
    #the code below may be party useful depending on what "depth" means for a given site
    # if('Level_m' %in% vd & ! 'Depth_m' %in% vd){ #use level if depth unavailable
    #     dd$Depth_m = dd$Level_m
    #     vd = c(vd, 'Depth_m')
    #     message('Using level data in place of missing depth data.')
    # }
    # if('Level_m' %in% vd & 'Depth_m' %in% vd){ #use col with more data if both
    #     na_cnt = sapply(dd[,c('Level_m','Depth_m')], function(x) sum(is.na(x)))
    #     level_or_depth = names(which.min(na_cnt))
    #     dd$Depth_m = dd[,level_or_depth]
    #     warning(paste0('Both level and depth data found. These measure the ',
    #         'same thing.\n\tUsing ', level_or_depth,
    #         ' because it has better coverage.'), call.=FALSE)
    # }

    #another test, now that depth has been acquired if available
    #UPDATE: turns out rating curves can directly relate level and discharge,
    #so depth is not needed at this point (though areal depth is needed, ultimately)
    #and sensor height is not needed as a result. This would only be useful
    #for correcting level to depth, which would be no better for buiding a rating curve
    missing_depth = ! 'Depth_m' %in% vd
    # missing_sens_height = !is.numeric(zq_curve$sensor_height)
    # if(using_zq_curve & missing_depth & missing_sens_height){
    # if(using_zq_curve & missing_depth){
    #     stop(paste0('Need either sensor_height or depth (level) time-series',
    #         ' data\n\tto compute discharge via rating curve.'),
    #         call.=FALSE)
    # }

    # check if desired interval is compatible with sample interval
    int_parts = strsplit(interval, ' ')[[1]] #get num and str components
    desired_int = as.difftime(as.numeric(int_parts[1]),
        units=paste0(int_parts[2], 's')) #get desired interval as difftime
    remainder = as.double(desired_int, units='mins') %% as.double(input_int)
    if(remainder != 0){
        message(paste0('Warning: Desired time interval (', interval,
            ') is not a multiple of sample interval (', #warning doesn't bubble
            as.character(as.numeric(input_int)), #up properly if imputation
            ' min).\n\tGaps will be introduced!')) #error occurs, so using message
    }

    #convert user-specified interval to minutes if it's in hours (because
    #fractional hours wont work with seq.POSIXct below)
    if(int_parts[2] == 'hour'){
        interval = paste(as.numeric(int_parts[1]) * 60, 'min')
    }

    #coerce to desired time interval. If multiple sample intervals are found...
    if(length(unique(ints_by_var$int)) > 1){

        #if some data are offset from the starting row, this may grab NAs instead
        #of data, so it iterates until finds the right starting row.
        na_props = 1
        starting_row = 1
        while(sum(na_props > 0.8) / length(na_props) > 0.4){ #heuristic test

            if(starting_row > 10){
                stop(paste0('Unable to coerce data to desired time interval.',
                    '\n\tTry specifying a different interval.'), call.=FALSE)
            }
            alldates = data.frame(DateTime_UTC=seq.POSIXt(dd[starting_row,1],
                dd[nrow(dd),1], by=interval))
            dd_temp = left_join(alldates, dd, by='DateTime_UTC')

            #get new NA proportions for each column
            na_props = apply(dd_temp[,-c(1:3)], 2,
                function(x){ sum(is.na(x)) / length(x) })
            starting_row = starting_row + 1
        }
        dd = dd_temp
    } else { #if just one sample interval is found, it's easy.
        alldates = data.frame(DateTime_UTC=seq.POSIXt(dd[1,1],
            dd[nrow(dd),1], by=interval))
        dd = left_join(alldates, dd, by='DateTime_UTC')
    }

    #acquire air pressure data if necessary
    missing_DOsat = all(! c('satDO_mgL','DOsat_pct') %in% vd)
    # missing_waterTemp = ! 'WaterTemp_C' %in% vd
    missing_airPres = ! 'AirPres_kPa' %in% vd
    need_airPres_for_DOsat = missing_DOsat & missing_airPres
    need_airPres_for_Q = using_zq_curve & missing_airPres & missing_depth #revisit "depth" here and elsewhere

    if(need_airPres_for_DOsat | need_airPres_for_Q){

        airpres = try(retrieve_air_pressure(md, dd), silent=TRUE)
        if(class(airpres) == 'try-error') {
            warning(paste('Failed to retrieve air pressure data.'), call.=FALSE)
        } else {
            dd = left_join(dd, airpres, by='DateTime_UTC')
            dd$AirPres_kPa = na.approx(dd$AirPres_kPa, na.rm=FALSE, rule=2)

            # linearly interpolate missing values for wind speed and air pressure
            # dd$wind_speed = approx(x=dd$wind_speed, xout=which(is.na(dd$wind_speed)))$y
            # dd$AirPres_kPa = approx(x=dd$AirPres_kPa,
            #     xout=which(is.na(dd$AirPres_kPa)))$y
        }

    }

    #correct any negative or 0 depth values (these break streamMetabolizer)
    if('Depth_m' %in% vd){
        if(any(na.omit(dd$Depth_m) <= 0)){
            warning('Depth values <= 0 detected. Replacing with 0.000001.',
                call.=FALSE)
            dd$Depth_m[dd$Depth_m <= 0] = 0.000001
        }
    }

    #estimate discharge from depth (or eventually level too) using Z-Q rating curve
    #if arguments have been supplied to zq_curve.
    if('Discharge_m3s' %in% vd & using_zq_curve){
        warning(paste0('Arguments supplied to zq_curve, so ignoring available',
            '\n\tdischarge time-series data.'), call.=FALSE)
    }
    if(zq_supplied){
        cat(paste0('Modeling discharge from rating curve.\n\tCurve will be ',
            'generated from supplied Z and Q samples.\n'))
        dd$Discharge_m3s = estimate_discharge(Z=Z, Q=Q, sh=sensor_height,
            dd=dd, fit=fit, plot=plot)
        vd = c(vd, 'Discharge_m3s')
    } else {
        if(ab_supplied){
            cat(paste0('Modeling discharge from rating curve determined by',
                '\n\tsupplied a and b parameters.\n'))
            dd$Discharge_m3s = estimate_discharge(a=a, b=b,
                sh=sensor_height, dd=dd, fit=fit, plot=plot)
            vd = c(vd, 'Discharge_m3s')
        }
        # if(ab_supplied & missing_depth){
        #     cat(paste0('Modeling discharge from rating curve determined by',
        #         '\n\tsupplied a and b parameters.\n\tEstimating depth or level',
        #         ' from water pressure, air pressure, and sensor height.\n'))
        #     dd$Discharge_m3s = estimate_discharge(a=zq_curve$a, b=zq_curve$b,
        #         sh=zq_curve$sensor_height, dd=dd)
        # } else {
        #     cat(paste0('Modeling discharge from rating curve determined by',
        #         '\n\tsupplied a and b parameters.\n\tUsing available ',
        #         'time-series data for depth.\n'))
        #     dd$Discharge_m3s = estimate_discharge(a=zq_curve$a, b=zq_curve$b,
        #         dd=dd)
        # }
    }

    # if(zq_supplied | ab_supplied){
    #     if('Discharge_m3s' %in% vd){
    #         warning(paste0('Discharge data available, so ignoring ',
    #             'zq_curve argument.'), call.=FALSE)
    #     } else {
    #         if(any(! c('AirPres_kPa', 'WaterPres_kPa', 'WaterTemp_C') %in% vd)){

    # if(! 'Discharge_m3s' %in% vd){

        # dd2 <<- dd
        # stop('a')
        # dd = dd2
    # }

    #correct any negative or 0 discharge values
    if('Discharge_m3s' %in% vd){
        if(any(na.omit(dd$Discharge_m3s) <= 0)){
            warning('Discharge values <= 0 detected. Replacing with 0.000001.',
                call.=FALSE)
            dd$Discharge_m3s[dd$Discharge_m3s <= 0] = 0.000001
        }
    }



    # convert UTC to solar time
    dd$solar.time = suppressWarnings(streamMetabolizer::convert_UTC_to_solartime(
        date.time=dd$DateTime_UTC, longitude=md$lon[1], time.type="mean solar"))

    # estimate par
    apparentsolartime = suppressWarnings(
        streamMetabolizer::convert_UTC_to_solartime(date.time=dd$DateTime_UTC,
            longitude=md$lon[1], time.type="apparent solar"))
    dd$light = suppressWarnings(streamMetabolizer::calc_solar_insolation(
        app.solar.time=apparentsolartime, latitude=md$lat[1], format="degrees"))
    if("Light_PAR" %in% vd){ # fill in with observations if available
        dd$light[!is.na(dd$Light_PAR)] = dd$Light_PAR[!is.na(dd$Light_PAR)]
    }else{
      cat("Estimating PAR based on latitude and time.\n")
    }

    # impute missing data. code found in gapfill_functions.R
    dd = select(dd, -c(region, site, DateTime_UTC))
    if(fillgaps != 'none') dd = gap_fill(dd, maxspan_days=5, knn=3,
        sint=desired_int, algorithm=fillgaps, ...)

    #rename variables
    if("DO_mgL" %in% vd) dd$DO.obs = dd$DO_mgL
    if("WaterTemp_C" %in% vd) dd$temp.water = dd$WaterTemp_C
    if('Discharge_m3s' %in% vd) dd$discharge = dd$Discharge_m3s

    #use discharge to estimate mean depth of area defined
    #by stream width and O2 turnover distance (if necessary or desired)
    if("Discharge_m3s" %in% vd & estimate_areal_depth){
        dd$depth = calc_depth(dd$discharge) #estimate mean areal depth
    } else {

        #otherwise just use depth directly.
        if("Depth_m" %in% vd){
            dd$depth = dd$Depth_m

            warning(paste0('Passing "Depth_m" values from StreamPULSE database',
                ' directly to streamMetabolizer.\n\tMetabolism estimates will',
                ' be best if "Depth_m" represents mean depth\n\tfor the ',
                'area defined by the width of the stream and the oxygen\n\t',
                'turnover distance. "Depth_m" may also represent depth-at-gage',
                ' or,\n\tworse, level-at-gage. These may result in poor ',
                'metabolism estimates.\n\tYou may want to use zq_curve.'),
                call.=FALSE)

        } else {

            if(estimate_areal_depth){
                stop(paste0('Missing discharge and depth data.\n\tNot enough ',
                    'information to proceed.\n\tMight parameter "zq_curve"',
                    ' be of service?'), call.=FALSE)
            } else {
                stop(paste0('Missing depth data.\n\tNot enough information to ',
                    'proceed.\n\tTry setting estimate_areal_depth to TRUE.'),
                    call.=FALSE)
            }
        }
    }

    # kPa to atm
    if("AirPres_kPa" %in% vd) dd$atmo.pressure = dd$AirPres_kPa / 101.325
    if(model=="streamMetabolizer"){
        if("satDO_mgL" %in% vd){
            dd$DO.sat = dd$satDO_mgL
        } else {
            if("DOsat_pct" %in% vd){

                # define multiplicative factor, to catch if variable is
                # fraction or percent... could do this on server first in future
                if(quantile(dd$DOsat_pct, 0.9, na.rm=TRUE) > 10){
                    ff=0.01
                } else { ff=1 }
                dd$DOsat_pct[dd$DOsat_pct == 0] = 0.000001 #prevent NaNs
                dd$DO.sat = dd$DO.obs / (dd$DOsat_pct*ff)
            } else {
                if(!all(c("temp.water","AirPres_kPa") %in% colnames(dd))){
                    stop(paste('Insufficient data to fit this model.\n\tNeed',
                        'either DO sat OR water temp and\n\tair ',
                        'pressure.'), call.=FALSE)
                }
                cat('Modeling DO.sat based on water temperature and',
                    'air pressure.\n')
                dd$DO.sat = LakeMetabolizer::o2.at.sat.base(temp=dd$temp.water,
                    baro=dd$AirPres_kPa*10, salinity=0, model='garcia-benson')
            }
        }
    }

    # Select variables for model
    if(model=="BASE"){
      model_variables = c("solar.time","DO.obs","temp.water","light",
          "atmo.pressure")
    }else{ # streamMetabolizer
        model_variables = c("solar.time","DO.obs","DO.sat","depth",
            "temp.water","light")
        if(type=="bayes") model_variables = c(model_variables,"discharge")
    }

    if(!all(model_variables %in% colnames(dd))){
        missing = model_variables[which(! model_variables %in% colnames(dd))]
        stop(paste0('Insufficient data to fit this model.\n\t',
            'Missing variable(s): ', paste0(missing, collapse=', ')),
            call.=FALSE)
    }

    # Structure data, add class for model name
    if(model=="BASE"){ # rename variables for BASE
        fitdata = dd %>% select_(.dots=model_variables) %>%
            mutate(Date=as.Date(solar.time),
                Time=strftime(solar.time, format="%H:%M:%S"), salinity=0) %>%
            rename(I=light, tempC=temp.water, DO.meas=DO.obs) %>%
            select(Date, Time, I, tempC, DO.meas, atmo.pressure, salinity)
        BASE = setClass("BASE", contains="data.frame")
        outdata = as(fitdata, "BASE")
    }else if(model=="streamMetabolizer"){
        fitdata = dplyr::select_(dd, .dots=model_variables)
        streamMetabolizer = setClass("streamMetabolizer",
            representation(type="character"), contains="data.frame")
        outdata = as(fitdata, "streamMetabolizer")
        outdata@type = type
    }
    outdata
}

fit_metabolism = function(fitdata){
    # check class of fitdata to determine which model to fit
    model = class(fitdata)
    if(model=="streamMetabolizer") model_type = fitdata@type
    # then, reset class of fitdata to data.frame, may not be necessary?
    class(fitdata) = "data.frame"

    if(model=="streamMetabolizer"){

        #choose appropriate model specifications based on model type
        if(model_type=='bayes'){
            engine = 'stan'; pool_K600='binned'; proc_err = TRUE
        } else {
            engine = 'nlm'; pool_K600='none'; proc_err = FALSE
        }

        # set most model specs
        modname = mm_name(type=model_type, pool_K600=pool_K600,
            err_obs_iid=TRUE, err_proc_acor=FALSE, err_proc_iid=proc_err,
            ode_method = "trapezoid", deficit_src="DO_mod", engine=engine)
        modspecs = specs(modname)

        # fitdata2 <<- fitdata
        # stop('a')
        # fitdata = fitdata2
        #get average log daily discharge and use it to parameterize k v. Q curve
        if(engine == 'stan'){
            addis = tapply(log(fitdata$discharge),
                substr(fitdata$solar.time,1,10), mean)
            # sum(is.infinite(log(fitdata$discharge))
            modspecs$K600_lnQ_nodes_centers = seq(from=min(addis),
                to=max(addis), length.out=7)
        }

        #fit model
        modfit = metab(specs=modspecs, data=fitdata)
        return(modfit)

    }else if(model=="BASE"){
        tmp = tempdir() # the temp dir for the data and model
        if(!dir.exists(tmp)) dir.create(tmp) # create if does not exist
        # create BASE directory
        directory = file.path(tmp,"BASE")
        if(!dir.exists(directory)){
            dir.create(directory)
            # add input folder
            dir.create(file.path(directory,"input"))
            # add output folder
            dir.create(file.path(directory,"output"))
            # - add instantaneous rates folder
            dir.create(file.path(directory,"output","instantaneous rates"))
            # - add validation plots folder
            dir.create(file.path(directory,"output","validation plots"))
        }
        # download BASE_metab_model_v2.2.txt
        download.file("https://raw.githubusercontent.com/streampulse/BASE/master/BASE/BASE_metab_model_v2.2.txt",
            file.path(directory,"BASE_metab_model_v2.2.txt"), quiet=TRUE)
        file.remove(list.files(file.path(directory,"input"),full.names=T)) # clear out input folder
        fitdata = split(fitdata, fitdata$Date)
        # write individual date base files
        lapply(fitdata, function(xx){
            if(nrow(xx)==96 && all(complete.cases(xx))){ # only full days
                date = unique(xx$Date)[1]
                write.csv(xx, file=paste0(directory,"/input/BASE_",date,".csv"), row.names=F)
            }
        })
        cat("Estimating metabolism with BASE.\n")
        # found in BASE_functions.R
        fit_BASE(directory=directory, interval=900, n.iter=30000, n.burnin=15000)
        structure(list(output_directory = directory), class="BASE") # return the BASE directory with class BASE
    }
}

predict_metabolism = function(model_fit){
    if(class(model_fit)=="BASE"){
      directory = model_fit$output_directory
      read.csv(paste0(directory,"/output/BASE_results.csv")) %>%
          separate(File, c("fileX", "date", "extX"), "_|\\.") %>%
          select(-fileX, -extX) %>% mutate(date=as.Date(date))
    }else{
      predict_metab(model_fit)
    }
}
