# # Source gapfilling functions
gapfill_functions = GET("https://raw.githubusercontent.com/streampulse/model/master/gapfill_functions.R")
eval(parse(text = content(gapfill_functions, as="text", encoding="UTF-8")),
    envir= .GlobalEnv)
#
# # Source BASE functions
BASE_functions = GET("https://raw.githubusercontent.com/streampulse/model/master/BASE_functions.R")
eval(parse(text = content(BASE_functions, as="text", encoding="UTF-8")),
    envir= .GlobalEnv)

request_data = function(sitecode, startdate=NULL, enddate=NULL, variables=NULL,
    flags=FALSE, token=NULL){
    # Download data from the streampulse platform

    # sitecode is a site name
    # startdate and enddate are YYYY-MM-DD strings, e.g., '1983-12-09'
    # variables is a vector of c('variable_one', ..., 'variable_n')
    # flags is logical, include flag data or not

    # Basic checks; make list of variables
    if(length(sitecode)>1) stop("Please only enter one site to model.")
    if(!is.null(startdate) & !is.null(enddate)){
        if(as.Date(enddate) < as.Date(startdate)){
            stop("Start date is after end date.")
        }
    }
    variables = c("DO_mgL","DOsat_pct","satDO_mgL","Level_m",
        "Depth_m","WaterTemp_C","Light_PAR","AirPres_kPa","Discharge_m3s")

    #assemble url based on user input
    u = paste0("http://data.streampulse.org/api?sitecode=",sitecode)
    if(!is.null(startdate)) u = paste0(u,"&startdate=",startdate)
    if(!is.null(enddate)) u = paste0(u,"&enddate=",enddate)
    if(!is.null(variables)) u = paste0(u,"&variables=",paste0(variables, collapse=","))
    if(flags) u = paste0(u,"&flags=true")
    cat(paste0('URL: ',u,'\n'))

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

    return(d)
}

sp_flags = function(d){

  flag_df = map_df(d$flags, data.frame) %>%
    as_tibble %>%
    select(id, region, site, variable, flag, startDate, endDate, comment) %>%
    mutate(startDate = as.POSIXct(startDate, format="%a, %d %b %Y %H:%M:%S", tz="UTC"),
           endDate = as.POSIXct(endDate, format="%a, %d %b %Y %H:%M:%S", tz="UTC"))

  return(flag_df)

}

#verify that datetimes from noaa are in utc and watch out for NA values (-9.96921e+36f)
# vars=c('windspeed', 'airpressure'); years = 2017
retrieve_pressure_wind = function(vars, years){

    # years = unique(c(substr(start_date,1,4), substr(end_date,1,4)))

    #format site data for use with geoknife package
    station = as.data.frame(t(d$sites[,c('lon','lat')]))
    station = simplegeom(station)

    if('windspeed' %in% vars){

        #setup and initializing
        components = c('vwnd', 'uwnd')
        datetime = .POSIXct(character()) #initialize empty POSIXct vector
        vwnd = uwnd = numeric()
        message(paste('Acquiring wind speed data for', length(years),
            'years. Each year takes a few minutes.'))

        #get data
        for(i in 1:length(years)){
            for(j in 1:length(components)){

                #get u and v components of windspeed from noaa
                fabric = webdata(url=paste0('https://www.esrl.noaa.gov/psd/th',
                    'redds/dodsC/Datasets/ncep.reanalysis/surface/',
                    components[j], '.sig995.', years[i], '.nc'),
                    variables=components[j])
                noaa_job = geoknife(stencil=station, fabric=fabric,
                    wait=TRUE)
                noaa_data = result(noaa_job, with.units=TRUE)

                if(j==1) datetime = append(datetime, noaa_data$DateTime)
                assign(components[j], append(get(components[j]), noaa_data$`1`))
                # -9.96921e+36f #missing data flagged with this value
            }
        }

        wnd = data.frame(datetime, vwnd, uwnd) %>%
            mutate(wind_speed=sqrt(uwnd^2 + vwnd^2)) %>%
            select(datetime, wind_speed)
        # attr(wnd$datetime, 'tzone') = 'UTC'
    }

    if('airpressure' %in% vars){ #same thing as above, but simpler

        # datetime2 = .POSIXct(character(), tz='GMT')
        # pres = numeric()
        pres = data.frame(datetime=.POSIXct(character()), pres=numeric())
        message(paste('Acquiring air pressure data for', length(years),
            'years. Each year takes a few minutes.'))

        for(i in 1:length(years)){

            fabric = webdata(url=paste0('https://www.esrl.noaa.gov/psd/th',
                'redds/dodsC/Datasets/ncep.reanalysis/surface/pres.sfc.',
                years[i], '.nc'), variables='pres')
            noaa_job = geoknife(stencil=station, fabric=fabric, wait=TRUE)
            noaa_data = result(noaa_job, with.units=TRUE)

            pres = rbind(pres, noaa_data[,c('DateTime','1')])
        }

        pres = data.frame(pres)
    }

    if('windspeed' %in% vars && 'airpressure' %in% vars){
        df_out = merge(pres, wnd, by.x='DateTime', by.y='datetime') %>%
            mutate(air_pressure=X1, DateTime_UTC=DateTime) %>% select(-X1)
    } else {
        if(vars == 'windspeed'){
            df_out = wnd
        } else df_out = pres %>%
                mutate(air_pressure = X1, DateTime_UTC=DateTime) %>% select(-X1)
    }

    return(df_out)
}

# d=streampulse_data; model="streamMetabolizer"; type="bayes"
# fillgaps=TRUE; interval='15 min'
# get_windspeed=TRUE; get_airpressure=TRUE
# d=streampulse_data; type='mle'; model='streamMetabolizer'; interval='15 min'; fillgaps='interpolation'
prep_metabolism = function(d, model="streamMetabolizer", type="bayes",
    interval='15 min', fillgaps='interpolation', ...){
    #, get_windspeed=FALSE,
    # get_airpressure=FALSE){

    #### format and prepare data for metabolism modeling

    # type is one of "bayes" or "mle"
    # model is one of "streamMetabolizer" or "BASE"
    # interval is the desired gap between successive observations. should be a
        # multiple of your sampling interval.
    # fillgaps must be one of the imputation methods available to
        # imputeTS::na.seasplit or 'none'
    # ... passes additional arguments to na.seasplit
    # get_windspeed and get_airpressure query NOAA's ESRL-PSD.
        # units are m/s and pascals, respectively

    # Basic checks
    if(model=="BASE") type="bayes" #can't use mle mode with BASE
    if(!grepl('\\d+ (min|hour)', interval, perl=TRUE)){ #correct interval format
        stop(paste('Interval must be of the form "length [space] unit" where',
            'length is numeric and unit is either "min" or "hour".'))
    }
    if(!fillgaps %in% c('interpolation','locf','mean','random','kalman','ma',
        'none')){
        stop(paste0("fillgaps must be one of 'interpolation', 'locf', 'mean', ",
            "'random', 'kalman', 'ma', or 'none'"))
    }

    #### Format data for models
    cat(paste("Formatting data for ",model,".\n", sep=""))
    dd = d$data

    # Use USGS level and discharge if missing local versions
    if("USGSLevel_m" %in% dd$variable && !"Level_m" %in% dd$variable){
        dd$variable[dd$variable=="USGSLevel_m"] = "Level_m"
    }
    if("USGSDischarge_m3s" %in% dd$variable && !"Discharge_m3s" %in% dd$variable){
        dd$variable[dd$variable=="USGSDischarge_m3s"] = "Discharge_m3s"
    }

    vd = unique(dd$variable) # variables
    dd = tidyr::spread(dd, variable, value) # spread out data
    md = d$sites # metadata

    # check for consistent sample interval
    run_lengths = rle(diff(as.numeric(dd$DateTime_UTC)))
    if(length(run_lengths$lengths) != 1){

        # if not consistent interval, set to mode
        uniqv = unique(run_lengths$values)
        input_int = uniqv[which.max(tabulate(match(run_lengths$values,
            uniqv)))] / 60
        warning(paste0('Sample interval is not consistent across dataset.',
            ' Gaps will be introduced! Using the most common interval: ',
            as.character(input_int), ' mins.'))
        input_int = as.difftime(input_int, unit='mins')
    } else {

        # if consistent, get the sample interval as a difftime object
        input_int = difftime(dd$DateTime_UTC[2],  dd$DateTime_UTC[1],
            units='mins')
    }

    # check if desired interval is compatible with sample interval
    int_parts = strsplit(interval, ' ')[[1]] #get num and str components
    desired_int = as.difftime(as.numeric(int_parts[1]),
        units=paste0(int_parts[2], 's')) #get desired interval as difftime
    remainder = as.double(desired_int, units='mins') %% as.double(input_int)
    if(remainder != 0){
        warning(paste0('Desired time interval (', interval,
            ') is not a multiple of sample interval (',
            as.character(as.numeric(input_int)),
            ' min). Gaps will be introduced!'))
    }

    # coerce to desired time interval
    alldates = data.frame(DateTime_UTC=seq.POSIXt(dd[1,1], dd[nrow(dd),1],
        by=interval))
    dd = left_join(alldates, dd, by='DateTime_UTC')

    # acquire additional variables if desired
    # if(get_windspeed || get_airpressure){
    #
    #     vars = character()
    #     if(get_windspeed) vars = append(vars, 'windspeed')
    #     if(get_airpressure) vars = append(vars, 'airpressure')
    #
    #     additional_vars = retrieve_pressure_wind(vars=vars, years=years)
    #     dd = left_join(dd, additional_vars, by='DateTime_UTC')
    #
    #     # linearly interpolate missing values for wind speed and air pressure
    #     dd$wind_speed = approx(x=dd$wind_speed, xout=which(is.na(dd$wind_speed)))$y
    #     dd$air_pressure = approx(x=dd$air_pressure,
    #         xout=which(is.na(dd$air_pressure)))$y
    # }

    # calculate/define model variables
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
      message("Estimating PAR based on location and date.")
    }

    # impute missing data. code found in gapfill_functions.R
    dd = select(dd, -c(region, site, DateTime_UTC))
    if(fillgaps != 'none') dd = gap_fill(dd, maxspan_days=5, knn=3,
        sint=desired_int, algorithm=fillgaps, ...)

    # Rename variables
    if("DO_mgL" %in% vd) dd$DO.obs = dd$DO_mgL
    if("WaterTemp_C" %in% vd) dd$temp.water = dd$WaterTemp_C
    if("Discharge_m3s" %in% vd) dd$discharge = dd$Discharge_m3s
    if("Depth_m" %in% vd){
        dd$depth = dd$Depth_m
    }else{ # get average depth from discharge
        dd$depth = calc_depth(dd$discharge)
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
                if(quantile(dd$DOsat_pct, 0.9, na.rm=T) > 10){
                    ff=0.01
                } else { ff=1 }
                dd$DO.sat = dd$DO.obs / (dd$DOsat_pct*ff)
            } else {
                if(!all(c("temp.water","AirPres_kPa") %in% colnames(dd))){
                    stop(paste('Insufficient data to fit this model. Need ',
                        'either DO sat (mgL) and water temp (C) OR air ',
                        'pressure (kPa).'))
                }
                message('Modeling DO.sat based on water temperature and',
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
        stop("Insufficient data to fit this model.")
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

        # streamMetabolizer functions
        modname = mm_name(type=model_type, pool_K600=pool_K600,
            err_obs_iid=TRUE, err_proc_acor=FALSE, err_proc_iid=proc_err,
            ode_method = "trapezoid", deficit_src="DO_mod", engine=engine)
        modspecs = specs(modname)
        modfit = metab(specs = modspecs, data = fitdata)
        modfit
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


