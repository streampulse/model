# Source gapfilling functions
gapfill_functions <- GET("https://raw.githubusercontent.com/streampulse/model/master/gapfill_functions.R")
eval(parse(text = content(gapfill_functions, as="text", encoding="UTF-8")),
    envir= .GlobalEnv)

# Source BASE functions
BASE_functions <- GET("https://raw.githubusercontent.com/streampulse/model/master/BASE_functions.R")
eval(parse(text = content(BASE_functions, as="text", encoding="UTF-8")),
    envir= .GlobalEnv)

retrieve_data <- function(sitecode, startdate=NULL, enddate=NULL, variables=NULL,
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
    variables <- c("DO_mgL","DOsat_pct","satDO_mgL","Level_m",
        "Depth_m","WaterTemp_C","Light_PAR","AirPres_kPa","Discharge_m3s")

    #assemble url based on user input
    u <- paste0("http://data.streampulse.org/api?sitecode=",sitecode)
    if(!is.null(startdate)) u <- paste0(u,"&startdate=",startdate)
    if(!is.null(enddate)) u <- paste0(u,"&enddate=",enddate)
    if(!is.null(variables)) u <- paste0(u,"&variables=",paste0(variables, collapse=","))
    if(flags) u <- paste0(u,"&flags=true")
    cat(paste0('URL: ',u,'\n'))

    #retrieve json; read into r object; format date
    if(is.null(token)){
        r <- httr::GET(u)
    }else{
        r <- httr::GET(u, httr::add_headers(Token = token))
    }
    json <- httr::content(r, as="text", encoding="UTF-8")
    d <- jsonlite::fromJSON(json)
    #d <- RJSONIO::fromJSON(json) # supposed to take care of NaN
    d$data$DateTime_UTC <- as.POSIXct(d$data$DateTime_UTC,tz="UTC")

    return(d)
}

sp_flags <- function(d){

  flag_df <- map_df(d$flags, data.frame) %>%
    as_tibble %>%
    select(id, region, site, variable, flag, startDate, endDate, comment) %>%
    mutate(startDate = as.POSIXct(startDate, format="%a, %d %b %Y %H:%M:%S", tz="UTC"),
           endDate = as.POSIXct(endDate, format="%a, %d %b %Y %H:%M:%S", tz="UTC"))

  return(flag_df)

}

# d=streampulse_data; model="streamMetabolizer"; type="bayes"
# fillgaps=TRUE; token=NULL
prep_metabolism <- function(d, model="streamMetabolizer", type="bayes",
    fillgaps=TRUE, token=NULL){
    #### Download and prepare data for metabolism modeling

    # type is one of "bayes" or "mle"
    # model is one of "streamMetabolizer" or "BASE"

    # Basic checks
    if(model=="BASE") type="bayes"

    #### Format data for models
    cat(paste("Formatting data for ",model,".\n", sep=""))
    dd <- d$data

    # Use USGS level and discharge if missing local versions
    if("USGSLevel_m" %in% dd$variable && !"Level_m" %in% dd$variable){
        dd$variable[dd$variable=="USGSLevel_m"] <- "Level_m"
    }
    if("USGSDischarge_m3s" %in% dd$variable && !"Discharge_m3s" %in% dd$variable){
        dd$variable[dd$variable=="USGSDischarge_m3s"] <- "Discharge_m3s"
    }

    vd <- unique(dd$variable) # variables
    dd <- tidyr::spread(dd, variable, value) # spread out data
    md <- d$sites # metadata

    # force into 15 minute intervals  ## MAKE THIS MORE VERSATILE
    alldates <- data.frame(DateTime_UTC=seq(dd[1,1],dd[nrow(dd),1],by="15 min"))
    dd <- full_join(alldates, dd, by="DateTime_UTC")

    # calculate/define model variables
    dd$solar.time <- suppressWarnings(streamMetabolizer::convert_UTC_to_solartime(
        date.time=dd$DateTime_UTC, longitude=md$lon[1], time.type="mean solar"))

    # estimate par
    apparentsolartime <- suppressWarnings(
        streamMetabolizer::convert_UTC_to_solartime(date.time=dd$DateTime_UTC,
            longitude=md$lon[1], time.type="apparent solar"))
    dd$light <- suppressWarnings(streamMetabolizer::calc_solar_insolation(
        app.solar.time=apparentsolartime, latitude=md$lat[1], format="degrees"))
    if("Light_PAR" %in% vd){ # fill in with observations if available
        dd$light[!is.na(dd$Light_PAR)] <- dd$Light_PAR[!is.na(dd$Light_PAR)]
    }else{
      cat("NOTE: Modeling PAR based on location and date.\n")
    }

    #### Gap filling
    # can specify maximum number of days to span
    # code found in gapfill_functions.R
    dd <- select(dd, -c(region, site, DateTime_UTC))
    if(fillgaps) dd <- gap_fill(dd, maxspan_days=5, knn=3)

    # Rename variables
    if("DO_mgL"%in%vd) dd$DO.obs <- dd$DO_mgL
    if("WaterTemp_C"%in%vd) dd$temp.water <- dd$WaterTemp_C
    if("Discharge_m3s"%in%vd) dd$discharge <- dd$Discharge_m3s
    if("Depth_m"%in%vd){
        dd$depth <- dd$Depth_m
    }else{ # get average depth from discharge
        dd$depth <- calc_depth(dd$discharge)
    }
    if("AirPres_kPa"%in%vd) dd$atmo.pressure <- dd$AirPres_kPa/101.325 # kPa to atm
    if(model=="streamMetabolizer"){
        if("satDO_mgL"%in%vd){
            dd$DO.sat <- dd$satDO_mgL
        }else{
            if("DOsat_pct"%in%vd){
                # define multiplicative factor, to catch if variable is fraction or percent... could do this on server first in future
                if(quantile(dd$DOsat_pct,0.9,na.rm=T)>10){ ff<-0.01 }else{ ff<-1 }
                dd$DO.sat <- dd$DO.obs/(dd$DOsat_pct*ff)
            }else{
                if(!all(c("temp.water","AirPres_kPa")%in%colnames(dd))){
                    stop("Insufficient data to fit this model.")
                }
                cat("NOTE: Modeling DO.sat based on water temperature and air pressure.\n")
                dd$DO.sat <- LakeMetabolizer::o2.at.sat.base(temp = dd$temp.water, baro = dd$AirPres_kPa*10, salinity = 0, model = 'garcia-benson')
            }
        }
    }

    # Select variables for model
    if(model=="BASE"){
      model_variables <- c("solar.time","DO.obs","temp.water","light","atmo.pressure")
    }else{ # streamMetabolizer
        model_variables <- c("solar.time","DO.obs","DO.sat","depth","temp.water","light")
        if(type=="bayes") model_variables <- c(model_variables,"discharge")
    }
    if(!all(model_variables%in%colnames(dd))){
        stop("Insufficient data to fit this model.")
    }

    # Structure data, add class for model name
    if(model=="BASE"){ # rename variables for BASE
        fitdata <- dd %>% select_(.dots=model_variables) %>%
            mutate(Date=as.Date(solar.time),Time=strftime(solar.time,format="%H:%M:%S"), salinity=0) %>%
            rename(I=light, tempC=temp.water, DO.meas=DO.obs) %>%
            select(Date, Time, I, tempC, DO.meas, atmo.pressure, salinity)
        BASE <- setClass("BASE", contains="data.frame")
        outdata <- as(fitdata, "BASE")
    }else if(model=="streamMetabolizer"){
        fitdata <- dplyr::select_(dd, .dots=model_variables)
        streamMetabolizer <- setClass("streamMetabolizer", representation(type="character"), contains="data.frame")
        outdata <- as(fitdata, "streamMetabolizer")
        outdata@type <- type
    }
    outdata
}

fit_metabolism <- function(fitdata){
    # check class of fitdata to determine which model to fit
    model <- class(fitdata)
    if(model=="streamMetabolizer") model_type <- fitdata@type
    # then, reset class of fitdata to data.frame, may not be necessary?
    class(fitdata) <- "data.frame"

    if(model=="streamMetabolizer"){
        # streamMetabolizer functions
        modname <- mm_name(type=model_type, pool_K600="binned",
            err_obs_iid=TRUE, err_proc_acor=FALSE, err_proc_iid=TRUE,
            ode_method = "trapezoid", deficit_src="DO_mod", engine="stan")
        modspecs <- specs(modname)
        modfit <- metab(specs = modspecs, data = fitdata)
        modfit
    }else if(model=="BASE"){
        tmp <- tempdir() # the temp dir for the data and model
        if(!dir.exists(tmp)) dir.create(tmp) # create if does not exist
        # create BASE directory
        directory <- file.path(tmp,"BASE")
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
        fitdata <- split(fitdata, fitdata$Date)
        # write individual date base files
        lapply(fitdata, function(xx){
            if(nrow(xx)==96 && all(complete.cases(xx))){ # only full days
                date <- unique(xx$Date)[1]
                write.csv(xx, file=paste0(directory,"/input/BASE_",date,".csv"), row.names=F)
            }
        })
        cat("Estimating metabolism with BASE.\n")
        # found in BASE_functions.R
        fit_BASE(directory=directory, interval=900, n.iter=30000, n.burnin=15000)
        structure(list(output_directory = directory), class="BASE") # return the BASE directory with class BASE
    }
}

predict_metabolism <- function(model_fit){
    if(class(model_fit)=="BASE"){
      directory <- model_fit$output_directory
      read.csv(paste0(directory,"/output/BASE_results.csv")) %>%
          separate(File, c("fileX", "date", "extX"), "_|\\.") %>%
          select(-fileX, -extX) %>% mutate(date=as.Date(date))
    }else{
      predict_metab(model_fit)
    }
}
