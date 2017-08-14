# na_fill <- function(x, tol){
#   # x is a vector of data
#   # tol is max number of steps missing (if greater, it retains NA)
#   ina = is.na(x)
#   csum = cumsum(!ina)
#   wg = as.numeric(names(which(table(csum) > tol))) # which gaps are too long
#   x[ina] = approx(x, xout=which(ina))$y
#   x[which(csum%in%wg)[-1]] = NA
#   return(x)
# }

sp_data <- function(sitecode, startdate=NULL, enddate=NULL, variables=NULL, flags=FALSE, token=NULL){
    # sitecode is a site name
    # startdate and enddate are YYYY-MM-DD strings, e.g., '1983-12-09'
    # variables is a vector of c('variable_one', ..., 'variable_n')
    # flags is logical, include flag data or not
    u <- paste0("http://data.streampulse.org/api?sitecode=",sitecode)
    if(!is.null(startdate)) u <- paste0(u,"&startdate=",startdate)
    if(!is.null(enddate)) u <- paste0(u,"&enddate=",enddate)
    if(!is.null(variables)) u <- paste0(u,"&variables=",paste0(variables, collapse=","))
    if(flags) u <- paste0(u,"&flags=true")
    cat(paste0('URL: ',u,'\n'))
    if(is.null(token)){
        r <- httr::GET(u)
    }else{
        r <- httr::GET(u, httr::add_headers(Token = token))
    }
    json <- httr::content(r, as="text", encoding="UTF-8")
    d <- jsonlite::fromJSON(json)
    d$data$DateTime_UTC <- as.POSIXct(d$data$DateTime_UTC,tz="UTC")
    return(d)
}

sp_prep_metab <- function(d, model="streamMetabolizer", fillgaps=T){
    dd <- d$data
    # rename USGSDepth_m and USGSDischarge_m3s
    if("USGSLevel_m"%in%dd$variable && !"Level_m"%in%dd$variable){
        dd$variable[dd$variable=="USGSLevel_m"] <- "Level_m"
    }
    if("USGSDischarge_m3s"%in%dd$variable && !"Discharge_m3s"%in%dd$variable){
        dd$variable[dd$variable=="USGSDischarge_m3s"] <- "Discharge_m3s"
    }
    vd <- unique(dd$variable)
    dd <- tidyr::spread(dd, variable, value) # need to reshape...
    # check if sufficient data
    md <- d$sites # metadata
    # force into 15 minute intervals
    alldates <- data.frame(DateTime_UTC=seq(dd[1,1],dd[nrow(dd),1],by="15 min"))
    dd <- full_join(alldates,dd, by="DateTime_UTC")
    # calculate/define model variables
    dd$solar.time <- suppressWarnings(streamMetabolizer::convert_UTC_to_solartime(date.time=dd$DateTime_UTC, longitude=md$lon[1], time.type="mean solar"))
    # estimate par
    apparentsolartime <- suppressWarnings(streamMetabolizer::convert_UTC_to_solartime(date.time=dd$DateTime_UTC, longitude=md$lon[1], time.type="apparent solar"))
    cat("NOTE: Modeling PAR based on location and date.\n")
    dd$light <- suppressWarnings(streamMetabolizer::calc_solar_insolation(app.solar.time=apparentsolartime, latitude=md$lat[1], format="degrees"))
    if("Light_PAR"%in%vd){ # fill in with observations
        dd$light[!is.na(dd$Light_PAR)] <- dd$Light_PAR[!is.na(dd$Light_PAR)]
    }
    # GAP FILLING
    if(fillgaps) dd <- gap_fill(dd)
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
    return(dd)
}

sp_data_metab <- function(sitecode, startdate=NULL, enddate=NULL, type="bayes", model="streamMetabolizer", fillgaps=TRUE, token=NULL){
    # return data for streamMetabolizer metab()
    # sitecode is a site name
    # startdate and enddate are strings "2016-12-11"
    # type is one of "bayes" or "mle"
    # model is one of "streamMetabolizer" or "BASE"
    if(model=="BASE") type="bayes"
    if(length(sitecode)>1) stop("Please only enter one site to model.")
    if(is.null(startdate)&is.null(enddate)){
        if(as.Date(enddate)<as.Date(startdate)) stop("Start date is after end date.")
    }
    # Add: check for type, decide on what variables to include
    variables <- c("DO_mgL","DOsat_pct","satDO_mgL","Level_m","Depth_m","WaterTemp_C","Light_PAR","AirPres_kPa","Discharge_m3s")
    #c("DO_mgL","satDO_mgL","Depth_m","WaterTemp_C","Light_PAR","AirPres_kPa")
    cat("Downloading data from StreamPULSE.\n")
    d <- sp_data(sitecode, startdate, enddate, variables, FALSE, token)
    cat(paste("Formatting data for ",model,".\n", sep=""))
    dd <- sp_prep_metab(d, model, fillgaps)
    if(model=="BASE"){
      model_variables <- c("solar.time","DO.obs","temp.water","light","atmo.pressure")
    }else{ # streamMetabolizer
        model_variables <- c("solar.time","DO.obs","DO.sat","depth","temp.water","light")
        if(type=="bayes") model_variables <- c(model_variables,"discharge")
    }
    if(!all(model_variables%in%colnames(dd))){
        stop("Insufficient data to fit this model.")
    }
    fitdata <- dplyr::select_(dd, .dots=model_variables)
    # gap fill linearly if less than 3h
    #fitdata = data.frame(solar.time=fitdata[,1],apply(fitdata[,2:ncol(fitdata)],2,na_fill,tol=12))
    return(fitdata)
}

fit_metabolism <- function(fitdata, model="streamMetabolizer", model_type="bayes"){
    # alternative, model="BASE"
    if(model=="BASE"){
        tmp <- tempdir() # the temp dir for the data and model
        if(!dir.exists(tmp)) dir.create(tmp) # create if does not exist
        # create BASE directory
        if(!dir.exists(file.path(tmp,"BASE"))){
          dir.create(file.path(tmp,"BASE"))
          # add input folder
          dir.create(file.path(tmp,"BASE","input"))
          # add output folder
          dir.create(file.path(tmp,"BASE","output"))
          # - add instantaneous rates folder
          dir.create(file.path(tmp,"BASE","output","instantaneous rates"))
          # - add validation plots folder
          dir.create(file.path(tmp,"BASE","output","validation plots"))
        }
        # add BASE_metab_model_v2.2.txt
        download.file("https://raw.githubusercontent.com/streampulse/BASE/master/BASE/BASE_metab_model_v2.2.txt",
            file.path(tmp,"BASE","BASE_metab_model_v2.2.txt"))
        # directory <- file.path(tmp,"BASE-master/")  # example
        prep_BASE(fitdata, tmp)
        fit_BASE(directory=tmp, interval=900)#, n.iter=30000, n.burnin=15000)
        preds <- predict_BASE(tmp)
    }else{
        # streamMetabolizer functions
        modname <- mm_name(type=model_type, pool_K600="binned",
            err_obs_iid=TRUE, err_proc_acor=FALSE, err_proc_iid=TRUE,
            ode_method = "trapezoid", deficit_src="DO_mod", engine="stan")
        modspecs <- specs(modname)
        modfit <- metab(specs = modspecs, data = fitdata)
        # Predictions and plots
        preds <- predict_metab(modfit)
    }
    preds
}
