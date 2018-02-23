FindNearestAirPres <- function(lat, long, start_datetime, end_datetime) {
  require(geosphere)
  require(dplyr)
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
  noaa.sites$dist_km <- dist
  tmp <- which((as.numeric(substr(noaa.sites$END,1,4)) >= as.numeric(substr(end_datetime, 1, 4))) & as.numeric(substr(noaa.sites$BEGIN,1,4)) <= as.numeric(substr(start_datetime, 1, 4)))
  noaa.sites <- noaa.sites[tmp,]
  noaa.sites <- noaa.sites[with(noaa.sites, order(dist_km)),]
  
  yrs <- seq(as.numeric(substr(start_datetime, 1, 4)),as.numeric(substr(end_datetime, 1, 4)), by = 1)
  for (i in 1:length(noaa.sites$dist_km)) {
    k <- i
    available <- vector(mode = 'logical', length = length(yrs))
    USAF <- as.character(noaa.sites$USAF[i])
    WBAN <- if(nchar(as.character(noaa.sites$WBAN[i])) == 5){
      as.character(noaa.sites$WBAN[i])
    } else {
      paste0(0,as.character(noaa.sites$WBAN[i]))
    }
    for(j in 1:length(yrs)){
      tf = tempfile()
      download.file(paste0("ftp://ftp.ncdc.noaa.gov/pub/data/noaa/isd-lite/",yrs[j],"/",USAF,"-",WBAN,"-",yrs[j],".gz"),tf,mode="wb")
      x = read.table(tf)
      x[x==-9999] = NA
      if(length(which(!is.na(x$V7))) >= 0.9 * length(x$V7)) {
        available[j] <- TRUE
      }else {
        break
      }
    }
    if(length(yrs) == length(which(available))){
      break
    }
  }
  
  return(noaa.sites[k,])
}

CollectNearestAirPres <- function(USAF, WBAN, start_datetime, end_datetime)  {
  require(readr)
  require(tibble)
  require(zoo)
  require(dplyr)
  yrs <- seq(as.numeric(substr(start_datetime, 1, 4)),as.numeric(substr(end_datetime, 1, 4)), by = 1)
  y <- as.data.frame(matrix(NA, nrow = 1, ncol = 12))
  for(j in 1:length(yrs)){
    tf = tempfile()
    download.file(paste0("ftp://ftp.ncdc.noaa.gov/pub/data/noaa/isd-lite/",yrs[j],"/",USAF,"-",WBAN,"-",yrs[j],".gz"),tf,mode="wb")
    x = read.table(tf)
    y <- rbind(x, y)
  }
  y <- y[!is.na(y$V1),]
  y[y==-9999] = NA
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
}

FindandCollect_airpres <- function(lat, long, start_datetime, end_datetime) {
  require(geosphere)
  require(readr)
  require(tibble)
  require(zoo)
  require(dplyr)
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
  print(noaa.sites[k,])
  return(select(xtmp, DateTime_UTC, air_kPa, air_temp))
}

