library(RMariaDB)
library(DBI)
library(stringr)
library(dplyr)
library(openxlsx)
library(tidyr)

#connect to MySQL
# setwd('/home/aaron/sp/scheduled_scripts/')
setwd('/home/mike/git/streampulse/server_copy/sp/scheduled_scripts/')

conf = readLines('../config.py')
extract_from_config = function(key){
    ind = which(lapply(conf, function(x) grepl(key, x)) == TRUE)
    val = str_match(conf[ind], '.*\\"(.*)\\"')[2]
    return(val)
}
pw = extract_from_config('MYSQL_PW')

con = dbConnect(RMariaDB::MariaDB(), dbname='sp',
    username='root', password=pw)

#read in and format daily discharge data; this will be converted to per second
setwd('/home/mike/Dropbox/streampulse/data/australia_data/')

Qnames = names(read.xlsx('daily_q_2014-17_LTIMsites.xlsx', 'Sheet1', rows=4))
Q = read.xlsx('daily_q_2014-17_LTIMsites.xlsx', 'Sheet1', startRow=5,
    detectDates=TRUE)
colnames(Q) = Qnames
Qdates = Q$LTIM.Site.Name

#this maps sites from series filenames to sites in discharge file
namemap = c('Hopwood'="Hopwood,.Windra Vale",
    'Windra Vale'="Hopwood,.Windra.Vale",'Widgee, Wakool River1'='Widgee',
    'Tralee1'='Tralee,.Cummins,.Llanos.Park',
    'Cummins'='Tralee,.Cummins,.Llanos.Park',
    'Llanos Park2'='Tralee,.Cummins,.Llanos.Park',
    'Barham Bridge'='Barham.Bridge',
    'Noorong2'='Noorong','Moss Road'="Moss.Rd/Day's.Rd",
    'McCoys Bridge'='McCoys','Loch Garry Gauge'='McCoys',
    "Darcy's Track"="Darcy's.Track",'WAL'='(Walbundry)',
    'LB'="Lane's.Bridge,.Cowl.Cowl",
    'WB'='Whealbah','CC'="Lane's.Bridge,.Cowl.Cowl")

dirs = list.dirs(recursive=FALSE)
# d=dirs[1]; f=files[4]
for(d in dirs){

    files = list.files(d)
    for(f in files){

        print(f)
        cur_time = Sys.time()
        attr(cur_time, 'tzone') = 'UTC'

        #add site to site table in db
        site_deets = stringr::str_match(f, '(.+)_(\\d+)_(.+)_(\\d+).csv')
        site_up_df = data.frame('region'='AU',
            'site'=paste0(site_deets[4], site_deets[3]),
            # 'site'=sitenm_dict[substr(f, 1, 6)][[1]],
            'name'=paste0(site_deets[2], ': ', site_deets[4], ' ', site_deets[3]),
            'latitude'=NA,
            'longitude'=NA,
            'usgs'=NA, 'addDate'=cur_time,
            'embargo'=0, 'by'=35,
            'contact'='Michael Grace', 'contactEmail'='michael.grace@monash.edu')

        # dbWriteTable(con, 'site', site_up_df, append=TRUE)

        #read in series data; deal with 3 possible date formats
        data = read.csv(paste0(d, '/', f), stringsAsFactors=FALSE)
        data$DateTime_UTC = as.POSIXct(paste(data$Date, data$Time),
            format='%d/%m/%Y %H:%M:%S', tz='UTC')

        if(is.na(data$DateTime_UTC[1])){
            data = read.csv(paste0(d, '/', f), stringsAsFactors=FALSE)
            data$DateTime_UTC = as.POSIXct(paste(data$Date, data$Time),
                format='%Y-%m-%dT %H:%M:%S', tz='UTC')
        }
        if(difftime(data$DateTime_UTC[1], data$DateTime_UTC[2]) != -10){
            data = read.csv(paste0(d, '/', f), stringsAsFactors=FALSE)
            data$Date = as.Date(data$Date)
            data$DateTime_UTC = as.POSIXct(paste(data$Date, data$Time),
                format='%Y-%m-%d %H:%M:%S', tz='UTC')
        }

        #conversions, adjustments, cleanup of series data
        data = mutate(data, 'AirPres_kPa'=atmo.pressure * 101.325) %>%
            select(-Date, -Time, -salinity, -atmo.pressure) %>%
            select(DateTime_UTC, 'Light_PAR'=I, 'WaterTemp_C'=tempC,
                'DO_mgL'=DO.meas, everything())
        data = data[! duplicated(data$DateTime_UTC),]

        #spread out daily discharge data so it covers all time points
        curcol = which(Qnames == namemap[site_deets[4]])
        ML_m3s_conv_fac = 1e6 / (1000 * 86400)
        discharge = Q[, curcol] * ML_m3s_conv_fac
        data$date = as.Date(data$DateTime_UTC)
        dateMatchInds = match(data$date, Qdates, nomatch=NA)
        data$Discharge_m3s = discharge[dateMatchInds]
        data = select(data, -date)

        #STILL GOTTA HANDLE TZ CONVERSION FROM LOCAL

        #convert to long format and insert into db
        data = gather(data, 'Variable', 'Value', Light_PAR:Discharge_m3s)
        # dbWriteTable(con, 'data', na_filt, append=TRUE)

        # OI!!! update site table with firstrecord, vars, etc


        # x = rle(as.numeric(data$DateTime_UTC))
        # if(! all(x$lengths == 1)){
        #     stop('not all same')
        # }
        # if(difftime(data$DateTime_UTC[1], data$DateTime_UTC[2]) != -10){
        #     stop('not 10')
        # }

    }
}

# site_data = data.frame('region'=rep(site_resp$data$stateCode, 2),
#     'site'=c(paste0(site_resp$data$siteCode, '-up'),
#         paste0(site_resp$data$siteCode, '-down')),
#     'name'=c(paste(site_resp$data$siteDescription, 'Upstream'),
#         paste(site_resp$data$siteDescription, 'Downstream')),
#     'latitude'=sensor_pos$referenceLatitude,
#     'longitude'=sensor_pos$referenceLongitude,
#     'usgs'=rep(NA, 2), 'addDate'=rep(cur_time, 2),
#     'embargo'=rep(0, 2), 'by'=rep(-900, 2),
#     'contact'=rep('NEON', 2), 'contactEmail'=rep(NA, 2))
#
# dbWriteTable(con, 'site', site_data, append=TRUE)

