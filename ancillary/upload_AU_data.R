library(RMariaDB)
library(DBI)
library(stringr)
library(dplyr)
library(openxlsx)
library(tidyr)
library(zoo)
# library(httr)
# library(jsonlite)

#connect to MySQL
setwd('/home/aaron/sp/')
# setwd('/home/mike/git/streampulse/server_copy/sp')

conf = readLines('config.py')
extract_from_config = function(key){
    ind = which(lapply(conf, function(x) grepl(key, x)) == TRUE)
    val = str_match(conf[ind], '.*\\"(.*)\\"')[2]
    return(val)
}
pw = extract_from_config('MYSQL_PW')

con = dbConnect(RMariaDB::MariaDB(), dbname='sp',
    username='root', password=pw)

#read in and format daily discharge data; this will be converted to per second
setwd('/home/aaron/misc_data/australia_import/')
# setwd('/home/mike/Dropbox/streampulse/data/australia_data/')

Qnames = names(read.xlsx('daily_q_2014-17_LTIMsites.xlsx', 'Sheet1', rows=4))
# Qnames[Qnames == 'Barham.Bridge,.Noorong'] = 'Barham.Bridge'
Q = read.xlsx('daily_q_2014-17_LTIMsites.xlsx', 'Sheet1', startRow=5,
    detectDates=TRUE)
colnames(Q) = Qnames
Qdates = Q$LTIM.Site.Name

#this maps sites from series filenames to sites in discharge file
namemap = c('Hopwood'="Hopwood,.Windra Vale",
    'Windra Vale'="Hopwood,.Windra.Vale",'Widgee, Wakool River1'='Widgee',
    'Tralee1'='Tralee,.Cummins,.Llanos.Park',
    # 'Cummins'='Tralee,.Cummins,.Llanos.Park', #no location data
    'Llanos Park2'='Tralee,.Cummins,.Llanos.Park',
    'Barham Bridge'='Barham.Bridge,.Noorong',
    'Noorong2'='Barham.Bridge,.Noorong','Moss Road'="Moss.Rd/Day's.Rd",
    'McCoys Bridge'='McCoys','Loch Garry Gauge'='McCoys',
    "Darcy's Track"="Darcy's.Track",'WAL'='(Walbundry)', # this is actually Wallanthery
    #no discharge data for loch garry's gauge
    'LB'="Lane's.Bridge,.Cowl.Cowl",
    'WB'='Whealbah','CC'="Lane's.Bridge,.Cowl.Cowl")

#read in location data
latlongs = read.csv('au_latlongs.csv', stringsAsFactors=FALSE)

# #query Ask Geo database for UTC offsets associated with each site's lat/long
# accnt_id = '2023'
# api_key = extract_from_config('ASKGEO_KEY')
# latlongs_f = paste(paste(latlongs$Latitude, latlongs$Longitude, sep='%2C'),
#      collapse='%3B')
# askGeoReq_base = paste0('https://api.askgeo.com/v1/', accnt_id, '/',
#     api_key, '/query.json?databases=Point%2CTimeZone&points=', latlongs_f)
# askGeoReq_summer = paste0(askGeoReq_base, '&dateTime=2018-06-21')
# askGeoReq_winter = paste0(askGeoReq_base, '&dateTime=2018-12-21')
#
# r = httr::GET(askGeoReq_summer)
# json = httr::content(r, as="text", encoding="UTF-8")
# d_summer = jsonlite::fromJSON(json)
#
# r = httr::GET(askGeoReq_winter)
# json = httr::content(r, as="text", encoding="UTF-8")
# d_winter = jsonlite::fromJSON(json)

# saveRDS(d_winter$data, 'winter_utc_offsets.rds')
# saveRDS(d_summer$data, 'summer_utc_offsets.rds')
# d_winter = d_summer = list()
# d_winter$data = readRDS('winter_utc_offsets.rds')
# d_summer$data = readRDS('summer_utc_offsets.rds')
#
# #convert offsets to hours, bind with latlongs, bind with latlong data
# d_summer = d_summer$data$TimeZone %>%
#     mutate(UTCOffset=CurrentOffsetMs / 1000 / 60 / 60) %>%
#     select(UTCOffset) %>% bind_cols(d_summer$data$Point)
# d_winter = d_winter$data$TimeZone %>%
#     mutate(UTCOffset=CurrentOffsetMs / 1000 / 60 / 60) %>%
#     select(UTCOffset) %>% bind_cols(d_winter$data$Point)
# latlongs = left_join(latlongs, d_summer, by=c('Latitude', 'Longitude'))
# latlongs = left_join(latlongs, d_winter, by=c('Latitude', 'Longitude'),
#     suffix=c('.summer', '.winter'))
#
# #?
# latlongs = latlongs[! duplicated(latlongs$Site_Name),]
# rownames(latlongs) = 1:nrow(latlongs)

#extract data from files and compile
dirs = list.dirs(recursive=FALSE)
dirs = dirs[dirs != './continuous_Q']
# d=dirs[2]; f=files[2]
dirs=dirs[2] ############REMOVE THIS
site_set = c()
for(d in dirs){

    files = list.files(d)
    files=files[2]#################REMOVE THIS
    for(f in files){

        cur_time = Sys.time()
        attr(cur_time, 'tzone') = 'UTC'

        #add site data to site table in db
        site_deets = stringr::str_match(f, '(.+)_(\\d+)_(.+)_(\\d+).csv')

        if(site_deets[,4] %in% latlongs$Site_Name){
            print(paste(site_deets[,4]))
        } else {
            cat(paste('\tno latlong data for', site_deets[,4], '\n'))
            next
        }

        latlong_ind = which(latlongs$Site_Name == site_deets[4])

        site_code = paste0(gsub(' ', '-', site_deets[4]), '-', site_deets[3])
        site_up_df = data.frame('region'='AU',
            'site'=site_code,
            'name'=paste0(site_deets[2], ': ', site_deets[4], ' ', site_deets[3]),
            'latitude'=latlongs$Latitude[latlong_ind],
            'longitude'=latlongs$Longitude[latlong_ind],
            'usgs'=NA, 'addDate'=cur_time,
            'embargo'=0, 'by'=-903,
            'contact'='Michael Grace', 'contactEmail'='michael.grace@monash.edu')


        if(! site_code %in% site_set){
            dbWriteTable(con, 'site', site_up_df, append=TRUE)
        }
        site_set = append(site_set, site_code)

        #read in series data; deal with 3 possible date formats; convert tz
        data = read.csv(paste0(d, '/', f), stringsAsFactors=FALSE)
        data$DateTime_UTC = as.POSIXct(paste(data$Date, data$Time),
            format='%d/%m/%Y %H:%M:%S', tz='Australia/NSW')

        if(is.na(data$DateTime_UTC[1])){
            data = read.csv(paste0(d, '/', f), stringsAsFactors=FALSE)
            data$DateTime_UTC = as.POSIXct(paste(data$Date, data$Time),
                format='%Y-%m-%dT %H:%M:%S', tz='Australia/NSW')
        }
        if(difftime(data$DateTime_UTC[1], data$DateTime_UTC[2]) != -10){
            data = read.csv(paste0(d, '/', f), stringsAsFactors=FALSE)
            data$Date = as.Date(data$Date)
            data$DateTime_UTC = as.POSIXct(paste(data$Date, data$Time),
                format='%Y-%m-%d %H:%M:%S', tz='Australia/NSW')
        }

        attr(data$DateTime_UTC, 'tzone') = 'UTC'

        #conversions, adjustments, cleanup of series data
        data = mutate(data, 'AirPres_kPa'=atmo.pressure * 101.325) %>%
            select(-Date, -Time, -salinity, -atmo.pressure) %>%
            select(DateTime_UTC, 'Light_PAR'=I, 'WaterTemp_C'=tempC,
                'DO_mgL'=DO.meas, everything())
        data = data[! duplicated(data$DateTime_UTC),]

        #spread out daily discharge data so it covers all time points
        # curcol = which(Qnames == gsub(' ', '.', namemap[site_deets[4]]))
        # ML_m3s_conv_fac = 1e6 / (1000 * 86400)
        # discharge = Q[, curcol] * ML_m3s_conv_fac
        # data$date = as.Date(data$DateTime_UTC)
        # dateMatchInds = match(data$date, Qdates, nomatch=NA)
        # data$Discharge_m3s = discharge[dateMatchInds]
        # data = select(data, -date)

        #merge continuous discharge data, conform intervals
        q_dat = read.csv('continuous_Q/405232.csv', stringsAsFactors=FALSE,
            skip=8)

        q_dat = mutate(q_dat, Datetime=as.POSIXct(Datetime,
                format='%d/%m/%Y %H:%M:%S')) %>%
            filter(QC != 255, QC.1 != 255, ! is.na(Datetime)) %>%
            arrange()

        fiveseq = data.frame(Datetime=seq(q_dat$Datetime[1],
                q_dat$Datetime[nrow(q_dat)], by='5 min'))

        q_dat = select(q_dat, Datetime, Level_m=Water.Level..m..Mean,
                Discharge_MLd=Discharge..Ml.d..Mean) %>%
            mutate(Discharge_m3s=Discharge_MLd * (1e6 / (1000 * 86400))) %>%
            select(-Discharge_MLd) %>%
            right_join(fiveseq, by='Datetime')

        q_dat$Level_m = zoo::na.approx(q_dat$Level_m,
            na.rm=FALSE, rule=2)
        q_dat$Discharge_m3s = zoo::na.approx(q_dat$Discharge_m3s,
            na.rm=FALSE, rule=2)

        data = left_join(data, q_dat, by=c('DateTime_UTC'='Datetime'))

        #convert to long format, add remaining fields, insert into db
        data = gather(data, 'Variable', 'Value', Light_PAR:Discharge_m3s)
        data$region = 'AU'
        data$site = site_code
        data$flag = NA
        data$upload_id = -903

        dbWriteTable(con, 'data', data, append=TRUE)
        cat('\tgood\n')

    }
}


# OI!!! update site table with firstrecord, vars, etc




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

# #make df of dates, weekdays, months from minyear to maxyear in set
# mindate = as.Date(paste0(substr(min(d$DateTimeUTC, na.rm=TRUE), 1, 4),
#     '-01-01'))
# maxdate = as.Date(paste0(substr(max(d$DateTimeUTC, na.rm=TRUE), 1, 4),
#     '-12-31'))
# tm = data.frame(unix=mindate:maxdate)
# tm$date = as.Date(tm$unix, origin='1970-01-01')
# tm$wkday = strftime(tm$date, format='%a')
# tm$month = strftime(tm$date, format='%m')

# if(any(tm$date <= as.Date('1985-12-31'))){
#     write(paste('Error: pre-1986 date encountered; update DST handers'),
#         logfile, append=TRUE)
#     writeLines('1', 'scheduled_scripts/cuahsi/last_run_error.log')
#     stop()
# }
#
# #between 1986 and 2007, DST was from 1st sunday in apr to last sunday in oct
# #get DST start and END for each year in this range
# dst_starts = dst_ends = vector()
# old_dst = filter(tm, date < as.Date('2007-03-09'))
# if(nrow(old_dst)){
#     apr_sundays = old_dst[old_dst$month %in% '04' & old_dst$wkday == 'Sun',]
#     oct_sundays = old_dst[old_dst$month %in% '10' & old_dst$wkday == 'Sun',]
#     apr_sundays$cnt = unlist(tapply(rep(1, nrow(apr_sundays)),
#         substr(apr_sundays$date, 1, 7), cumsum))
#     oct_sundays$cnt = unlist(tapply(rep(1, nrow(oct_sundays)),
#         substr(oct_sundays$date, 1, 7), cumsum))
#     dst_starts = append(dst_starts,
#         apr_sundays$date[which(apr_sundays$cnt == 1)])
#     dst_ends = append(dst_ends,
#         oct_sundays$date[which(oct_sundays$cnt %in% 4:5 &
#                 substr(oct_sundays$date, 9, 10) %in% as.character(24:30))])
# }
#
# #as of 2007, DST is from 2nd sunday in march to 2nd sunday in nov
# #get DST start and END for each year in this range
# new_dst = filter(tm, date > as.Date('2007-03-09'))
# mar_sundays = new_dst[new_dst$month %in% '03' & new_dst$wkday == 'Sun',]
# nov_sundays = new_dst[new_dst$month %in% '11' & new_dst$wkday == 'Sun',]
# mar_sundays$cnt = unlist(tapply(rep(1, nrow(mar_sundays)),
#     substr(mar_sundays$date, 1, 7), cumsum))
# nov_sundays$cnt = unlist(tapply(rep(1, nrow(nov_sundays)),
#     substr(nov_sundays$date, 1, 7), cumsum))
# dst_starts = append(dst_starts,
#     mar_sundays$date[which(mar_sundays$cnt == 2)])
# dst_ends = append(dst_ends,
#     nov_sundays$date[which(nov_sundays$cnt == 2)])
#
# #determine utc offsets and local time for each date in dataset, join to set
# d = left_join(d, d_time, by=c('SiteCode'='site'))
# dst_starts_exp = paste0("as.POSIXct('", dst_starts, "')")
# dst_ends_exp = paste0("as.POSIXct('", dst_ends, "')")
# bool_exp = paste(
#     paste('(d$DateTimeUTC', dst_starts_exp, sep=' > '),
#     paste(dst_ends_exp, 'd$DateTimeUTC)', sep=' > '),
#     sep=' & ', collapse=' | ')
# d$UTCOffset = ifelse(eval(parse(text=bool_exp)),
#     d$UTCOffset.summer, d$UTCOffset.winter)
# d$LocalDateTime = d$DateTimeUTC + as.difftime(d$UTCOffset, units='hours')


