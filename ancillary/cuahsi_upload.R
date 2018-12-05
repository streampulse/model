library(openxlsx)
library(RMariaDB)
library(DBI)
library(stringr)
library(dplyr)
library(httr)
library(jsonlite)
setwd('/home/mike/Dropbox/streampulse/data/cuahsi_connection/hydroserver_templates')
date = as.character(Sys.Date())

#connect to database
# conf = readLines('/home/aaron/sp/scheduled_scripts/')
conf = readLines('/home/mike/git/streampulse/server_copy/sp/config.py')
extract_from_config = function(key){
    ind = which(lapply(conf, function(x) grepl(key, x)) == TRUE)
    val = str_match(conf[ind], '.*\\"(.*)\\"')[2]
    return(val)
}
pw = extract_from_config('MYSQL_PW')
con = dbConnect(RMariaDB::MariaDB(), dbname='sp', username='root', password=pw)

#read in all site data from site table
res = dbSendQuery(con, paste("SELECT CONCAT(region, '_', site) AS site,",
    "ROUND(latitude, 5) AS lat, ROUND(longitude, 6) AS lon FROM site;"))
sites = dbFetch(res)
# sites = sites[-c(5, 55, 62:64),]
dbClearResult(res)

#query Ask Geo database for time zones associated with each lat/long
accnt_id = '2023'
api_key = extract_from_config('ASKGEO_KEY')
latlongs = paste(paste(sites$lat, sites$lon, sep='%2C'), collapse='%3B')
askGeoReq_base = paste0('https://api.askgeo.com/v1/', accnt_id, '/',
    api_key, '/query.json?databases=Point%2CTimeZone&points=', latlongs)
    # api_key, '/query.json?databases=Point&points=', latlongs)
askGeoReq_summer = paste0(askGeoReq_base, '&dateTime=2018-06-21')
askGeoReq_winter = paste0(askGeoReq_base, '&dateTime=2018-12-21')

r = httr::GET(askGeoReq_summer)
json = httr::content(r, as="text", encoding="UTF-8")
d_summer = try(jsonlite::fromJSON(json), silent=TRUE)

r = httr::GET(askGeoReq_winter)
json = httr::content(r, as="text", encoding="UTF-8")
d_winter = try(jsonlite::fromJSON(json), silent=TRUE)

#load excel worksheets into dataframes and export as csv
wb = loadWorkbook('NC_Advanced.xlsx')
shtnames = names(wb)[-(1:2)]
shtnames = shtnames[-which(shtnames == 'DataValues')]
for(s in shtnames){
    dat = read.xlsx('NC_Advanced.xlsx', s)
    dat = dat[-(1:4),-1]
    dir.create(date)
    write.csv(dat, paste0(date, '/', s, '.csv'))
}

#read data from database
res = dbSendQuery(con, paste("SELECT data.value AS DataValue,",
    "data.DateTime_UTC AS DateTimeUTC,",
    "CONCAT(data.region, '_', data.site) AS SiteCode, data.variable AS VariableCode,",
    "flag.flag AS QualityControlLevelCode FROM data LEFT JOIN flag ON",
    "data.flag=flag.id WHERE data.region='NC';"))
resout = dbFetch(res)
dbClearResult(res)
# head(resout)
# d_backup = d
# d = d_backup
# d_summer = d

#add utc offsets and local time to upload dataset
d_summer = d_summer$data$TimeZone %>%
    mutate(UTCOffset=CurrentOffsetMs / 1000 / 60 / 60) %>%
    select(UTCOffset) %>% bind_cols(d_summer$data$Point)
d_winter = d_winter$data$TimeZone %>%
    mutate(UTCOffset=CurrentOffsetMs / 1000 / 60 / 60) %>%
    select(UTCOffset) %>% bind_cols(d_winter$data$Point)
d = left_join(sites, d_summer, by=c('lat'='Latitude', 'lon'='Longitude'))
d = left_join(d, d_winter, by=c('lat'='Latitude', 'lon'='Longitude'),
    suffix=c('.summer', '.winter'))
old_dst = resout %>% filter(DateTimeUTC < as.POSIXct('2007-03-09 02:00:00'))
new_dst = resout %>% filter(DateTimeUTC > as.POSIXct('2007-03-09 02:00:00'))
###GET DF OF ST EN DST FOR EACH OF THE ABOVE (OR BOTH TOGETHER)
###THEN ASSEMBLE HUGE BOOLEAN EXPRESSION THAT INCORPORATES THEM ALL


# resout %>% filter(LocalDateTime=DateTimeUTC
x = left_join(res)
dl = c('2007-04-01', '2006-04-02', '2005-04-03', '2004-04-04', '2003-04-06', '2002-04-07')
dl2 = c('2007-10-29', '2006-10-29', '2005-10-30', '2004-10-31', '2003-10-26', '2002-10-27')
strftime(as.Date(dl), '%A')
strftime(as.Date(dl2), '%A')

#read in datavalues worksheet as dataframe
dat = read.xlsx('NC_Advanced.xlsx', 'DataValues')
dat = dat[-(1:4),-1]


resout$UTCOffset =
select(resout, value, )
colnames(resout)
