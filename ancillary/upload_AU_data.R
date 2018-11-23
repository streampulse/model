library(RMariaDB)
library(DBI)
library(stringr)
library(dplyr)

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

#setup
setwd('/home/mike/Dropbox/streampulse/data/australia_data/')

# sitenm_dict = list('Edward'='EW', 'Goulbu'='Goulburn',
#     'Lachla'='Lachlan')

dirs = list.dirs(recursive=FALSE)
d=dirs[1]; f=files[1]
for(d in dirs){
    files = list.files(d)
    for(f in files){
        data = read.csv(paste0(d, '/', f), stringsAsFactors=FALSE)
        data$DateTime_UTC = as.POSIXct(paste(data$Date, data$Time),
            format='%d/%m/%Y %H:%M:%S', tz='UTC')
        data = data %>% select(-Date, -Time, -salinity) %>%
            select(DateTime_UTC, everything())

        cur_time = Sys.time()
        attr(cur_time, 'tzone') = 'UTC'
        site_deets = stringr::str_match(f, '(.+)_(\\d+)_(.+)_(\\d+).csv')
        sql_up_df = data.frame('region'='AU',
            'site'=paste0(site_deets[4], site_deets[3]),
            # 'site'=sitenm_dict[substr(f, 1, 6)][[1]],
            'name'=paste0(site_deets[2], ': ', site_deets[4], ' ', site_deets[3])
            'latitude'=,
            'longitude'=,
            'usgs'=NA, 'addDate'=cur_time,
            'embargo'=0, 'by'=35,
            'contact'=, 'contactEmail'=)

        dbWriteTable(con, 'data', na_filt, append=TRUE)
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

