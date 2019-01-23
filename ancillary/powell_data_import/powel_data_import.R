# install_github("USGS-R/sbtools")
library(sbtools)
library(sourcetools)
library(stringr)
library(RMariaDB)
library(xml2)
library(rvest)
library(dplyr)
library(tidyr)
library(streamMetabolizer)

setwd('/home/mike/git/streampulse/model/ancillary/powell_data_import')
cur_time = Sys.time()
attr(cur_time, 'tzone') = 'UTC'

#load sb credentials and local directory locations (specific to my machine)
conf = read_lines('/home/mike/git/streampulse/server_copy/sp/config.py')
extract_from_config = function(key){
    ind = which(lapply(conf, function(x) grepl(key, x)) == TRUE)
    val = str_match(conf[ind], '.*\\"(.*)\\"')[2]
    return(val)
}

#log in to sciencebase and get ids for data products
sb_usr = extract_from_config('SB_USER')
sb_pass = extract_from_config('SB_PASS')
authenticate_sb(sb_usr, sb_pass)

data_products = item_list_children('59bff507e4b091459a5e0982',
    fields='id', limit=99999)

#connect to mysql via mariadb
pw = extract_from_config('MYSQL_PW')
con = dbConnect(RMariaDB::MariaDB(), dbname='sp', username='root', password=pw)

#site data ####
dir.create('site_data')

#download and extract site data
site_obj = item_get('59bff64be4b091459a5e098b')
site_files = item_file_download(site_obj, dest_dir='site_data')
tsv = which(grepl('.tsv', site_files))
site_data = read.table(site_files[tsv], header=TRUE,
    sep='\t', stringsAsFactors=FALSE, quote='')

#remove superfluous double quotes
char_cols = lapply(site_data, class) == 'character'
site_data[char_cols] = apply(site_data[,char_cols], 2, str_replace_all,
    pattern='"', replacement='')

#scrape table of state names and abbreviations; turn into mappings df
u = 'http://www.printabledirect.com/list-of-all-50-states-abbreviations-chart.htm'
u_html = read_html(u)
# write_html(u_html, '~/Desktop/hax/utilities/state_abbrevs.html')

elements = html_nodes(u_html, xpath="//tr//font[not(strong)]") %>%
    html_text(trim=TRUE)
elements = elements[elements != '']
elements = str_replace_all(elements, '^([a-zA-Z]*)\\s+?([a-zA-Z]*)$', '\\1\\2')
elements = toupper(elements)
state_abb_map = data.frame(state=elements[seq(1, length(elements), 2)],
    abb=elements[seq(2, length(elements), 2)], stringsAsFactors=FALSE)

#use mappings to determine state abbreviation where possible, document where not
site_data$region = str_match(site_data$X.long_name., ' ([A-Z]{2})$')[,2]
site_data$region[! site_data$region %in% state_abb_map$abb] = NA
full_state_names = str_match(site_data$X.long_name., ',\\s?([a-zA-Z]{3,})$')[,2]
full_state_names = toupper(full_state_names)
abbs_from_full_names = match(full_state_names, state_abb_map$state, nomatch=NA)
fullname_inds = ! is.na(abbs_from_full_names)
site_data$region[which(fullname_inds)] =
    state_abb_map$abb[abbs_from_full_names[fullname_inds]]

still_missing = is.na(site_data$region)
missing_state_latlon = site_data[still_missing, c('X.lat.', 'X.lon.')]

#get the rest from AskGeo, via lat and long (commented to avoid accidental query)
# api_key = extract_from_config('ASKGEO_KEY')
# accnt_id = '2023'
# latlongs = paste(paste(missing_state_latlon$X.lat., missing_state_latlon$X.lon.,
#     sep='%2C'), collapse='%3B')
# askGeoReq = paste0('https://api.askgeo.com/v1/', accnt_id, '/',
#     api_key, '/query.json?databases=UsState2010&points=', latlongs)
#
# r = httr::GET(askGeoReq)
# json = httr::content(r, as="text", encoding="UTF-8")
# askGeoResp = jsonlite::fromJSON(json)

#load saved askgeo results
askGeoResp = readRDS('~/Desktop/askGeoResp.rds')

state_lookups = toupper(askGeoResp$data$UsState2010$CensusAreaName)
state_lookups = str_remove_all(state_lookups, '\\s')
state_lookups_abb = state_abb_map$abb[match(state_lookups, state_abb_map$state)]
state_lookups_abb[is.na(state_lookups_abb)] = c('DC', 'PR', 'PR')

site_data$region[still_missing] = state_lookups_abb

#insert site data into db
nsites = nrow(site_data)
site_data_db = data.frame('region'=site_data$region,
    'site'=site_data$X.site_name.,
    'name'=site_data$X.long_name.,
    'latitude'=site_data$X.lat.,
    'longitude'=site_data$X.lon.,
    'usgs'=site_data$X.nwis_id., 'addDate'=rep(cur_time, nsites),
    'embargo'=rep(0, nsites), 'by'=rep(-902, nsites),
    'contact'=rep('https://doi.org/10.5066/F70864KX', nsites),
    'contactEmail'=rep(NA, nsites),
    stringsAsFactors=FALSE)

dbWriteTable(con, 'site', site_data_db, append=TRUE)

#model inputs ####
dir.create('modIn_data')

#download and extract model input data per site
modIn_children = item_list_children('59eb9b9de4b0026a55ffe37c',
    fields='id', limit=99999)
modIn_ids = sapply(modIn_children, function(x) x$id)

for(i in 1:length(modIn_ids)){

    #download zip, extract, read into R, convert to long format
    modIn_obj = item_get(modIn_ids[i])
    modIn_zip = item_file_download(modIn_obj, dest_dir='modIn_data')
    unzipped = unzip(zipfile=paste0(modIn_zip), exdir=paste0('modIn_data'))
    modIn_data = read.table(unzipped, header=TRUE,
        sep='\t', stringsAsFactors=FALSE, quote='')
    modIn_data = gather(modIn_data, 'variable', 'value', -'solar.time')

    #get site data from above and populate rest of db data table
    sitename = str_match(modIn_zip, 'modIn_data/(nwis_[0-9]+)_.*')[,2]
    site_deets = site_data_db[site_data_db$site == sitename,]

    modIn_data$region = site_deets$region
    modIn_data$site = site_deets$site
    modIn_data$flag = NA
    modIn_data$upload_id = -902
    modIn_data$DateTime_UTC = as.POSIXct(modIn_data$DateTime_UTC,
        format='%Y-%m-%dT%H:%M:%SZ', tz='UTC')
    modIn_data$DateTime_UTC = convert_solartime_to_UTC(modIn_data$DateTime_UTC,
        site_deets$longitude, time.type='mean solar')

    colnames(modIn_data)[colnames(modIn_data) == 'solar.time'] = 'DateTime_UTC'

    dbWriteTable(con, 'data', modIn_data, append=TRUE)
}

#model outputs ####
dir.create('modOut')

#download and extract model input data per site
modOut_children = item_list_children('59eb9ba6e4b0026a55ffe37f',
    fields='id', limit=99999)
modOut_ids = sapply(modOut_children, function(x) x$id)

for(i in 1:length(modOut_ids)){

    modOut_obj = item_get(modOut_ids[i])
    modOut_zip = item_file_download(modOut_obj, dest_dir='modOut')
    unzipped = unzip(zipfile=paste0(modOut_zip), exdir=paste0('modOut'))
    items = str_match(unzipped, '^modOut/(.*).tsv$')[,2]
    modOut = list()
    invisible(lapply(items, function(x) {
        modOut[[x]] <<- read.table(paste0('modOut/', x, '.tsv'), header=TRUE,
            sep='\t', stringsAsFactors=FALSE, quote='')
        # assign(x, read.table(paste0('modOut/', x, '.tsv'), header=TRUE,
        #     sep='\t', stringsAsFactors=FALSE, quote=''), pos=.GlobalEnv)
    }))
    modOut = gather(modOut, 'variable', 'value', -'solar.time')

    #get site data from above and populate rest of db data table
    sitename = str_match(modOut_zip, 'modOut_data/(nwis_[0-9]+)_.*')[,2]
    site_deets = site_data_db[site_data_db$site == sitename,]

}
az = readRDS('~/git/streampulse/server_copy/sp/shiny/data/modOut_AZ_LV_2018.rds')
azp = readRDS('~/git/streampulse/server_copy/sp/shiny/data/predictions_AZ_LV_2018.rds')
str(az$data_daily)
str(az$data)
str(az$fit$daily)

str(modOut$daily)
str(modOut$KQ_overall)
str(modOut$overall)
str(modOut$KQ_binned)
