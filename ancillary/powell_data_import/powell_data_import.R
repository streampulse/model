    #FILTER SITE DATA SO THAT IT ONLY INCLUDES SITES THAT WERE MODELED

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

setwd('/home/aaron/mike_files/powell_data_import')
# setwd('/home/mike/git/streampulse/model/ancillary/powell_data_import')
cur_time = Sys.time()
attr(cur_time, 'tzone') = 'UTC'

#load sb credentials and local directory locations (specific to my machine)
# conf = read_lines('/home/mike/git/streampulse/server_copy/sp/config.py')
conf = read_lines('/home/aaron/sp/config.py')
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
site_files = item_file_download(site_obj, dest_dir='site_data',
    overwrite_file=TRUE)
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
# read_html('~/Desktop/hax/utilities/state_abbrevs.html')

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
api_key = extract_from_config('ASKGEO_KEY')
accnt_id = '2023'
latlongs = paste(paste(missing_state_latlon$X.lat., missing_state_latlon$X.lon.,
    sep='%2C'), collapse='%3B')
askGeoReq = paste0('https://api.askgeo.com/v1/', accnt_id, '/',
    api_key, '/query.json?databases=UsState2010&points=', latlongs)

r = httr::GET(askGeoReq)
json = httr::content(r, as="text", encoding="UTF-8")
askGeoResp = jsonlite::fromJSON(json)

#load saved askgeo results
askGeoResp = readRDS('~/Desktop/askGeoResp.rds')

state_lookups = toupper(askGeoResp$data$UsState2010$CensusAreaName)
state_lookups = str_remove_all(state_lookups, '\\s')
state_lookups_abb = state_abb_map$abb[match(state_lookups, state_abb_map$state)]
state_lookups_abb[is.na(state_lookups_abb)] = c('DC', 'PR', 'PR')

site_data$region[still_missing] = state_lookups_abb

#modify sitenames for portal compatibility
site_data$X.site_name. = sub('_', '-', site_data$X.site_name.)

#prepare site data for db insertion (later)
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

#model inputs ####
dir.create('modIn_data')
dir.create('RDS_components')
dir.create('RDS_components/outData')
dir.create('RDS_components/inData')
dir.create('RDS_components/outExtra')

#download and extract model input data per site
modIn_children = item_list_children('59eb9b9de4b0026a55ffe37c',
    fields='id', limit=99999)
modIn_ids = sapply(modIn_children, function(x) x$id)

for(i in 1:length(modIn_ids)){

    #download zip, extract, read in
    modIn_obj = item_get(modIn_ids[i])
    modIn_zip = item_file_download(modIn_obj, dest_dir='modIn_data',
        overwrite_file=TRUE)
    unzipped = unzip(zipfile=paste0(modIn_zip), exdir=paste0('modIn_data'))
    modIn_data = read.table(unzipped, header=TRUE,
        sep='\t', stringsAsFactors=FALSE, quote='')

    #get site data from above
    sitename = str_match(modIn_zip, 'modIn_data/(nwis_[0-9]+)_.*')[,2]
    sitename = sub('_', '-', sitename)
    site_deets = site_data_db[site_data_db$site == sitename,]

    #save df that will become the "data" component of the RDS obj saved on server
    modIn_data$date = as.Date(modIn_data$solar.time, format='%Y-%m-%dT%H:%M:%SZ')
    modIn_data$solar.time = as.POSIXct(modIn_data$solar.time,
        format='%Y-%m-%dT%H:%M:%SZ', tz='UTC')
    sitecode = str_match(sitename, '[0-9]+')[[1]]
    modyear = substr(modIn_data$date[nrow(modIn_data) / 2], 1, 4)

    #save input object
    saveRDS(modIn_data, paste0('RDS_components/inData_', site_deets$region, '_',
        sitecode, '_', modyear, '.rds'))

    modIn_data$date = NULL

    #convert to long format, populate rest of db data table
    modIn_data$DateTime_UTC = convert_solartime_to_UTC(modIn_data$solar.time,
        site_deets$longitude, time.type='mean solar')
    modIn_data$solar.time = NULL
    colname_map = c('depth'='Depth_m', 'discharge'='Discharge_m3s',
        'DO.obs'='DO_mgL', 'DO.sat'='satDO_mgL', 'light'='Light_PAR',
        'temp.water'='WaterTemp_C', 'DateTime_UTC'='DateTime_UTC')
    newnames = unname(colname_map[match(colnames(modIn_data), names(colname_map))])
    colnames(modIn_data) = newnames
    # colnamefac = factor(colnames(modIn_data))
    # levels(colnamefac)

    modIn_data = gather(modIn_data, 'variable', 'value', -'DateTime_UTC')

    modIn_data$region = site_deets$region
    modIn_data$site = site_deets$site
    modIn_data$flag = NA
    modIn_data$upload_id = -902

    dbWriteTable(con, 'powell', modIn_data, append=TRUE)

    # spread(modIn_data[1:10,], 'variable', 'value') %>%
    #     select(DateTime_UTC
}


# temp: read modIn data, split by year, rewrite ####
#didn't do this above, during acquisition from sciencebase, so doing now
rdsdir = '/home/mike/git/streampulse/model/ancillary/powell_data_import/RDS_components/temp'
rr = dir(rdsdir, pattern='inData')
for(r in rr){
    print(r)
    f = readRDS(paste0(rdsdir, '/', r))
    datebounds = range(f$date)
    datebounds = as.character(datebounds)
    if(substr(datebounds[1], 6, 7) == '12'){
        datebounds[1] = as.numeric(substr(datebounds[1], 1, 4)) + 1
    }
    if(substr(datebounds[2], 6, 7) == '01'){
        datebounds[2] = as.numeric(substr(datebounds[2], 1, 4)) - 1
    }
    datebounds = as.numeric(substr(datebounds, 1, 4))
    yrs = seq(datebounds[1], datebounds[2], 1)
    for(y in yrs){
        ff = f[substr(f$date, 1, 4) == y,]
        if(nrow(ff) > 0){
            # print(paste(y, '-', nrow(ff)))
            substr(r, nchar(r) - 7, nchar(r)) = paste0(y, '.rds')
            saveRDS(ff, paste0('RDS_components/inData/', r))
            # substr(r, nchar(r) - 7, nchar(r)) = paste0(y, '.csv')
            # write.csv(ff, paste0('RDS_components/', r))
        }
    }
}

#model outputs ####
dir.create('modOut')

#download and extract model output series per site; save to disk for later
modOut_children = item_list_children('59eb9ba6e4b0026a55ffe37f',
    fields='id', limit=99999)
modOut_ids = sapply(modOut_children, function(x) x$id)

for(i in 1:length(modOut_ids)){

    #acquire object details from sciencebase, reauthenticate if necessary
    modOut_obj = try(item_get(modOut_ids[i]))
    if(class(modOut_obj) == 'try-error'){
        print('reauthenticating')
        authenticate_sb(sb_usr, sb_pass)
        modOut_obj = item_get(modOut_ids[i])
    }

    #download, unzip, identify components
    modOut_zip = item_file_download(modOut_obj, dest_dir='modOut',
        overwrite_file=TRUE)

    if(length(modOut_zip) > 1){
        finer_res = which.min(as.numeric(str_match(modOut_zip,
            '([0-9]{2,3})min')[,2]))
        modOut_zip = modOut_zip[finer_res]
    }

    sitename = str_match(modOut_zip, 'modOut/(nwis_[0-9]+)?_.*')[,2]
    try(dir.create(paste0('modOut/', sitename)))
    unzipped = unzip(zipfile=paste0(modOut_zip),
        exdir=paste0('modOut/', sitename))
    items = str_match(unzipped, '^modOut/nwis_[0-9]+/(.*).tsv$')[,2]
    items = na.omit(items)

    #get site details for savefile names
    site_deets = site_data_db[site_data_db$site == sitename,]
    sitecode = str_match(sitename, '[0-9]+')[[1]]

    #some components are cumulative; read these tsvs into a list and save
    items = items[items != 'daily']
    modOut = list()
    invisible(lapply(items, function(x) {
        modOut[[x]] <<- read.table(paste0('modOut/', sitename, '/', x, '.tsv'),
            header=TRUE, sep='\t', stringsAsFactors=FALSE, quote='')
    }))

    #write output extras object
    saveRDS(modOut, paste0('RDS_components/outExtra/outExtra_',
        site_deets$region, '_', sitecode, '.rds'))

    #read in core output data, separate by year, save
    outdf = read.table(paste0('modOut/', sitename, '/daily.tsv'),
        header=TRUE, sep='\t', stringsAsFactors=FALSE, quote='')

    datebounds = range(outdf$date)
    datebounds = as.character(datebounds)
    if(substr(datebounds[1], 6, 7) == '12'){
        datebounds[1] = as.numeric(substr(datebounds[1], 1, 4)) + 1
    }
    if(substr(datebounds[2], 6, 7) == '01'){
        datebounds[2] = as.numeric(substr(datebounds[2], 1, 4)) - 1
    }
    datebounds = as.numeric(substr(datebounds, 1, 4))
    yrs = seq(datebounds[1], datebounds[2], 1)

    for(y in yrs){
        ff = outdf[substr(outdf$date, 1, 4) == y,]
        if(nrow(ff) > 0){

            #write output object
            saveRDS(ff, paste0('RDS_components/outData/outData_',
                site_deets$region, '_', sitecode, '_', y, '.rds'))
        }
    }

    print(paste(sitecode, '-', length(yrs)))

    # sitecode = str_match(sitename, '[0-9]+')[[1]]
    # modyear = substr(modIn_data$date[nrow(modIn_data) / 2], 1, 4)
    # saveRDS(modOut, paste0('RDS_components/', sitename, '_modOut.rds'))
}


#model estimates and predictors ####
dir.create('estpred')

#download and extract estpred data
estpred_obj = item_get('59eb9c0ae4b0026a55ffe389')
estpred_files = item_file_download(estpred_obj, dest_dir='estpred',
    overwrite_file=TRUE)
zp = which(grepl('.zip', estpred_files))
unzipped = unzip(zipfile=paste0(estpred_files[zp]), exdir='estpred')
estpred = read.table(unzipped, header=TRUE,
    sep='\t', stringsAsFactors=FALSE, quote='')

ep_slice = estpred[estpred$site_name == 'nwis_06461500',]

saveRDS(estpred, paste0('RDS_components/estpred.rds'))

#model config ####
dir.create('config')

#download and extract config data
config_obj = item_get('59eb9b86e4b0026a55ffe379')
config_files = item_file_download(config_obj, dest_dir='config',
    overwrite_file=TRUE)
zp = which(grepl('.zip', config_files))
unzipped = unzip(zipfile=paste0(config_files[zp]), exdir='config')
config = read.table(unzipped, header=TRUE,
    sep='\t', stringsAsFactors=FALSE, quote='')

saveRDS(config, paste0('RDS_components/config.rds'))

#model diagnostics ####
dir.create('diag')

#download and extract diag data
diag_obj = item_get('59eb9bafe4b0026a55ffe382')
diag_files = item_file_download(diag_obj, dest_dir='diag',
    overwrite_file=TRUE)
zp = which(grepl('.zip', diag_files))
unzipped = unzip(zipfile=paste0(diag_files[zp]), exdir='diag')
mod_diag = read.table(unzipped, header=TRUE,
    sep='\t', stringsAsFactors=FALSE, quote='')

str(mod_diag)
saveRDS(mod_diag, paste0('RDS_components/diag.rds'))


#other ####

#filter sites so that only those that were modeled get inserted into database
outd_files = dir('RDS_components/outData/', pattern='outData')
sitecodes = str_match(outd_files, '([0-9]+)_[0-9]{4}.rds$')[,2]
all_sitecodes = str_match(site_data_db$site, '[0-9]+')
unmodeled = which(! all_sitecodes %in% sitecodes)
site_data_db = site_data_db[-unmodeled,]

#insert
dbWriteTable(con, 'site', site_data_db, append=TRUE)


#testing####
# testpath2 = '/home/mike/git/streampulse/model/ancillary/powell_data_import/RDS_components/nwis_06461500_inData.rds'
# rdsdata = readRDS(testpath2)
#
# str(rdsdata)
#
# testfit = readRDS('~/Desktop/test_modelfit.rds')
# az = readRDS('~/git/streampulse/server_copy/sp/shiny/data/modOut_AZ_LV_2018.rds')
# azp = readRDS('~/git/streampulse/server_copy/sp/shiny/data/predictions_AZ_LV_2018.rds')
# str(az$data_daily) #missing, but not needed?
# str(az$data) #all but DO.mod
# str(az$fit$daily)
#
# str(modOut$daily)
# str(modOut$KQ_overall)
# str(modOut$overall)
# str(modOut$KQ_binned)
#
# rdsdata = readRDS(paste0('RDS_components/', sitename, '_inData.rds'))
# str(rdsdata)
# testpath = '/home/mike/git/streampulse/model/ancillary/powell_data_import/RDS_components/nwis_11463000_modOut.rds'
#
#
#
#
# authenticate_sb(sb_usr, sb_pass)
# data_products
