# install_github("USGS-R/sbtools")
library(sbtools)
library(sourcetools)
library(stringr)
library(RMariaDB)
library(rvest)
library(xml2)
library(dplyr)

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

#download and extract site data
site_obj = item_get('59bff64be4b091459a5e098b')
site_files = item_file_download(site_obj, dest_dir='.')
tsv = which(grepl('.tsv', site_files))
site_data = read.table(site_files[tsv], header=TRUE, sep='\t',
    stringsAsFactors=FALSE, quote='')

#remove superfluous double quotes
char_cols = lapply(site_data, class) == 'character'
site_data[char_cols] = apply(site_data[,char_cols], 2, str_replace_all,
    pattern='"', replacement='')

#scrape table of state names and abbreviations; turn into mappings df
u = 'http://www.printabledirect.com/list-of-all-50-states-abbreviations-chart.htm'
u_html = read_html(u)

html_node("Verdana, Arial, Helvetica, sans-serif") %>%
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

#get the rest from AskGeo, via lat and long
api_key = extract_from_config('ASKGEO_KEY')
accnt_id = '2023'
latlongs = paste(paste(missing_state_latlon$X.lat., missing_state_latlon$X.lon.,
    sep='%2C'), collapse='%3B')
askGeoReq = paste0('https://api.askgeo.com/v1/', accnt_id, '/',
    api_key, '/query.json?databases=UsState2010&points=', latlongs)

r = httr::GET(askGeoReq)
json = httr::content(r, as="text", encoding="UTF-8")
askGeoResp = jsonlite::fromJSON(json)
# saveRDS(askGeoResp, '~/Desktop/askGeoResp.rds')

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
    'contactEmail'=rep(NA, nsites))

dbWriteTable(con, 'site', site_data_db, append=TRUE)


