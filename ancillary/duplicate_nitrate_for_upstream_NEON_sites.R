library(RMariaDB)
library(DBI)

# setwd('/home/mike/git/streampulse/server_copy/sp/scheduled_scripts/')
setwd('/home/aaron/sp/scheduled_scripts/')

conf = readLines('../config.py')
extract_from_config = function(key){
    ind = which(lapply(conf, function(x) grepl(key, x)) == TRUE)
    val = str_match(conf[ind], '.*\\"(.*)\\"')[2]
    return(val)
}
pw = extract_from_config('MYSQL_PW')

con = dbConnect(RMariaDB::MariaDB(), dbname='sp',
    username='root', password=pw)

r = dbSendQuery(con, paste("SELECT *",
    "FROM data WHERE upload_id=-900 and variable='Nitrate_mgL'"))
res = dbFetch(r)
dbClearResult(r)

res$site = paste0(substr(res$site,1,4), '-up')
res$id = NULL

dbWriteTable(con, 'data', res, append=TRUE)

dbDisconnect(con)
