#script for generating gradient boosted classification tree model to predict
#problematic datapoints in data uploaded to StreamPULSE data portal

#acquire training and testing data####

# install.packages('RMySQL')
library(RMySQL)

# MySQL() #should do the same thing as dbDriver
drv <- dbDriver('MySQL')

#first try it the usual way
# sp = dbConnect(drv, user='root', password='',
#     dbname='sp', host='45.55.47.104')#, port=3306)

#if that fails, try port forwarding (connect to local 3307 as if it's remote 3306).
#in terminal, forward local port 3307 to remote port 3306
#ssh -L 3307:45.55.47.104:3306 aaron@45.55.47.104
# sp = dbConnect(drv, user='root', password='',
#     dbname='sp', host='localhost', port=3307)

#if that fails, remote connection is probably not configured, server side.
#next option is to sftp a zipped database dump and reconstitute the db locally.
pw = Sys.getenv('SPDB') #rstudio must be spawned from terminal
sp = dbConnect(drv, user='root', password=pw,
    dbname='sp', host='localhost')

#roll time
qaqc_sites = dbGetQuery(sp, 'select distinct site from flag;')
qaqc_sites = paste0('"',
    paste(as.vector(as.matrix(qaqc_sites)), collapse='","'), '"')
q = dbGetQuery(sp, paste('select * from data where site in (',
    qaqc_sites, ');'))
# saveRDS(q, '/home/mike/Desktop/untracked/anomaly_set.rds')

#for changing the db itself (i.e. not just select statements)
# q = dbSendQuery(...)
# data = fetch(q, n=-1) #n=-1 fetches all pending records

#engineer features ####

#setup
# install.packages("drat", repos="https://cran.rstudio.com")
# drat:::addRepo("dmlc")
# install.packages("xgboost", repos="http://dmlc.ml/drat/", type = "source")
library(xgboost)
library(tidyr)
library(data.table)

d = readRDS('/home/mike/Desktop/untracked/anomaly_set.rds')
d = spread(d, variable, value)
m = gather(d,
head(m)
# d$DateTimeUTC = as.POSIXct(d$DateTimeUTC)#, format='%Y-%m-%d %H:%M:%S')
# tm = ts(x$DO_mg.l, deltat = 1/288, start=188)
# tm = na.seadec(tm, algorithm='interpolation')

#throw away site WS1500 for now (not enough flags for ML)
d = d[!d$site == 'WS1500',]

#create labels
labels = d$flag
labels[!is.na(labels)] = 1 #flagged
labels[is.na(labels)] = 0 #not flagged

#get some deets on flaggage per site
sites = unique(d$site)
for(i in sites){
    print(i) #site name
    print(sum(d$site == i)) #number of obs for site
    print(sum(labels & d$site == i)) #number of flags for site
}

#create lags and leads of target variable. add them to feature matrix.
features = matrix(NA, nrow=length(labels), ncol=21+length(sites))
features[,1] = scale(d$DO_mgL)
lags = shift(features[,1], n=1:10, fill=NA, type='lag', give.names=TRUE)
features[,2:11] = do.call(cbind, lags)
leads = shift(features[,1], n=1:10, fill=NA, type='lead', give.names=TRUE)
features[,12:21] = do.call(cbind, leads)

#convert site column to one-hot matrix
library('caret')
features[,22:ncol(features)] = class2ind(factor(d$site))

#split training and testing sets
intrain = createDataPartition(labels, p=0.75, list=FALSE)
features_train = features[intrain,]
features_test = features[-intrain,]
labels_train = labels[intrain]
labels_test = labels[-intrain]

#run model

#input data as dgCMatrix (sparse)
# bstSparse = xgboost(data=train$data, label=train$label, max_depth=2, eta=1,
#     nthread=2, nrounds=2, objective="binary:logistic")

#input as base (dense) matrix
bst = xgboost(data=features_train, label=labels_train, max_depth=2,
    eta=1, nthread=2, nrounds=2, objective="binary:logistic")

#input as some kind of turbo matrix that unlocks sweet xgboost hax if ur legit
# dtrain = xgb.DMatrix(data=train$data, label=train$label)
# bstDMatrix = xgboost(data=dtrain, max_depth=2, eta=1, nthread=2,
#     nrounds=2, objective="binary:logistic", verbose=2)

#predict and convert to binary classification
pred = predict(bst, features_test)
prediction = as.numeric(pred > 0.5)

#average error
err = mean(as.numeric(pred > 0.5) != test$label)
print(paste("test-error=", err))

probabilityVectorPreviouslyComputed != test$label
mean(vectorOfErrors)
