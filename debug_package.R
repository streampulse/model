detach('package:StreamPULSE', unload=TRUE)
setwd('~/git/streampulse/model/StreamPULSE/R/')
fs = list.files()
for(i in fs) source(i)
