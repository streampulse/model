m18 = read.csv('~/Downloads/SE_M18_2016-10-22.csv', stringsAsFactors=FALSE)
m6 = read.csv('~/Downloads/SE_M6_2016-10-23.csv', stringsAsFactors=FALSE)

m18$solar.time = as.POSIXct(m18$solar.time, tz='CET')

library(streamMetabolizer)

str(m6)

m18
68.3429
18.9539
m6
68.3067
18.9148

m6$solar.time[1]
