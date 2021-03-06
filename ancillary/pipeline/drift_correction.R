x = read.csv('~/Desktop/NC_NHC_sensorData.csv', stringsAsFactors = F)
plot(x$DO_mgL, type='l', xlim=c(900,1500), ylim=c(4,10))

install.packages('driftR')
library(driftR)
data(sondeRaw)
data(sondeCal)
data(sondeClean)
str(sondeRaw)
plot(sondeRaw$DO, type='l')
lines(sondeClean$DO_Corr, col='red')
