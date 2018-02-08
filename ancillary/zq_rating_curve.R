d = streampulse_data$data
wp = d[d$variable == 'WaterPres_kPa', c('DateTime_UTC','value')]
ap = d[d$variable == 'AirPres_kPa', c('DateTime_UTC','value')]
plot(wp, type='l')
plot(ap, type='l')
dim(wp)
dim(ap)
library(dplyr)
p = full_join(wp, ap, by='DateTime_UTC')
dim(p)
colnames(p)[2:3] = c('wat','air')
plot(p$air, type='l')
p$w_minus_a = p$wat - p$air
