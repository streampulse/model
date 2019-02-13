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
par(mfrow=c(2,1), mar=c(0,0,0,0))
plot(p$air, type='l')
plot(p$wat, type='l')

#remove air pressure so all pressure is from hydraulic head
p$w_minus_a = p$wat - p$air
#compute fluid density
p$p =  (999.83952 + 16.945176 p$ - 7.9870401e-03 T ref 2 - 46.170461e-06 T
    ref 3 +  105.56302e-09 T ref 4  - 280.54253e-12 T ref 5 ) / (1 + 16.879850e-03 T
    ref )

#2####

#remove air pressure so all pressure is from hydraulic head ("downwell pressure")
dw_pres = dd$WaterPres_kPa - dd$AirPres_kPa

#compute fluid density
wat_temp = dd$WaterTemp_C
T1 = 16.945176 * wat_temp
T1b = 16.879850e-03 * wat_temp
T2 = 7.9870401e-03 * wat_temp^2
T3 = 46.170461e-06 * wat_temp^3
T4 = 105.56302e-09 * wat_temp^4
T5 = 280.54253e-12 * wat_temp^5
fl_dens = (999.83952 + T1 - T2 - T3 + T4 - T5) / (1 + T1b)

#convert fluid density to lb/ft^3
fl_dens = 0.0624279606 * fl_dens
# mean(fl_dens[!is.na(fl_dens)]) #bogus

#3####

#convert downwell pressure to density dependent fluid depth
ft_to_m = 0.3048
kPa_to_psi = 0.1450377
psi_to_psf = 144.0
fl_dens = 62.408
fl_depth = ft_to_m * (kPa_to_psi * psi_to_psf * dw_pres) / fl_dens

fl_depth = sens_depth
mean(sens_depth[!is.na(sens_depth)])
plot(density((sens_depth[!is.na(sens_depth)])))

ref_lvl - (dw_pres - fl_depth)


#ugh. no reference time. let's just not use a ref water level at all ####

#compute hydraulic pressure (downwell pressure minus air pressure)
dd = dd2
hyd_pres = dd$WaterPres_kPa - dd$AirPres_kPa

#convert hydraulic pressure to sensor depth
ft_to_m = 0.3048
kPa_to_psi = 0.1450377
psi_to_psf = 144.0

fl_dens = 62.408 #lb/ft^3 at 50F (10C)
sens_depth = ft_to_m * (kPa_to_psi * psi_to_psf * hyd_pres) / fl_dens
mean(sens_depth[!is.na(sens_depth)])
plot(density((sens_depth[!is.na(sens_depth)])))

sens_height = .21 #in m
wat_lvl = sens_depth + sens_height
mean(wat_lvl[!is.na(wat_lvl)])
max(wat_lvl[!is.na(wat_lvl)])
plot(density((wat_lvl[!is.na(wat_lvl)])))

#convert lvl to meters

#run level through site-specific rating curve to get discharge
disch = 0.3610671493 * wat_lvl^8.8560174769
plot(wat_lvl, disch, ylim=c(0, 100))

zq = read.csv('~/Dropbox/streampulse/data/ZQ_data.csv')
zq = zq[order(z),]
z = zq$level_m
q = zq$discharge_cms
mod = nls(q ~ (a * z^b), start=list(a=0.1, b=1))
plot(z, q)
lines(z, predict(mod, list(z=sort(z))))
# class(mod$m$getPars())
params = summary(mod)$parameters
params[1,1]
params[2,1]
