library(dplyr)

site_deets = read.csv('~/git/streampulse/model/site_deets.csv',
    stringsAsFactors=FALSE)
site_deets = filter(site_deets, substr(site_code, 1, 2) == 'NC') %>%
    filter(! duplicated(site_code))

V = site_deets$velocity_ms
S = site_deets$slope_prop
D = site_deets$depth_m
Q = site_deets$discharge_m3s

#get min and max for all terms in raymond eqn 7
t1a = (V * S)^(0.86 + 0.016)
t1b = (V * S)^(0.86 - 0.016)
t2a = Q^(-0.14 + 0.012)
t2b = Q^(-0.14 - 0.012)
t3a = D^(0.66 + 0.029)
t3b = D^(0.66 - 0.029)
t1min = apply(cbind(t1a, t1b), 1, min)
t1max = apply(cbind(t1a, t1b), 1, max)
t2min = apply(cbind(t2a, t2b), 1, min)
t2max = apply(cbind(t2a, t2b), 1, max)
t3min = apply(cbind(t3a, t3b), 1, min)
t3max = apply(cbind(t3a, t3b), 1, max)

#calculate k600 (m/d)
site_deets$raymondk600min = (4725 - 445) * t1min * t2min * t3min
site_deets$raymondk600max = (4725 + 445) * t1max * t2max * t3max
site_deets$raymondk600mean = 4725 * (V * S)^(0.86) * Q^(-0.14) * D^(0.66)

#convert to K600 (1/d)
K600s = data.frame(
    apply(select(site_deets, starts_with('raymond')), 2, function(x) x / D))
colnames(K600s) = sub('k', 'K', colnames(K600s))
site_deets = cbind(site_deets, K600s)

#visualize
png(width=6, height=7, units='in', type='cairo', res=300,
    filename='~/Dropbox/streampulse/figs/raymond_K600.png')
par(mfrow=c(2, 1), mar=c(0, 4, 4, 1))
plot(1:7, site_deets$raymondk600mean, xaxt='n', xlab='', ylab='k600 (m/d)',
    pch=20, ylim=range(c(site_deets$raymondk600min, site_deets$raymondk600max)))
segments(1:7, site_deets$raymondk600min, 1:7, site_deets$raymondk600max)
mtext('Raymond eqn 7 gas exchange', 3, line=0)

par(mar=c(4, 4, 0, 1))
plot(1:7, site_deets$raymondK600mean, xaxt='n', xlab='', ylab='K600 (1/d)',
    pch=20, ylim=range(c(site_deets$raymondK600min, site_deets$raymondK600max)))
axis(1, substr(site_deets$site_code, 4, nchar(site_deets$site_code)), at=1:7)
segments(1:7, site_deets$raymondK600min, 1:7, site_deets$raymondK600max)
dev.off()

#save data
savedata = select(site_deets, -start_date, -end_date, -tok, -int, -skip)
write.csv(savedata, '~/Dropbox/streampulse/data/raymond_K600.csv')
