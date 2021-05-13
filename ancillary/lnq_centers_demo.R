#make some fake datetime and discharge data
dt = seq(as.POSIXct('2017-01-01 00:00:00'), as.POSIXct('2017-12-31 00:00:00'),
    '15 min')
Q = rnorm(length(dt), 30, 10)
Q[Q <= 0] = 0.001
fitdata = data.frame(solar.time=dt, discharge=Q)

#determine centers (same way it's done in the package)
n_KQ_nodes = 7
addis = tapply(log(fitdata$discharge),
    substr(fitdata$solar.time,1,10), mean)
K600_lnQ_nodes_centers = seq(from=min(addis, na.rm=TRUE),
    to=max(addis, na.rm=TRUE), length.out=n_KQ_nodes)

#verify
plot(as.Date(names(addis)), addis)
abline(h=K600_lnQ_nodes_centers, col='red')
