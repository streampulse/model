#sudo apt install libgsl0-dev
# install.packages('tsoutliers')
library(tsoutliers)
library(imputeTS)
list.files('~/Dropbox/streampulse/data')
#v long
x = read.csv('/home/mike/Dropbox/streampulse/data/WI_BEC_2012-02-16_XX.csv')
#short
x = read.csv('/home/mike/Dropbox/streampulse/data/NC_UEno_2017-01-04_HL_test.csv')
z = x$Abs.Pres..kPa..LGR.S.N..10838018..SEN.S.N..10838018.
z[is.na(z)] = 10
#med
x = read.csv('/home/mike/git/streampulse/server_copy/spuploads/WI_BEC_2014-09-24_XX.csv',
    stringsAsFactors=FALSE)
head(x)
str(x)
x$dateTimeUTC = as.POSIXct(x$dateTimeUTC, format='%Y-%m-%d %H:%M:%S')
tm = ts(x$DO_mg.l, deltat = 1/288, start=188)
tm = na.seadec(tm, algorithm='interpolation')
z = data.frame(date=x$dateTimeUTC, y=tm)
head(z)

out = tso(tm, types='AO', tsmethod='stsm')
plot(out)
# print(Sys.time())


#gangster mode
tsoutl <- function(x,plot=FALSE)
{
    x <- as.ts(x)
    if(frequency(x)>1)
        resid <- stl(x,s.window="periodic",robust=TRUE)$time.series[,3]
    else
    {
        tt <- 1:length(x)
        resid <- residuals(loess(x ~ tt))
    }
    resid.q <- quantile(resid,prob=c(0.25,0.75))
    iqr <- diff(resid.q)
    limits <- resid.q + 1.5*iqr*c(-1,1)
    score <- abs(pmin((resid-limits[1])/iqr,0) + pmax((resid - limits[2])/iqr,0))
    if(plot)
    {
        plot(x)
        x2 <- ts(rep(NA,length(x)))
        x2[score>0] <- x[score>0]
        tsp(x2) <- tsp(x)
        points(x2,pch=19,col="red", cex=0.5)
        return(invisible(score))
    }
    else
        return(score)
}
x$DO_mg.l[is.na(x$DO_mg.l)] = 10
tsoutl(ts(x$DO_mg.l, deltat=1/288, start=188), TRUE) #for short
    # deltat=1/96)[-(771:772)], TRUE) #for short
tsoutl(ts(x$Abs.Pres..kPa..LGR.S.N..10838018..SEN.S.N..10838018.), TRUE)

#twitter
# devtools::install_github("twitter/AnomalyDetection")
library(AnomalyDetection)
# help(AnomalyDetectionTs) #for ts obj
# help(AnomalyDetectionVec) #for vector
res = AnomalyDetectionTs(z, max_anoms=0.02, direction='both', plot=TRUE,
    longterm=TRUE)
res = AnomalyDetectionVec(z, max_anoms=0.02, period=96,
    direction='both', plot=TRUE)
res$plot

#holt-winters
A = B = Q = seq(0.1, 0.9, 0.1)
best_A = best_B = best_Q = NA
# err = matrix(NA, nrow=length(A), ncol=4)
err = Inf
for(a in A){
    for(b in B){
        for(q in Q){
            hw = HoltWinters(tm, alpha=a, beta=b, gamma=q)
            print(paste(hw$alpha, hw$beta, hw$gamma, hw$SSE))
            if(!is.nan(hw$SSE) & hw$SSE < err){
                err = hw$SSE
                best_A = hw$alpha; best_B = hw$beta; best_Q = hw$gamma
            }
        }
    }
}
best_A; best_B; best_Q
hw = HoltWinters(tm, alpha=best_A, beta=best_B, gamma=best_Q)
x11()
plot(hw$fitted[,1], col='red', xlim=c(200,300), ylim=c(10,11))
par(new=T)
plot(hw$x, xlim=c(200,300), ylim=c(10,11))

#anomalous (identify anomalous series, rather than points)
devtools::install_github("robjhyndman/anomalous")
z <- ts(matrix(rnorm(3000),ncol=100),freq=4)
y <- tsmeasures(z)
biplot.features(y)
anomaly(y)


#take two (four-week segments) ####
setwd('/home/mike/Dropbox/streampulse/data/NC_download/')
site_files = dir()
d = read.csv(site_files[1], stringsAsFactors=FALSE,
    colClasses=c('DateTime_UTC'='POSIXct'))
sitevars = unique(d$variable)
d = d[d$variable == sitevars[1],]
head(d)
d = d[d$DateTime_UTC > as.POSIXct('2016-10-01 00:00:00') & d$DateTime_UTC <
        as.POSIXct('2016-10-29 00:00:00'),]
plot(d$value)

#twitter
# series_in = d[,c('DateTime_UTC', 'value')]
# devtools::install_github("twitter/AnomalyDetection")
library(AnomalyDetection)
# help(AnomalyDetectionTs) #for ts obj
# help(AnomalyDetectionVec) #for vector
res = AnomalyDetectionVec(d$value, max_anoms=0.02, direction='both',
    plot=TRUE, period=96)
res$plot

#tsoutliers
library(tsoutliers)
library(imputeTS)
tm = ts(d$value, deltat = 1/96)
tm = na.seadec(tm, algorithm='interpolation')
# z = data.frame(date=d$DateTime_UTC, y=tm)

#omg so slow
out = tso(tm, types='AO', tsmethod='auto.arima')
plot(out)


#hot-winters
A = B = Q = seq(0.1, 0.9, 0.1)
best_A = best_B = best_Q = NA
# err = matrix(NA, nrow=length(A), ncol=4)
err = Inf
for(a in A){
    for(b in B){
        for(q in Q){
            hw = HoltWinters(tm, alpha=a, beta=b, gamma=q)
            print(paste(hw$alpha, hw$beta, hw$gamma, hw$SSE))
            if(!is.nan(hw$SSE) & hw$SSE < err){
                err = hw$SSE
                best_A = hw$alpha; best_B = hw$beta; best_Q = hw$gamma
            }
        }
    }
}
best_A; best_B; best_Q
hw = HoltWinters(tm, alpha=best_A, beta=best_B, gamma=best_Q)
x11()
plot(hw$fitted[,1], col='red', xlim=c(0,30))#, ylim=c(10,11))
par(new=T)
plot(hw$x, xlim=c(0,30))#, ylim=c(10,11))


#gangster mode
tsoutl = function(x, plot=FALSE)
{
    # x <- as.ts(x)
    # if(frequency(x)>1)
    #     resid <- stl(x,s.window="periodic",robust=TRUE)$time.series[,3]
    # else
    # {
    tt <- 1:length(x)
    resid <- residuals(loess(x ~ tt))
    # }
    resid.q <- quantile(resid,prob=c(0.25,0.75))
    iqr <- diff(resid.q)
    limits <- resid.q + 1.5*iqr*c(-1,1)
    score <- abs(pmin((resid-limits[1])/iqr,0) + pmax((resid - limits[2])/iqr,0))
    if(plot)
    {
        plot(x)
        x2 <- ts(rep(NA,length(x)))
        x2[score>0] <- x[score>0]
        tsp(x2) <- tsp(x)
        points(x2,pch=19,col="red", cex=0.5)
        return(invisible(score))
    }
    else
        return(score)
}
# x$DO_mg.l[is.na(x$DO_mg.l)] = 10
# tsoutl(ts(x$DO_mg.l, deltat=1/288, start=188), TRUE) #for short
# deltat=1/96)[-(771:772)], TRUE) #for short
# tsoutl(ts(x$Abs.Pres..kPa..LGR.S.N..10838018..SEN.S.N..10838018.), TRUE)
tsoutl(tm, TRUE) #ugh


#winging it
diffs = diff(tm)

x11()
par(mfcol=c(2,1), mar=c(0,0,0,0), oma=c(5,5,0,0))
plot(tm, col='darkgreen', lwd=2, xaxt='n')
mtext('var', 2, 3, outer=FALSE, font=2)
# par(new=TRUE)
plot(diffs, col='steelblue3', lwd=2)
mtext('diff', 2, 3, outer=FALSE, font=2)
mtext('day', 1, 3, font=2)

u = mean(diffs, na.rm=TRUE)
sd = sd(diffs, na.rm=TRUE)

abline(h=c(-2*sd, 2*sd), lty=3, lwd=2, col='gray60')

pos_jumps = which(diffs > 2*sd)
neg_jumps = which(diffs < -2*sd)

# x = intersect(pos_jumps-1, neg_jumps)
# x = append(x, intersect(pos_jumps, neg_jumps-1))
# x = append(x, intersect(pos_jumps-2, neg_jumps))
# x = append(x, intersect(pos_jumps, neg_jumps-2))
# x = append(x, intersect(pos_jumps-3, neg_jumps))
# x = append(x, intersect(pos_jumps, neg_jumps-3))
# x = append(x, intersect(pos_jumps-4, neg_jumps))
# x = append(x, intersect(pos_jumps, neg_jumps-4))
# x = append(x, intersect(pos_jumps-5, neg_jumps))
# x = append(x, intersect(pos_jumps, neg_jumps-5))
#
# points(time(tm)[x], diffs[x], col='red', pch=20)
# plot(tm, xlim=c(26,28), type='b')


points(time(tm)[pos_jumps], diffs[pos_jumps], col='red', pch=20)
points(time(tm)[neg_jumps], diffs[neg_jumps], col='purple', pch=20)
# first_big_jump = min(c(pos_jumps, neg_jumps))
jump_inds = sort(c(pos_jumps, neg_jumps))
# jump_inds_s = sort(jump_inds)
library(accelerometry)
runs = rle2(as.numeric(jump_inds %in% pos_jumps), indices=TRUE)
lr = runs[,'lengths'] > 3
long_runs = runs[lr, 2:3]
# long_run_inds = jump_inds[unique(as.vector(long_runs))]

keep = numeric(length=nrow(long_runs))
for(i in 1:nrow(long_runs)){
    r = long_runs[i,]
    j = jump_inds[unique(c(r[1], r[1] + 1, r[2] - 1, r[2]))]
    t = time(tm)[j + 1] #+1 to convert from diff indices to original ts indices
    not_outlier = which.min(c(t[2] - t[1], t[length(t)] - t[length(t)-1]))
    keep[i] = r[-not_outlier]
}

posNeg_jump_pairs = sort(unique(c(keep, as.vector(runs[!lr,2:3]))))
outlier_inds = jump_inds[posNeg_jump_pairs]
points(time(tm)[outlier_inds], diffs[outlier_inds], col='orange', pch=20, cex=2)

# pairset = diff(sort(unique(r))) < 25
# pair_locator_1 = c(pairset[1], pairset)
# pair_locator_2 = c(pairset, pairset[length(pairset)])
# pairs = unique(r)[pair_locator_1 & pair_locator_2]
# pair_neighbs = c(posNeg_jump_pairs + 1, posNeg_jump_pairs - 1)
# pair_neighbs = pair_neighbs[pair_neighbs > 0 & pair_neighbs < length(jump_inds)]
# jump_inds[posNeg_jump_pairs] !%in%
# jump_inds[pair_neighbs]
# keep = numeric()
# for(i in posNeg_jump_pairs){
#     if(
# }

if(length(outlier_inds) %% 2 == 0){
    mapply(function(x, y), outlier_inds[seq(1, length(outlier_inds), 2)],
        outlier_inds[seq(2, length(outlier_inds), 2)])

plot(diffs, col='steelblue3', lwd=2, xlim=c(6.5,7))
points(time(tm)[outlier_inds+1], diffs[outlier_inds], col='orange', pch=20, cex=2)
plot(tm, col='darkgreen', lwd=2, xaxt='n', xlim=c(6.5,7))
points(time(tm)[outlier_inds], tm[outlier_inds], col='orange', pch=20, cex=2)

plot(diffs, col='steelblue3', lwd=2, xlim=c(17,18))
points(time(tm)[outlier_inds+1], diffs[outlier_inds], col='orange', pch=20, cex=2)
plot(tm, col='darkgreen', lwd=2, xaxt='n', xlim=c(17,18))
points(time(tm)[outlier_inds], tm[outlier_inds], col='orange', pch=20, cex=2)

plot(diffs, col='steelblue3', lwd=2, xlim=c(27,28))
points(time(tm)[outlier_inds+1], diffs[outlier_inds], col='orange', pch=20, cex=2)
plot(tm, col='darkgreen', lwd=2, xaxt='n', xlim=c(27,28))
points(time(tm)[outlier_inds], tm[outlier_inds], col='orange', pch=20, cex=2)
