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
