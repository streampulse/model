#take one (ful series uploads) ####

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


#winging it (four weeks) ####

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

library(imputeTS)
tm = ts(d$value, deltat = 1/96)
tm = na.seadec(tm, algorithm='interpolation')


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
    outlier_ts = mapply(function(x, y){ seq(x+1, y, 1) },
        outlier_inds[seq(1, length(outlier_inds), 2)],
        outlier_inds[seq(2, length(outlier_inds), 2)],
        SIMPLIFY=FALSE)
    outlier_ts = unlist(outlier_ts)
}

# plot(diffs, col='steelblue3', lwd=2)#, xlim=c(6.5,7))
# points(time(tm)[outlier_ts+1], diffs[outlier_ts], col='orange', pch=20, cex=2)
plot(tm, col='darkgreen', lwd=2, xaxt='n')#, xlim=c(6.5,7))
points(time(tm)[outlier_ts], tm[outlier_ts], col='orange', pch=20, cex=2)

# plot(diffs, col='steelblue3', lwd=2, xlim=c(6.5,7))
# points(time(tm)[outlier_inds+1], diffs[outlier_inds], col='orange', pch=20, cex=2)
plot(tm, col='darkgreen', lwd=2, xaxt='n', xlim=c(6.5,7))
points(time(tm)[outlier_ts], tm[outlier_ts], col='orange', pch=20, cex=2)

# plot(diffs, col='steelblue3', lwd=2, xlim=c(17,18))
# points(time(tm)[outlier_inds+1], diffs[outlier_inds], col='orange', pch=20, cex=2)
plot(tm, col='darkgreen', lwd=2, xaxt='n', xlim=c(17,18))
points(time(tm)[outlier_ts], tm[outlier_ts], col='orange', pch=20, cex=2)

# plot(diffs, col='steelblue3', lwd=2, xlim=c(27,28))
# points(time(tm)[outlier_inds+1], diffs[outlier_inds], col='orange', pch=20, cex=2)
plot(tm, col='darkgreen', lwd=2, xaxt='n', xlim=c(27,28))
points(time(tm)[outlier_ts], tm[outlier_ts], col='orange', pch=20, cex=2)

# testing left- and right-border outlier cutoffs ####
#also pos/neg initial jump (relative to first non-outlier jump, etc.)

library(imputeTS)
setwd('/home/mike/Dropbox/streampulse/data/NC_download/')
site_files = dir()
d_chili = read.csv(site_files[1], stringsAsFactors=FALSE,
    colClasses=c('DateTime_UTC'='POSIXct'))
sitevars = unique(d_chili$variable)

#mess with input series here
v=7
# d = d[d$DateTime_UTC >= as.POSIXct('2016-10-05 17:45:00') & d$DateTime_UTC <
#         as.POSIXct('2016-10-11 17:45:00'),]
d = d_chili[d_chili$variable == sitevars[v],]

#get xvals from previous plot
library(plotrix)
loc = locator()
locg = rescale(c(1,loc$x[1],loc$x[2],max(time(tm))), c(1,nrow(d)))
loc2 = locator()
locg = rescale(c(1,loc2$x[1],loc2$x[2],max(time(tm))), c(1,nrow(d)))
d = d[round(locg[2]):round(locg[3]),]

rm(list=c('sd_scaler','big_jump_prop','pos_jumps','neg_jumps','outlier_inds','outlier_ts',
    'posNeg_jump_pairs','lr','diffs'))


#change values
# derp = which(d$DateTime_UTC %in% chili)
# d[derp,]
# d$value[derp] = 110

tm = ts(d$value, deltat = 1/96)

ts_u = mean(tm, na.rm=TRUE)
ts_sd = sd(tm, na.rm=TRUE)
# mad = sum(abs(tm - ts_u)) / sum(!is.na(tm))
big_outliers_h = which(tm > ts_u + (4*ts_sd))
big_outliers_l = which(tm < ts_u - (4*ts_sd))
big_outliers = unique(c(big_outliers_h, big_outliers_l))

tm = na.seadec(tm, algorithm='interpolation')

#real stuff
diffs = diff(tm)
# x11()
par(mfcol=c(3,1), mar=c(0,0,0,0), oma=c(5,5,0,0))
plot(tm, col='darkgreen', lwd=2, xaxt='n')
# abline(h=ts_u, lty=2, col='black')
abline(h=c(ts_u - (4*ts_sd), ts_u + (4*ts_sd)), lty=2, col='black')
points(time(tm)[big_outliers], tm[big_outliers], col='red')
mtext('var', 2, 3, outer=FALSE, font=2)
# par(new=TRUE)
plot(diffs, col='steelblue3', lwd=2)
mtext('diff', 2, 3, outer=FALSE, font=2)
mtext('day', 1, 3, font=2)

u = mean(diffs, na.rm=TRUE)
sd = sd(diffs, na.rm=TRUE)

sd_scaler = 1.8
big_jump_prop = Inf
while(big_jump_prop > 0.03){
    sd_scaler = sd_scaler + 0.2
    abline(h=c(-sd_scaler*sd, sd_scaler*sd), lty=3, lwd=2, col='gray60')
    big_jump_prop = sum(diffs > sd_scaler * sd) / length(diffs)
}

pos_jumps = which(diffs > sd_scaler * sd)
neg_jumps = which(diffs < -sd_scaler * sd)

points(time(tm)[pos_jumps], diffs[pos_jumps], col='red', pch=20)
points(time(tm)[neg_jumps], diffs[neg_jumps], col='purple', pch=20)
# first_big_jump = min(c(pos_jumps, neg_jumps))
jump_inds = sort(c(pos_jumps, neg_jumps))
# jump_inds_s = sort(jump_inds)
library(accelerometry)

runs = rle2(as.numeric(jump_inds %in% pos_jumps), indices=TRUE)
# run_start_inds = jump_inds[runs[,'starts']]
lr = runs[,'lengths'] > 3
long_runs = runs[lr, 2:3, drop=FALSE]
# long_run_inds = jump_inds[unique(as.vector(long_runs))]

keep = numeric()
if(length(long_runs)){
    for(i in 1:nrow(long_runs)){
        r = long_runs[i,]
        j = jump_inds[unique(c(r[1], r[1] + 1, r[2] - 1, r[2]))]
        if(abs(j[2]-j[1]) < 8 & abs(j[length(j)]-j[length(j)-1]) < 9) next
        t = time(tm)[j + 1] #+1 to convert from diff indices to original ts indices
        left_jump_interv = t[2] - t[1]
        right_jump_interv = t[length(t)] - t[length(t)-1]
        not_outlier = which.min(c(left_jump_interv, right_jump_interv))
        keep = append(keep, r[-not_outlier])
    }
}

posNeg_jump_pairs = sort(unique(c(keep, as.vector(runs[!lr,2:3]))))
outlier_inds = jump_inds[posNeg_jump_pairs]
points(time(tm)[outlier_inds], diffs[outlier_inds], col='orange', pch=20, cex=2)

# inds = outlier_inds
# inds = multijumps
# inds=inds[1:2,]
seq_map = function(inds, x_advance=0, y_advance=0){
    seq_list = mapply(function(x, y){ seq(x+x_advance, y+y_advance, 1) },
        # inds[seq(1, length(inds), 2)],
        # inds[seq(2, length(inds), 2)],
        inds[,1], inds[,2],
        SIMPLIFY=FALSE)
    seq_vec = unlist(seq_list)

    return(seq_vec)
}
# if(length(outlier_inds) == 1){
#     if(length(tm) - outlier_inds <= 15){
#         outlier_ts = seq(outlier_inds+1, length(tm), 1)
#     } else {
#         if(outlier_inds - 1 <= 15){
#             outlier_ts = seq(1, outlier_inds, 1)
#         } else {
#             outlier_ts = NULL
#         }
#     }

n_outlier_pieces = Inf
counter = 0
if(length(outlier_inds) == 1){
    outlier_ts = NULL
} else {
    while(length(outlier_inds) > 1 & n_outlier_pieces > 50 & counter < 5){
        outdif = diff(outlier_inds)
        rm_multjump = rm_oneway = NULL

        short_jumps = outdif < 15
        if(!short_jumps[1]){
            outlier_inds = outlier_inds[-1]
            counter = counter + 1
            next
        }

        jump_runs = rle2(as.numeric(short_jumps), indices=TRUE,
            return.list=FALSE)
        multijumps = jump_runs[jump_runs[,'lengths'] > 1, 2:3, drop=FALSE]
        if(length(multijumps)){
            multijumps = multijumps[short_jumps[multijumps[,'starts']], ,
                drop=FALSE]
            rm_multjump = seq_map(multijumps, y_advance=1)
            # outlier_inds = outlier_inds[-rm_multjump]
        }

        big_outdif = outdif > 15
        if(big_outdif[length(big_outdif)]){
            outlier_inds = outlier_inds[-length(outlier_inds)]
            counter = counter + 1
            next
        }

        if(all(big_outdif)){
            outlier_ts = NULL
            break
        } else {
        # if(big_outdif[1]){
            # rm_oneway = outlier_inds[1]
        # } else {
            same_sign_runs = rle2(as.numeric(big_outdif), indices=TRUE,
                return.list=TRUE)

            if(length(same_sign_runs)){
                l = same_sign_runs$lengths
                one_way_jumps = numeric()
                for(i in 1:length(l)){
                    if(l[i] > 1){
                        s = same_sign_runs$starts[i]
                        one_way_jumps = append(one_way_jumps,
                            seq(s, s + (l[i] - 2), 1))
                            # same_sign_runs$starts[i]:(i + (l[i] - 2))])
                    }
                }
            }

            # if(same_sign_runs$lengths > 1)
            # one_way_jumps = same_sign_runs$starts[]
            if(length(one_way_jumps)){
                one_way_jumps = one_way_jumps[big_outdif[one_way_jumps]]
                outdif_filt = outdif
                outdif_filt[-one_way_jumps] = 0
                rm_oneway = which(outdif_filt > 0) + 1
                # rm_oneway = which.max(outdif_filt) + 1

                # if(!all(big_outdiff)){
                    # outlier_inds = outlier_inds[-1]
                # } else {
                    # outlier_ts = NULL
                # }
            }
        }
        removals = unique(c(rm_multjump, rm_oneway))
        if(!is.null(removals)){
            outlier_inds = outlier_inds[-removals]
        }

        if(length(outlier_inds) %% 2 == 0){
            outlier_inds = matrix(outlier_inds, ncol=2, byrow=TRUE)
            outlier_ts = seq_map(outlier_inds, x_advance=1)
        } else {
            smallest_diff = which.min(abs(diffs[outlier_inds]))
            outlier_inds = outlier_inds[-smallest_diff]
            counter = counter + 1
            next
        }



        # } else {
            # min_outl = which.min(abs(diffs[outlier_inds]))
            # outlier_inds = outlier_inds[-min_outl]
        # }
        # outlier_inds[which(rle2(outdif)$lengths == 2) + 1]
        # rle(sign(diff(outlier_inds, differences=2)))$lengths
        # rm_dd = ifelse(dd[1] - dd[length(dd)] > 0, 1, length(outlier_inds))
        # outlier_inds = outlier_inds[-rm_dd]
        # outlier_ts = seq_map(outlier_inds, x_advance=1)
        # }
        n_outlier_pieces = length(outlier_ts)
        counter = counter + 1
    }
}
if(n_outlier_pieces > 500){
    outlier_ts = NULL
    print('argh2')
}

outlier_ts = unique(c(outlier_ts, big_outliers))
#         if(rm_dd == 1 & outlier_inds[1] <= 15){
#             left_bounded_outlier_seq = seq(1, outlier_inds[1], 1)
#             outlier_inds = outlier_inds[-1]
#             outlier_ts = seq_map(outlier_inds)
#             outlier_ts = c(outlier_ts, left_bounded_outlier_seq)
#         } else {
#             right_bounded_outlier_seq = seq(outlier_inds[length(outlier_inds)]+1,
#                 length(tm), 1)
#             if(rm_dd == 2 & length(right_bounded_outlier_seq) <= 15) {
#                 outlier_inds = outlier_inds[-length(outlier_inds)]
#                 outlier_ts = seq_map(outlier_inds)
#                 outlier_ts = c(outlier_ts, right_bounded_outlier_seq)
#             }
#         }
#     }
# }
# if(length(outlier_ts) > 50) outlier_ts = NULL
# outlier_ts = outlier_inds
# plot(diffs, col='steelblue3', lwd=2)#, xlim=c(6.5,7))
# points(time(tm)[outlier_ts+1], diffs[outlier_ts], col='orange', pch=20, cex=2)
plot(tm, col='darkgreen', lwd=2, xaxt='n')#, xlim=c(6.5,7))
points(time(tm)[outlier_ts], tm[outlier_ts], col='orange', pch=20, cex=2)
(chili = d$DateTime_UTC[outlier_ts])
print(paste('sd_scaler =', sd_scaler))
print(paste('big_jump_prop =', big_jump_prop))

# plot(tm, col='darkgreen', lwd=2, xaxt='n', xlim=c(16.5,16.8))
# points(time(tm)[outlier_ts], tm[outlier_ts], col='orange', pch=20, cex=2)
# axis(1, time(tm)[outlier_ts], tm[outlier_ts], las=2)


# plot(tm, col='darkgreen', lwd=2, xaxt='n', xlim=c(21.5,22.5))
# # plot(tm, col='darkgreen', lwd=2, xaxt='n', xlim=c(.8,1.5))
# points(time(tm)[outlier_ts], tm[outlier_ts], col='orange', pch=20, cex=2)
# axis(1, time(tm), d$DateTime_UTC, las=2)

#testing the function version used in app.py ####

# should copy/paste the newest version before using this

#! /usr/bin/Rscript

#TODO: make this work with any sample interval

# df = commandArgs(trailingOnly=TRUE)
# print(head(df))


setwd('/home/mike/Dropbox/streampulse/data/NC_download/')
# site_files = dir()
# d_chili = read.csv(site_files[1], stringsAsFactors=FALSE,
#     colClasses=c('DateTime_UTC'='POSIXct'))
# sitevars = unique(d_chili$variable)
#
# v=1
# d = d_chili[d_chili$variable == sitevars[v],]

# df = read.csv('../test_outl.csv', stringsAsFactors=FALSE)
df = read.csv('../test_outl2.csv', stringsAsFactors=FALSE)



# seq_map = function(inds, x_advance=0, y_advance=0){
#     seq_list = mapply(function(x, y){ seq(x+x_advance, y+y_advance, 1) },
#         inds[,1], inds[,2],
#         SIMPLIFY=FALSE)
#     seq_vec = unlist(seq_list)
#
#     return(seq_vec)
# }

find_outliers = function(df){

    library(imputeTS)
    library(plotrix)
    library(accelerometry)

    df = subset(df, select=-c(DateTime_UTC)) #remove datetime col

    outlier_list = list()

    for(col in 1:ncol(df)){

        # print(colnames(df)[col])

        if(sum(is.na(df[,col])) / nrow(df) > 0.98){ #if almost all NA
            outlier_list[[col]] = 'NONE'
            print(paste(1, names(outlier_list), colnames(df)[col]))
            names(outlier_list)[col] = colnames(df)[col]
            next
        }

        tm = ts(df[,col], deltat = 1/96)
        tm = na.seadec(tm, algorithm='interpolation')

        #real stuff
        diffs = diff(tm)
        # x11()
        u = mean(diffs, na.rm=TRUE)
        sd = sd(diffs, na.rm=TRUE)

        sd_scaler = 1.8
        big_jump_prop = Inf
        while(big_jump_prop > 0.03){
            sd_scaler = sd_scaler + 0.2
            big_jump_prop = sum(diffs > sd_scaler * sd) / length(diffs)
        }

        pos_jumps = which(diffs > sd_scaler * sd)
        neg_jumps = which(diffs < -sd_scaler * sd)

        jump_inds = sort(c(pos_jumps, neg_jumps))

        runs = rle2(as.numeric(jump_inds %in% pos_jumps), indices=TRUE)
        lr = runs[,'lengths'] > 3
        long_runs = runs[lr, 2:3, drop=FALSE]

        keep = numeric()
        if(length(long_runs)){
            for(i in 1:nrow(long_runs)){
                r = long_runs[i,]
                j = jump_inds[unique(c(r[1], r[1] + 1, r[2] - 1, r[2]))]
                if(abs(j[2]-j[1]) < 8 & abs(j[length(j)]-j[length(j)-1]) < 9) next
                t = time(tm)[j + 1] #+1 to convert from diff indices to original ts indices
                left_jump_interv = t[2] - t[1]
                right_jump_interv = t[length(t)] - t[length(t)-1]
                not_outlier = which.min(c(left_jump_interv, right_jump_interv))
                keep = append(keep, r[-not_outlier])
            }
        }

        posNeg_jump_pairs = sort(unique(c(keep, as.vector(runs[!lr,2:3]))))
        outlier_inds = jump_inds[posNeg_jump_pairs]

        n_outlier_pieces = Inf
        counter = 0
        if(length(outlier_inds) == 1){
            outlier_ts = 'NONE'
        } else {
            while(length(outlier_inds) > 1 & n_outlier_pieces > 50 & counter < 5){
                outdif = diff(outlier_inds)
                rm_multjump = rm_oneway = NULL

                short_jumps = outdif < 15
                if(!short_jumps[1]){
                    outlier_inds = outlier_inds[-1]
                    counter = counter + 1
                    next
                }

                jump_runs = rle2(as.numeric(short_jumps), indices=TRUE,
                    return.list=FALSE)
                multijumps = jump_runs[jump_runs[,'lengths'] > 1, 2:3, drop=FALSE]
                if(length(multijumps)){
                    multijumps = multijumps[short_jumps[multijumps[,'starts']], ,
                        drop=FALSE]
                    seq_list = mapply(function(x, y){ seq(x, y+1, 1) },
                        multijumps[,1], multijumps[,2],
                        SIMPLIFY=FALSE)
                    rm_multjump = unlist(seq_list)
                    # rm_multjump = seq_map(multijumps, y_advance=1)
                }

                big_outdif = outdif > 15
                if(big_outdif[length(big_outdif)]){
                    outlier_inds = outlier_inds[-length(outlier_inds)]
                    counter = counter + 1
                    next
                }

                if(all(big_outdif)){
                    outlier_ts = 'NONE'
                    break
                } else {
                    same_sign_runs = rle2(as.numeric(big_outdif), indices=TRUE,
                        return.list=TRUE)

                    if(length(same_sign_runs)){
                        l = same_sign_runs$lengths
                        one_way_jumps = numeric()
                        for(i in 1:length(l)){
                            if(l[i] > 1){
                                s = same_sign_runs$starts[i]
                                one_way_jumps = append(one_way_jumps,
                                    seq(s, s + (l[i] - 2), 1))
                            }
                        }
                    }

                    if(length(one_way_jumps)){
                        one_way_jumps = one_way_jumps[big_outdif[one_way_jumps]]
                        outdif_filt = outdif
                        outdif_filt[-one_way_jumps] = 0
                        rm_oneway = which(outdif_filt > 0) + 1
                    }
                }
                removals = unique(c(rm_multjump, rm_oneway))
                if(!is.null(removals)){
                    outlier_inds = outlier_inds[-removals]
                }

                if(length(outlier_inds) %% 2 == 0){
                    outlier_inds = matrix(outlier_inds, ncol=2, byrow=TRUE)
                    seq_list = mapply(function(x, y){ seq(x+1, y, 1) },
                        outlier_inds[,1], outlier_inds[,2],
                        SIMPLIFY=FALSE)
                    outlier_ts = unlist(seq_list)
                    # outlier_ts = seq_map(outlier_inds, x_advance=1)
                } else {
                    smallest_diff = which.min(abs(diffs[outlier_inds]))
                    outlier_inds = outlier_inds[-smallest_diff]
                    counter = counter + 1
                    next
                }
                n_outlier_pieces = length(outlier_ts)
                counter = counter + 1
            }
        }
        if(n_outlier_pieces > 50 | is.null(outlier_ts)){
            outlier_ts = 'NONE'
            print('argh2')
        }
        # print(paste('sd_scaler =', sd_scaler))
        # print(paste('big_jump_prop =', big_jump_prop))

        outlier_list[[col]] = outlier_ts
        print(paste(2, names(outlier_list), colnames(df)[col]))
        names(outlier_list)[col] = colnames(df)[col]
    }

    return(outlier_list)
}

find_outliers(df)
