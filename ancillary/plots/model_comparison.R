#this code comes from the compare_models function in sp_internals, with
#plotting code added

#define a function for determining whether to accept a new model if its
#coverage is lower than the existing best model (if penalty is
#proportionately lower, accept)
x1 = c(-365, 0, 365); y1 = c(1, 0, -1)
m1 = lm(y1 ~ x1)
a1 = m1$coefficients[2]
b1 = m1$coefficients[1]
above_line1 = pen_dif > a1 * coverage_dif + b1

#define a function for determining whether to accept a new model if its
#penalty is higher than the existing best model (if coverage is
#2X proportionately higher, accept)
x2 = c(-365, 0, 365); y2 = c(0.5, 0, -0.5)
m2 = lm(y2 ~ x2)
a2 = m2$coefficients[2]
b2 = m2$coefficients[1]
above_line2 = pen_dif > a2 * coverage_dif + b2

#plot ER-K600 correlation penalty and Kmax penalty ####
pdf(width=8, height=4, file='~/Dropbox/streampulse/figs/model_comparison1-2.pdf',
    compress=FALSE)

par(mfrow=c(1,2))

x=c(0.5, 1); y=0:1
m = lm(y ~ x)
plot(c(0,1), c(0,1), type='n', xaxs='i', yaxs='i', bty='l',
    xlab='ER x K600 Pearson Corr. (Abs)', ylab='ER-K penalty')
s = seq(0,90,0.1); p = predict(m, data.frame(x=s))
lines(c(s), c(unname(p)), lty=1, col='blue')

# dev.off()

#plot kmax penalty

# pdf(width=4, height=4, file='~/Dropbox/streampulse/figs/model_comparison2.pdf',
#     compress=FALSE)

x = c(45, 60, 80, 90); y = c(0, .5, .9, 1)
m2 = lm(y ~ x + I(x^2) + I(x^3))
plot(c(0,100), c(0,1), type='n', xaxs='i', yaxs='i', bty='l',
    xlab='Max daily K600 (Kmax)', ylab='Kmax penalty')
s = seq(0,90,0.1); p = predict(m2, data.frame(x=s))
lines(c(s,100), c(unname(p),1), lty=1, col='blue')

dev.off()


#plot coverage_dif x penalty_dif ####
pdf(width=4, height=4, file='~/Dropbox/streampulse/figs/model_comparison3.pdf',
    compress=FALSE)

plot(x1, y1, type='n', xaxs='i', yaxs='i', bty='l',
    xlab='Coverage Diff', ylab='Score Diff')
# lines(x1[1:2], fitted(m1)[1:2], lty=3)
# lines(x1[2:3], fitted(m2)[2:3], lty=3)
abline(h=0, lty=3, col='red')
abline(v=0, lty=3, col='red')
abline(h=-0.5, lty=3, col='red')
polygon(x=c(x1, rev(x1)),
    y=c(fitted(m1)[1:2], fitted(m2)[3], rep(2,3)),
    col='gray', border='gray')
points(c(-365,0,365), c(1,0,-0.5), pch=20, cex=.8, col='blue')

dev.off()
