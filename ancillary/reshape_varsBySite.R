d = read.csv('~/varsBySite_20180416.csv', stringsAsFactors=FALSE, header=FALSE)
head(d)
str(d)
colnames(d) = c('region', 'site', 'variable')

u = unique(d[,1:2])
u = u[order(u$region, u$site),]
v = sort(unique(d$variable))
# sites = paste(u$region, u$site, sep='_')
out = matrix(NA, nrow=nrow(u), ncol=length(v),
    dimnames=list(1:nrow(u), 1:length(v)))
for(i in 1:nrow(u)){
    reg = u$region[i]
    site = u$site[i]
    rownames(out)[i] = paste(reg, site, sep='_')
    for(j in 1:length(v)){
        if(v[j] %in% d$variable[which(d$region == reg & d$site == site)]){
            out[i,j] = v[j]
        } else {
            out[i,j] = ''
        }
    }
}
colnames(out) = v
write.csv(out, '~/variablesBySite_spread_20180416.csv')


lengths = vector()
for(i in 1:nrow(u)){
    reg = u$region[i]
    site = u$site[i]
    lengths = append(lengths,
        length(d$variable[which(d$region == reg & d$site == site)]))
}
maxl = max(lengths)

out2 = matrix(NA, nrow=nrow(u), ncol=maxl,
    dimnames=list(paste(u$region, u$site, sep='_'), 1:maxl))
for(i in 1:nrow(u)){
    reg = u$region[i]
    site = u$site[i]
    vars = sort(d$variable[which(d$region == reg & d$site == site)])
    vars = append(vars, rep('', maxl - length(vars)))
    out2[i,] = vars
}

write.csv(out2, '~/variablesBySite_compact_20180416.csv')
