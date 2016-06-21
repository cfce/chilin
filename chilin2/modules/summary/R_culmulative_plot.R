%# value

col=c('#FFB5C5','#5CACEE','#7CFC00','#FFD700','#8B475D','#8E388E','#FF6347','#FF83FA','#EEB422','#CD7054')
pch=c(21,22,24,25,21,22,24,25,21,22,24,25,21,22,24,25)

value = c(\VAR{current_data | join(", ")})
name = c(\VAR{ ids | surround_by_quote |join(',')})

pdf('\VAR{pdf}',height=10,width=10)
par(cex=1.8)

historicData<-c(\VAR{historic_data | join(", ")})
ecdf_func<-ecdf(historicData)

density <- ecdf_func(historicData)
value_quantile_pair <- data.frame(historicData,density)
value_quantile_pair_sorted<-value_quantile_pair[order(value_quantile_pair[,1]),]

plot(value_quantile_pair_sorted[,1],smooth(value_quantile_pair_sorted[,2]),type='l',pch=18,col='blue',main='\VAR{section}',xlab='\VAR{section}',ylab='fn(\VAR{section})')

xx_bad = value_quantile_pair_sorted[value_quantile_pair_sorted$historicData<=\VAR{cutoff},]
yy_bad = c(0,xx_bad[,2],0)
xx_bad =c(0,xx_bad[,1],\VAR{cutoff})
polygon(xx_bad,yy_bad, col = 'lightpink')

mid = quantile(historicData,0.5)
xx_warn = value_quantile_pair_sorted[value_quantile_pair_sorted$historicData<=mid,]
xx_warn = xx_warn[xx_warn[,1]>=\VAR{cutoff},]
yy_warn = c(0,xx_warn[,2],0)
xx_warn =c(\VAR{cutoff},xx_warn[,1],mid)
polygon(xx_warn,yy_warn, col = 'lightgoldenrod1')


xx_good = value_quantile_pair_sorted[value_quantile_pair_sorted$historicData>=mid,]
yy_good = c(0,xx_good[,2],0)
xx_good =c(mid,xx_good[,1],max(historicData))
polygon(xx_good,yy_good, col = 'palegreen')

len = length(value)


for (i in seq(len)){
points(value[i],jitter(ecdf_func(value[i])),pch=pch[i],bg=col[i])
}


legend('topleft',name,pch=pch[1:len],pt.bg=col[1:len])

dev.off()

