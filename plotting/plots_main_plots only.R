# 
# 
cols = colorRampPalette(colors = c("black", "red", "yellow", "white"))(100)
colbr_h = seq(0,2.5,length.out=101)
colbr_r = seq(0,100,length.out=101)

hor_cut = length(kivec)
layout(rbind(c(1,2),c(1,2),c(1,2),c(1,2),c(3,4)))
par(mar=c(4,6,1,1), oma=c(5,1,1,1), cex.axis=1.8)
image(t(hdiffmat[,1:hor_cut]), xlab="", ylab="", xaxt="no", yaxt="no", col = cols, breaks = colbr_h)
axis(side=1, at=seq(0,length(kivec[1:hor_cut]), length.out=5)/hor_cut, labels=sprintf("%.4g", kivec[1:hor_cut])[as.integer(seq(1,hor_cut, length.out=5))])
axis(side=2, at=seq(0,length(rivec), length.out=5)/25, labels=sprintf("%.2g", rivec)[as.integer(seq(1,length(chvec), length.out=5))])
mtext("Imitation kernel size", side = 1, line = 3, outer = F, cex=1.2)
mtext("Imitation rate", side = 2, line = 3, outer = F, cex=1.2)
# mtext("Spatial correlation length of \nresource growth rate", side = 2, line = 3, outer = F, cex=1.2)

image(t(percapmat[,1:hor_cut]), xlab="", ylab="", xaxt="no", yaxt="no", col = cols, breaks = colbr_r)
axis(side=1, at=seq(0,length(kivec[1:hor_cut]), length.out=5)/hor_cut, labels=sprintf("%.4g", kivec[1:hor_cut])[as.integer(seq(1,hor_cut, length.out=5))])
axis(side=2, at=seq(0,length(rivec), length.out=5)/25, labels=sprintf("%.2g", rivec)[as.integer(seq(1,length(chvec), length.out=5))])
mtext("Imitation kernel size", side = 1, line = 3, outer = F, cex=1.2)
mtext("Imitation rate", side = 2, line = 3, outer = F, cex=1.2)
# mtext("Spatial correlation length of \nresource growth rate", side = 2, line = 3, outer = F, cex=1.2)

image(as.matrix(colbr_h), y=1, x=as.matrix(colbr_h), col=cols, yaxt="n", xlab="", ylab="")
mtext("Difference between harvesting\n rates of milkers and killers", side = 1, line = 5, outer = F, cex=1.2)
image(as.matrix(colbr_r), y=1, x=as.matrix(colbr_r), col=cols, yaxt="n", xlab="", ylab="")
mtext("Average per capita \nresource consumed", side = 1, line = 5, outer = F, cex=1.2)


# Plot timeseries

# id = seq(0,7500, by=10)
# id[1]=1
# layout(rbind(c(1,1,1,1,4,4), c(2,2,2,2,5,5), c(3,3,3,3,6,6)))
# par(cex.axis=1.8, cex.lab=1.8, mar=c(4,7,1,1), oma=c(1,1,1,1))
# cols = colorRampPalette(colors = c("black", "blue", "cyan", "white"))(100)
# image(x = id*100/1e5, y= c(brks_h[1:nbins],brks_h[nbins]+0.01), log(hall_h[id,]+1), col = cols, xlab="Time", ylab="Harvesting \n rate")
# image(x = id*100/1e5, y= brks_kd[1:(nbins/2.5)], log(hall_k[id,1:(nbins/2.5)]+1), col = cols, xlab="Time", ylab="Dispersal \n distance")
# image(x = id*100/1e5, y= c(brks_rc[1:nbins],brks_rc[nbins]+0.01), log(hall_r[id,]+1), col = cols, xlab="Time", ylab="Resource \n consumed")
# # points(x=times, y=r_tot/1.0125e7*600, type="l", col="white")
# # points(x=id*100, y=filter(avgrc*6, filter=rep(1,20)/20)[id], type="l", col="green")
# 
# # par(cex.axis=1.8, cex.lab=1.8, mar=c(4,0,1,1), oma=c(1,1,1,1))
# mids_h=c((brks_h[1:99]+brks_h[1:99+1])/2,2.5)
# mids_kd=c((brks_kd[1:99]+brks_kd[1:99+1])/2,25)
# mids_rc=c((brks_rc[1:99]+brks_rc[1:99+1])/2,600)
# plot(mids_h~h_scan[ib,], type="l", ylab="", yaxt="n", xlab="Frequency", col="blue", lwd=2)
# plot(mids_kd[1:(nbins/2.5)]~kd_scan[ib,1:(nbins/2.5)], type="l", ylab="",  yaxt="n", xlab="Frequency", col="blue", lwd=2)
# plot(mids_rc~rc_scan[ib,], type="l", ylab="",  yaxt="n", xlab="Frequency", col="blue", lwd=2)
