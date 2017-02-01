setwd("/home/jaideep/austria_project/figures_final_data/map_bxch")

h_lo_mat = read.csv(header = T, file="h_lo.csv", row.names = 1)
h_hi_mat = read.csv(header = T, file="h_hi.csv", row.names = 1)
h_avg_mat = read.csv(header = T, file="h_avg.csv", row.names = 1)
h_diff_mat = h_hi_mat - h_lo_mat

kd_lo_mat = read.csv(header = T, file="kd_lo.csv", row.names = 1)
kd_hi_mat = read.csv(header = T, file="kd_hi.csv", row.names = 1)
kd_avg_mat = read.csv(header = T, file="kd_avg.csv", row.names = 1)
kd_diff_mat = kd_hi_mat - kd_lo_mat

ld_lo_mat = read.csv(header = T, file="ld_lo.csv", row.names = 1)
ld_hi_mat = read.csv(header = T, file="ld_hi.csv", row.names = 1)
ld_avg_mat = read.csv(header = T, file="ld_avg.csv", row.names = 1)
ld_diff_mat = ld_hi_mat - ld_lo_mat

rc_avg_mat = read.csv(header = T, file="rc_avg.csv", row.names = 1)
r_env_mat = read.csv(header = T, file="r_env.csv", row.names = 1)

image.natural = function(x, xlabels, ylabels, xlab, ylab, cols, col_brks, 
                         xfmt = "%.1g", yfmt = "%.1g", hor_cut=-1, 
                         title=""){
  
  layout(rbind(c(1),c(1),c(1),c(1),c(2)))
  par(mar=c(4,6,1,1), oma=c(5,1,1,1), cex.axis=1.8)
  
  if(hor_cut==-1) hor_cut = dim(x)[2]
  image(t(x[,1:hor_cut]), xlab="", ylab="", xaxt="no", yaxt="no", col = cols, breaks = col_brks)
  axis(side=1, at=seq(0,length(xlabels[1:hor_cut]), length.out=5)/hor_cut, labels=sprintf(xfmt, xlabels[1:hor_cut])[as.integer(seq(1,hor_cut, length.out=5))])
  axis(side=2, at=seq(0,length(ylabels), length.out=5)/25, labels=sprintf(yfmt, ylabels)[as.integer(seq(1,length(ylabels), length.out=5))])
  mtext(xlab, side = 1, line = 3, outer = F, cex=1.2)
  mtext(ylab, side = 2, line = 3, outer = F, cex=1.2)
  
  image(as.matrix(col_brks), y=1, x=as.matrix(col_brks), col=cols, yaxt="n", xlab="", ylab="")
  mtext(title, side = 1, line = 5, outer = F, cex=1.2)
  
}

# 
# 
cols = colorRampPalette(colors = c("black", "red", "yellow", "white"))(100)
colbr_h = seq(0,2.5,length.out=101)
colbr_rc = seq(0,86.6,length.out=101)
colbr_ld = seq(0,6,length.out=101)
colbr_r = seq(0,50,length.out=101)


bvec = exp(seq(log(0.0002), log(0.2), length.out=50))[(1:50)%%2 != 0]/cd
kivec = exp(seq(log(10), log(1000), length.out=25))/L#[-3]
rivec = exp(seq(log(0.001), log(1), length.out=25))/r
muvec = exp(seq(log(0.01), log(1), length.out=25))
chvec = exp(seq(log(0.001), log(10), length.out=25))/cd

npar = 25


hor_cut = 25

setEPS(width=3.76, height=4.69)
postscript("ld_diff.eps")
image.natural(x = ld_diff_mat, xlabels = bvec, ylabels = chvec, 
              cols = cols, 
              xlab = "Benefit of harvesting", 
              ylab = "Cost of harvesting", 
              col_brks = colbr_ld,
              yfmt = "%.3g",
              title = "Difference between dispersal\ndistance of milkers and killers")
dev.off()




image(t(rc_avg_mat[,1:hor_cut]), xlab="", ylab="", xaxt="no", yaxt="no", col = cols, breaks = colbr_r)
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

id = seq(0,7500, by=10)
id[1]=1
layout(rbind(c(1,1,1,1,4,4), c(2,2,2,2,5,5), c(3,3,3,3,6,6)))
par(cex.axis=1.8, cex.lab=1.8, mar=c(4,7,1,1), oma=c(1,1,1,1))
cols = colorRampPalette(colors = c("black", "blue", "cyan", "white"))(100)
image(x = id*100/1e5, y= c(brks_h[1:nbins],brks_h[nbins]+0.01), log(hall_h[id,]+1), col = cols, xlab="Time", ylab="Harvesting \n rate")
image(x = id*100/1e5, y= brks_kd[1:(nbins/2.5)], log(hall_k[id,1:(nbins/2.5)]+1), col = cols, xlab="Time", ylab="Dispersal \n distance")
image(x = id*100/1e5, y= c(brks_rc[1:nbins],brks_rc[nbins]+0.01), log(hall_r[id,]+1), col = cols, xlab="Time", ylab="Resource \n consumed")
# points(x=times, y=r_tot/1.0125e7*600, type="l", col="white")
# points(x=id*100, y=filter(avgrc*6, filter=rep(1,20)/20)[id], type="l", col="green")

# par(cex.axis=1.8, cex.lab=1.8, mar=c(4,0,1,1), oma=c(1,1,1,1))
mids_h=c((brks_h[1:99]+brks_h[1:99+1])/2,2.5)
mids_kd=c((brks_kd[1:99]+brks_kd[1:99+1])/2,25)
mids_rc=c((brks_rc[1:99]+brks_rc[1:99+1])/2,600)
plot(mids_h~h_scan[ib,], type="l", ylab="", yaxt="n", xlab="Frequency", col="blue", lwd=2)
plot(mids_kd[1:(nbins/2.5)]~kd_scan[ib,1:(nbins/2.5)], type="l", ylab="",  yaxt="n", xlab="Frequency", col="blue", lwd=2)
plot(mids_rc~rc_scan[ib,], type="l", ylab="",  yaxt="n", xlab="Frequency", col="blue", lwd=2)
