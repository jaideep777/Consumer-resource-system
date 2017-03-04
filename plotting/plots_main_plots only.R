setwd("/home/jaideep/austria_project/figures_final_data/map_muxkI_b0.00144_rI0.1")

r = 0.2
K = 50
MSY = 89
sigE = 4
cD = 0.1


h_lo_mat = read.csv(header = T, file="h_lo.csv", row.names = 1)/r
h_hi_mat = read.csv(header = T, file="h_hi.csv", row.names = 1)/r
h_avg_mat = read.csv(header = T, file="h_avg.csv", row.names = 1)/r
h_diff_mat = h_hi_mat - h_lo_mat

kd_lo_mat = read.csv(header = T, file="kd_lo.csv", row.names = 1)/sigE
kd_hi_mat = read.csv(header = T, file="kd_hi.csv", row.names = 1)/sigE
kd_avg_mat = read.csv(header = T, file="kd_avg.csv", row.names = 1)/sigE
kd_diff_mat = kd_hi_mat - kd_lo_mat

ld_lo_mat = read.csv(header = T, file="ld_lo.csv", row.names = 1)/r/sigE
ld_hi_mat = read.csv(header = T, file="ld_hi.csv", row.names = 1)/r/sigE
ld_avg_mat = read.csv(header = T, file="ld_avg.csv", row.names = 1)/r/sigE
ld_diff_mat = ld_hi_mat - ld_lo_mat

rc_avg_mat = read.csv(header = T, file="rc_avg.csv", row.names = 1)/MSY
r_env_mat = read.csv(header = T, file="r_env.csv", row.names = 1)/450/450/K

image.natural = function(x, xlabels, ylabels, xlab, ylab, cols, col_brks, 
                         xfmt = "%.3g", yfmt = "%.1g", hor_cut=-1, ver_cut=-1,
                         title=""){
  
  layout(rbind(c(1),c(1),c(1),c(1),c(2)))
  par(mar=c(4,6,1,1), oma=c(5,1,1,1), cex.axis=1.8)
  
  if(hor_cut==-1) hor_cut = dim(x)[2]
  if(ver_cut==-1) ver_cut = 1
  image(t(x[ver_cut:25,1:hor_cut]), xlab="", ylab="", xaxt="no", yaxt="no", col = cols, breaks = col_brks)
  axis(side=1, at=seq(0,length(xlabels[1:hor_cut]), length.out=5)/hor_cut, labels=sprintf(xfmt, xlabels[1:hor_cut])[as.integer(seq(1,hor_cut, length.out=5))])
  axis(side=2, at=seq(0,length(ylabels[ver_cut:25]), length.out=5)/(26-ver_cut), labels=sprintf(yfmt, ylabels[ver_cut:25])[as.integer(seq(1,26-ver_cut, length.out=5))])
  mtext(xlab, side = 1, line = 3, outer = F, cex=1.2)
  mtext(ylab, side = 2, line = 3, outer = F, cex=1.2)
  
  image(as.matrix(col_brks), y=1, x=as.matrix(col_brks), col=cols, yaxt="n", xlab="", ylab="")
  mtext(title, side = 1, line = 5, outer = F, cex=1.2)
  
}

image.natural.nonlinear = function(x, col_data, xlabels, ylabels, xlab, ylab, cols, col_brks, 
                         xfmt = "%.3g", yfmt = "%.1g", hor_cut=-1, ver_cut=-1,
                         title=""){
  
  layout(rbind(c(1),c(1),c(1),c(1),c(2)))
  par(mar=c(4,6,1,1), oma=c(5,1,1,1), cex.axis=1.8)
  
  if(hor_cut==-1) hor_cut = dim(x)[2]
  if(ver_cut==-1) ver_cut = 1
  image(t(x[ver_cut:25,1:hor_cut]), xlab="", ylab="", xaxt="no", yaxt="no", col = cols, breaks = col_brks)
  axis(side=1, at=seq(0,length(xlabels[1:hor_cut]), length.out=5)/hor_cut, labels=sprintf(xfmt, xlabels[1:hor_cut])[as.integer(seq(1,hor_cut, length.out=5))])
  axis(side=2, at=seq(0,length(ylabels[ver_cut:25]), length.out=5)/(26-ver_cut), labels=sprintf(yfmt, ylabels[ver_cut:25])[as.integer(seq(1,26-ver_cut, length.out=5))])
  mtext(xlab, side = 1, line = 3, outer = F, cex=1.2)
  mtext(ylab, side = 2, line = 3, outer = F, cex=1.2)
  
  image(as.matrix(col_brks), y=1, x=as.matrix(col_data), col=cols, yaxt="n", xlab="", ylab="")
  mtext(title, side = 1, line = 5, outer = F, cex=1.2)
  
}


# 
# 
cols = colorRampPalette(colors = c("black", "red", "yellow", "white"))(100)
colbr_hdiff = seq(0,2.5/r,length.out=101)
colbr_havg = seq(0,.8/r,length.out=101)
colbr_rc = seq(.5,1,length.out=101)
colbr_ldavg = seq(0,2/r/sigE,length.out=101)
colbr_lddiff = seq(0,16/r/sigE,length.out=101)
colbr_r = seq(0,50/K,length.out=101)
colbr_kddiff = seq(0,25/sigE,length.out=101)
colbr_kdavg = seq(0,25/sigE,length.out=101)

cols = colorRampPalette(colors = c("black", "red", "yellow", "white"))(100)
colbr_hdiff = seq(0,2.5/r,length.out=101)
colbr_havg = seq(0,1.5/r,length.out=101)
colbr_rc = seq(0,1,length.out=101)
colbr_ldavg = seq(0,4.8/r/sigE,length.out=101)
colbr_lddiff = seq(0,6.4/r/sigE,length.out=101)
colbr_r = seq(0,50/K,length.out=101)
colbr_kddiff = seq(0,25/sigE,length.out=101)
colbr_kdavg = seq(0,25/sigE,length.out=101)



bvec = exp(seq(log(0.0002), log(0.2), length.out=50))[(1:50)%%2 != 0]/cD
kivec = exp(seq(log(10), log(1000), length.out=25))/sigE#[-3]
rivec = exp(seq(log(0.001), log(1), length.out=25))/r
muvec = exp(seq(log(0.01), log(1), length.out=25))
chvec = exp(seq(log(0.001), log(10), length.out=25))/cD
corrlenvec = c(1.387867,1.939169,2.490472,3.041774,3.593077,4.144379,4.695682,5.246984,5.798287,6.349589,6.900892,7.452194,8.003497,8.5548,9.106102,9.657405,10.208707,10.76001,11.311312,11.862615,12.413917,12.96522,13.516522,14.067825,14.619127)*2.25/sigE

npar = 25


xlabel = "Imitation kernel size"
ylabel = "Spatial correlation length\nof spatial heterogeneity"
xaxvec = kivec
yaxvec = corrlenvec


hor_cut = 25

setEPS(width=3.76, height=4.69)
postscript("ld_diff.eps")
image.natural(x = ld_diff_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = seq(0,8,length.out=101), #colbr_lddiff,
              yfmt = "%.3g",
              ver_cut = 1,
              title = "Difference between dispersal\ndistance of milkers and killers")
dev.off()


setEPS(width=3.76, height=4.69)
postscript("ld_avg.eps")
image.natural.nonlinear(x = ld_avg_mat, xlabels = xaxvec, ylabels = yaxvec, 
                        col_data = seq(0,6,length.out=101), cols = cols, 
                        xlab = xlabel, 
                        ylab = ylabel, 
                        col_brks = seq(0,2,length.out=101), #colbr_ldavg,
                        yfmt = "%.3g",
                        title = "Average distance dispersed\nper unit time")
dev.off()

setEPS(width=3.76, height=4.69)
postscript("h_hi.eps")
image.natural(x = h_hi_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = colbr_hdiff,
              yfmt = "%.3g",
              title = "Harvesting rate\nof killers")
dev.off()

setEPS(width=3.76, height=4.69)
postscript("h_diff.eps")
image.natural(x = h_diff_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = colbr_hdiff,
              yfmt = "%.3g",
              title = "Difference between harvesting\nrates of milkers and killers")

dev.off()

setEPS(width=3.76, height=4.69)
postscript("h_avg.eps")
image.natural(x = h_avg_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = seq(0,4,length.out=101), #colbr_havg,
              yfmt = "%.3g",
              title = "Average harvesting\nrate")
dev.off()



setEPS(width=3.76, height=4.69)
postscript("kd_diff.eps")
image.natural(x = kd_diff_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = seq(0,2.5,length.out=101), #colbr_kddiff,
              yfmt = "%.3g",
              title = "Difference between dispersal\nkernel size of milkers and killers")
dev.off()


setEPS(width=3.76, height=4.69)
postscript("kd_avg.eps")
image.natural(x = kd_avg_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = seq(0,1.5,length.out=101), #
              yfmt = "%.3g",
              title = "Average dispersal\nkernel size")
dev.off()

setEPS(width=3.76, height=4.69)
postscript("rc_avg.eps")
image.natural(x = rc_avg_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = seq(0.5,1,length.out=101),
              yfmt = "%.3g",
              ver_cut=1,
              title = "Average resource\nextraction rate")

dev.off()


setEPS(width=3.76, height=4.69)
postscript("r_env.eps")
image.natural(x = r_env_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = seq(0.2,.8,length.out=101),
              yfmt = "%.3g",
              ver_cut = 1,
              title = "Resource left\nin environment")
dev.off()











########################################################################
# 
# image(t(rc_avg_mat[,1:hor_cut]), xlab="", ylab="", xaxt="no", yaxt="no", col = cols, breaks = colbr_r)
# axis(side=1, at=seq(0,length(kivec[1:hor_cut]), length.out=5)/hor_cut, labels=sprintf("%.4g", kivec[1:hor_cut])[as.integer(seq(1,hor_cut, length.out=5))])
# axis(side=2, at=seq(0,length(rivec), length.out=5)/25, labels=sprintf("%.2g", rivec)[as.integer(seq(1,length(chvec), length.out=5))])
# mtext("Imitation kernel size", side = 1, line = 3, outer = F, cex=1.2)
# mtext("Imitation rate", side = 2, line = 3, outer = F, cex=1.2)
# # mtext("Spatial correlation length of \nresource growth rate", side = 2, line = 3, outer = F, cex=1.2)
# 
# image(as.matrix(colbr_h), y=1, x=as.matrix(colbr_h), col=cols, yaxt="n", xlab="", ylab="")
# mtext("Difference between harvesting\n rates of milkers and killers", side = 1, line = 5, outer = F, cex=1.2)
# image(as.matrix(colbr_r), y=1, x=as.matrix(colbr_r), col=cols, yaxt="n", xlab="", ylab="")
# mtext("Average per capita \nresource consumed", side = 1, line = 5, outer = F, cex=1.2)


# Plot timeseries
png(filename = "~/kd_random_exp.png", width = 625*3, height = 825*3, res = 300)

# layout(rbind(c(1,1,1,1,2,2), c(3,3,3,3,4,4), c(5,5,5,5,6,6), c(7,7,7,7,8,8), c(9,9,9,9,10,10)))

id = seq(0,7500, by=10)
id[1]=1
layout(rbind(c(1,1,1,1,4,4), c(2,2,2,2,5,5), c(3,3,3,3,6,6)))
par(cex.axis=1.8, cex.lab=1.8, mar=c(4,7,1,1), oma=c(1,1,1,1))
cols = colorRampPalette(colors = c("black", "blue", "cyan", "white"))(100)
image(x = id*100*.1/5*.1, y= c(brks_h[1:nbins],brks_h[nbins]+0.01)/.2, log(dist_ts_h[id,]+1), col = cols, xlab="Time (imitation steps)", ylab="Harvesting\nrate")
image(x = id*100*.1/5*.1, y= brks_kd[1:(nbins/2.5)]/4, log(dist_ts_kd[id,1:(nbins/2.5)]+1), col = cols, xlab="Time (imitation steps)", ylab="Dispersal\nkernel size")
image(x = id*100*.1/5*.1, y= c(brks_rc[1:nbins],brks_rc[nbins]+0.01)/89, log(dist_ts_rc[id,]+1), col = cols, xlab="Time (imitation steps)", ylab="Resource\nextraction rate")
# points(x=times, y=r_tot/1.0125e7*600, type="l", col="white")
# points(x=id*100, y=filter(avgrc*6, filter=rep(1,20)/20)[id], type="l", col="green")

# par(cex.axis=1.8, cex.lab=1.8, mar=c(4,0,1,1), oma=c(1,1,1,1))
mids_h=c((brks_h[1:99]+brks_h[1:99+1])/2,2.5)
mids_kd=c((brks_kd[1:99]+brks_kd[1:99+1])/2,25)
mids_rc=c((brks_rc[1:99]+brks_rc[1:99+1])/2,600)
plot(mids_h~dist_avg_h, type="l", ylab="", yaxt="n", xlab="Frequency", col="blue", lwd=2)
plot(mids_kd[1:(nbins/2.5)]~dist_avg_kd[1:(nbins/2.5)], type="l", ylab="",  yaxt="n", xlab="Frequency", col="blue", lwd=2)
plot(mids_rc~dist_avg_rc, type="l", ylab="",  yaxt="n", xlab="Frequency", col="blue", lwd=2)

dev.off()
