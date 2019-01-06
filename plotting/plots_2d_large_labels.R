# setwd("/home/jaideep/austria_project/figures_final_data/map_muxkI_b0.00144_rI0.1")
setwd("E:/austria_project/figures/figures_2D_scans_data/map_muxkI_b0.00144_rI0.1")
# setwd("/media/jaideep/WorkData/austria_project/figures/final_manuscript/optimal")
# setwd("~/bxch_results")
setwd("E:/austria_project/figures/figures_2D_scans_data/map_rixki_b0.004")
setwd("/media/jaideep/WorkData/austria_project/figures/figures_2D_scans_data/map_bxch_costKd")

r = 0.2
K = 50
MSY = 89.3
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


h_opt_mat = read.csv(header = T, file="hoptmat.csv", row.names = 1)/r
kd_opt_mat = read.csv(header = T, file="kdoptmat.csv", row.names = 1)/sigE
rc_opt_mat = read.csv(header = T, file="rcoptmat.csv", row.names = 1)/MSY
ld_opt_mat = read.csv(header = T, file="ldoptmat.csv", row.names = 1)/r/sigE


image.natural = function(x, xlabels, ylabels, xlab, ylab, cols, col_brks, 
                         xfmt = "%.3g", yfmt = "%.1g", hor_cut=-1, ver_cut=-1,
                         title="", axlabels = NA){
  
  axsiz = 1.3
  
  layout(cbind(c(1,1,1,1,1,1,2,2)))
  par(mar=c(5.5,6,1.5,1), oma=c(5,1,1,1), cex.axis=1.8*axsiz)
  
  nlab = 3
  if(hor_cut==-1) hor_cut = dim(x)[2]
  if(ver_cut==-1) ver_cut = 1
  image(t(x[ver_cut:25,1:hor_cut]), xlab="", ylab="", xaxt="no", yaxt="no", col = cols, breaks = col_brks)
  axis(padj = .55, side=1, at=seq(0,length(xlabels[1:hor_cut]), length.out=nlab)/hor_cut, labels=sprintf(xfmt, xlabels[1:hor_cut])[as.integer(seq(1,hor_cut, length.out=nlab))])
  axis(padj = -.2, side=2, at=seq(0,length(ylabels[ver_cut:25]), length.out=nlab)/(26-ver_cut), labels=sprintf(yfmt, ylabels[ver_cut:25])[as.integer(seq(1,26-ver_cut, length.out=nlab))])
  mtext(xlab, side = 1, line = 4.7, outer = F, cex=1.2*axsiz)
  mtext(ylab, side = 2, line = 4, outer = F, cex=1.2*axsiz)
  
  image(as.matrix(col_brks), y=1, x=as.matrix(col_brks), col=cols, xaxt="n", yaxt="n", xlab="", ylab="", padj=.5)
  mtext(title, side = 1, line = 7, outer = F, cex=1.2*axsiz)
  
  if (is.na(axlabels)){
    lab_pos = as.integer(seq(1,length(col_brks), length.out=3))
    lab_names = sprintf("%.3g", col_brks[as.integer(seq(1,length(col_brks), length.out=3))])
  }
  else{
    lab_pos = NULL
    lab_names = NULL
    for (i in 1:length(axlabels)){
      lab_pos = c(lab_pos, (which(col_brks >= axlabels[i]))[1])
      lab_names = c(lab_names, sprintf("%.3g", axlabels[i]))
    }
  }
  
  axis(side=1, at=col_brks[lab_pos], 
       labels=lab_names, padj=.5 )
  
}

image.natural.nonlinear = function(x, col_data, col_labels, xlabels, ylabels, xlab, ylab, cols, col_brks, 
                                   xfmt = "%1.0g", yfmt = "%.1g", hor_cut=-1, ver_cut=-1,
                                   title=""){
  
  axsiz = 1.3
  
  layout(cbind(c(1,1,1,1,1,1,2,2)))
  par(mar=c(5.5,6,1.5,1), oma=c(5,1,1,1), cex.axis=1.8*axsiz)
  
  nlab = 3
  if(hor_cut==-1) hor_cut = dim(x)[2]
  if(ver_cut==-1) ver_cut = 1
  image(t(x[ver_cut:25,1:hor_cut]), xlab="", ylab="", xaxt="no", yaxt="no", col = cols, breaks = col_brks)
  axis(padj = .55, side=1, at=seq(0,length(xlabels[1:hor_cut]), length.out=nlab)/hor_cut, labels=sprintf(xfmt, xlabels[1:hor_cut])[as.integer(seq(1,hor_cut, length.out=nlab))])
  axis(padj = -.2, side=2, at=seq(0,length(ylabels[ver_cut:25]), length.out=nlab)/(26-ver_cut), labels=sprintf(yfmt, ylabels[ver_cut:25])[as.integer(seq(1,26-ver_cut, length.out=nlab))])
  mtext(xlab, side = 1, line = 4.7, outer = F, cex=1.2*axsiz)
  mtext(ylab, side = 2, line = 4, outer = F, cex=1.2*axsiz)
  
  image(as.matrix(col_data), y=1, x=as.matrix(col_labels), col=cols, xaxt="n", yaxt="n", xlab="", ylab="", padj=.5)
  axis(side=1, at=col_labels[as.integer(seq(1,length(col_labels), length.out=3))], 
       labels=sprintf("%.2g", col_labels[as.integer(seq(1,length(col_labels), length.out=3))]), padj=.5 )
  mtext(title, side = 1, line = 7, outer = F, cex=1.2*axsiz)
  
}

createPalette <- function(cols, values, n=100){
  nval = length(values)
  dv = values[2:nval] - values[2:nval-1]
  ncols = round(dv*n/sum(dv))
  cols1 = c()
  for (i in 1:(length(dv))){
    cols_sec = colorRampPalette(colors=c(cols[i],cols[i+1]))(ncols[i]+1)
    cols_sec = cols_sec[-length(cols_sec)]
    cols1 = c(cols1,  cols_sec)
  }
  if (sum(ncols) < n) cols1[n] = cols1[n-1]
  if (sum(ncols) > n) cols1 = cols1[1:n]
  cols1
}


col1 = c(colorRampPalette(colors=c("black","blue"))(10),
          colorRampPalette(colors=c("blue","green"))(40),
          colorRampPalette(colors=c("green","yellow"))(25),
          colorRampPalette(colors=c("yellow","red"))(25))



bvec = exp(seq(log(0.0002), log(0.2), length.out=50))[(1:50)%%2 != 0]/cD
kivec = exp(seq(log(10), log(1000), length.out=25))/sigE#[-3]
rivec = exp(seq(log(0.001), log(1), length.out=25))/r
muvec = exp(seq(log(0.01), log(1), length.out=25))
chvec = exp(seq(log(0.001), log(10), length.out=25))/cD
corrlenvec = c(1.387867,1.939169,2.490472,3.041774,3.593077,4.144379,4.695682,5.246984,5.798287,6.349589,6.900892,7.452194,8.003497,8.5548,9.106102,9.657405,10.208707,10.76001,11.311312,11.862615,12.413917,12.96522,13.516522,14.067825,14.619127)*2.25/sigE
irvvec = exp(seq(log(1),log(1000),length.out = 25))
drvvec = exp(seq(log(0.1),log(10),length.out = 25))

ir_dV95 = 1/irvvec*log(.95/(1-.95))
ir_dV5  = 1/irvvec*log(.05/(1-.05))
irwidth = ir_dV95-ir_dV5

dr_dV95 = 1/drvvec*log(.95/(1-.95))
dr_dV5  = 1/drvvec*log(.05/(1-.05))
drwidth = dr_dV95-dr_dV5

dvexpwidth = 0.5876795 - -0.6097295 
rexpwidth = 16.9478  - 0

irwidth = irwidth/dvexpwidth
drwidth = drwidth/rexpwidth


npar = 25


xlabel = "Imitation radius"
ylabel = "Imitation rate"
xaxvec = kivec
yaxvec = rivec
# xlabel = "Benefit of harvesting"
# ylabel = "Cost of harvesting"
# xaxvec = bvec
# yaxvec = chvec

xfmt = "%.3g"
yfmt = "%.2g"

hor_cut = 25



cols = colorRampPalette(colors = c("black", "red", "yellow", "white"))(100)
cols = createPalette(c("black", "blue","green3","yellow","red"),c(0,.75,2,3,12.5))
# cols = createPalette(c("black", "blue","pink","red"),c(0,.66,2.2,12.5))
image(x=seq(0,12.5,length.out=100), z=as.matrix(seq(0,1,length.out=100)), col=cols)

colbr_hdiff = seq(0,2.5/r,length.out=101)
colbr_havg = seq(0,.8/r,length.out=101)
colbr_rc = seq(.5,1,length.out=101)
colbr_ldavg = seq(0,2/r/sigE,length.out=101)
colbr_lddiff = seq(0,16/r/sigE,length.out=101)
colbr_r = seq(0,50/K,length.out=101)
colbr_kddiff = seq(0,25/sigE,length.out=101)
colbr_kdavg = seq(0,25/sigE,length.out=101)


cols = createPalette(c("black", "blue","green3","yellow","white"),c(0,.1,1,4,8))
# cols = createPalette(c("black", "blue","pink","red"),c(0,.66,2.2,12.5))
# image(x=seq(0,8,length.out=100), z=as.matrix(seq(0,1,length.out=100)), col=cols)


setEPS(width=3.49, height=5.16)
postscript("ld_diff_col.eps")
image.natural(x = ld_diff_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = seq(0,8,length.out=101), #colbr_lddiff,
              xfmt = xfmt,
              yfmt = yfmt,
              ver_cut = 1,
              title = "Difference in\ndistance dispersed",
              axlabels = c(0,0.1,1,4,8))
dev.off()

setEPS(width=3.49, height=5.16)
postscript("ld_avg_col.eps")
image.natural(x = ld_avg_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = seq(0,8,length.out=101), #colbr_lddiff,
              xfmt = xfmt,
              yfmt = yfmt,
              ver_cut = 1,
              title = "Average\ndistance dispersed",
              axlabels = c(0,0.1,1,4,8))
dev.off()

setEPS(width=3.49, height=5.16)
postscript("ld_hi_col.eps")
image.natural(x = ld_hi_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = seq(0,8,length.out=101), #colbr_lddiff,
              xfmt = xfmt,
              yfmt = yfmt,
              ver_cut = 1,
              title = "Dispersal distance of\nmobile consumers",
              axlabels = c(0,0.1,1,4,8))
dev.off()


setEPS(width=3.49, height=5.16)
postscript("ld_lo_col.eps")
image.natural(x = ld_lo_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = seq(0,8,length.out=101), #colbr_lddiff,
              xfmt = xfmt,
              yfmt = yfmt,
              ver_cut = 1,
              title = "Dispersal distance of\nsedentary consumers",
              axlabels = c(0,0.1,1,4,8))
dev.off()



cols = createPalette(c("black", "blue","green3","yellow","red"),c(0,.75,2,3,12.5))


setEPS(width=3.49, height=5.16)
postscript("h_diff_col.eps")
image.natural(x = h_diff_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = colbr_hdiff,
              yfmt = yfmt, xfmt=xfmt,
              title = "Difference in\nharvesting rate",
              axlabels = c(0,0.75, 2, 6, 12))

dev.off()

setEPS(width=3.49, height=5.16)
postscript("h_avg_col.eps")
image.natural(x = h_avg_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = seq(0,12.5,length.out=101), #colbr_havg,
              yfmt = yfmt, xfmt=xfmt,
              title = "Average\nharvesting rate",
              axlabels = c(0,0.75, 2, 6, 12))
dev.off()

setEPS(width=3.49, height=5.16)
postscript("h_hi_col.eps")
image.natural(x = h_hi_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = seq(0,12.5,length.out=101), #colbr_havg,
              xfmt = xfmt,
              yfmt = yfmt,
              title = "Harvesting rate of\nmobile consumers",
              axlabels = c(0,0.75, 2, 6, 12))
dev.off()

setEPS(width=3.49, height=5.16)
postscript("h_lo_col.eps")
image.natural(x = h_lo_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = seq(0,12.5,length.out=101), #colbr_havg,
              yfmt = yfmt, xfmt=xfmt,
              title = "Harvesting rate of\nsedentary consumers",
              axlabels = c(0,0.75, 2, 6, 12))
dev.off()


cols = createPalette(c("black", "blue","green3","yellow","white"),c(0,.1,1,2,4))
# cols = createPalette(c("black", "blue","pink","red"),c(0,.66,2.2,12.5))
# image(x=seq(0,8,length.out=100), z=as.matrix(seq(0,1,length.out=100)), col=cols)


setEPS(width=3.49, height=5.16)
postscript("kd_diff.eps")
image.natural(x = kd_diff_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = seq(0,4,length.out=101), #colbr_kddiff,
              yfmt = yfmt, xfmt=xfmt,
              title = "Difference in\ndispersal radius")
dev.off()


setEPS(width=3.49, height=5.16)
postscript("kd_avg.eps")
image.natural(x = kd_avg_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = seq(0,4,length.out=101), #
              yfmt = yfmt, xfmt=xfmt,
              title = "Average\ndispersal radius")
dev.off()

setEPS(width=3.49, height=5.16)
postscript("kd_hi.eps")
image.natural(x = kd_hi_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = seq(0,4,length.out=101), #
              yfmt = yfmt, xfmt=xfmt,
              title = "Average\ndispersal radius")
dev.off()

setEPS(width=3.49, height=5.16)
postscript("kd_lo.eps")
image.natural(x = kd_lo_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = seq(0,4,length.out=101), #
              yfmt = yfmt, xfmt=xfmt,
              title = "Average\ndispersal radius")
dev.off()



cols = createPalette(cols = c("black", "blue", "green"), values = c(0, .42, 1.004))

setEPS(width=3.49, height=5.16)
postscript("rc_avg_col.eps")
image.natural(x = rc_avg_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = seq(0,1.004,length.out=101),
              yfmt = yfmt, xfmt=xfmt,
              ver_cut=1,
              title = "Average resource\nextraction rate",
              axlabels = c(0,0.4,1))

dev.off()

cols = createPalette(cols = c("black", "orange", "yellow", "green", "darkgreen"), values = c(0, .25, .5, .75, 1))

setEPS(width=3.49, height=5.16)
postscript("r_env_col.eps")
image.natural(x = r_env_mat, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = seq(0,1,length.out=101),
              yfmt = yfmt, xfmt=xfmt,
              ver_cut = 1,
              title = "Resource left\nin environment")
dev.off()

par(mar=c(6,6,2,1))
Nc = c(128, 256, 512, 1024)/3.164
MSYNc = c(300.4, 165.7, 86.6, 43.4)/K/4/225/225/r/0.1*100
plot(MSYNc~Nc, log="xy", xlab = "Mobile consumers/nper 1000 units area", ylab = "per capita mMSY", type="o", cex.lab=1.8, cex.axis=1.8, cex=2, lwd=2)
# mtext()


setEPS(width=3.49, height=5.16)
postscript("havg_vs_opt.eps")
cols = createPalette(cols = c("blue", "white", "red"), values = c(0, 1, 3))
image.natural(x = hoptvsavg, xlabels = xaxvec, ylabels = yaxvec, 
              cols = cols, 
              xlab = xlabel, 
              ylab = ylabel, 
              col_brks = seq(0,3,length.out=101), #colbr_havg,
              yfmt = "%4.3g", xfmt="%.1g",
              title = "Average\nharvesting rate",
              axlabels = c(0,1, 2, 3))
dev.off()
