homdir = "/home/jaideep/austria_project/gpu_codes/output"
outdir = "smooth_imit_response_full" 
  
irvvec = exp(seq(log(1), log(1000), length.out=25))
bvec = exp(seq(log(0.0002), log(0.2), length.out=50))
bvec = bvec[(1:50)%%2 != 0]
kivec = exp(seq(log(10), log(1000), length.out=25))
rivec = exp(seq(log(0.001), log(10), length.out=25))
#bvec = bvec[1:20]
#bvec[1] = bvec[2]
#bvec = exp(seq(log(0.001), log(1), length.out=50))
#bvec = exp(seq(log(1), log(10), length.out=17))
#bvec = bvec[-1]
#bvec = bvec[seq(1,50,by=2)]
# bvec = exp(seq(log(0.0002), log(0.2), length.out=50))
npar = length(bvec)
nbins = 100
times = 1:7500
cols = rgb(colorRamp(c("black","red","yellow","white"), space = "rgb", interpolate = "spline")(0:75/75)/255)

h_scan = matrix(data = 0, nrow = npar, ncol=nbins)
kd_scan = matrix(data = 0, nrow = npar, ncol=nbins)
rc_scan = matrix(data = 0, nrow = npar, ncol=nbins)
r_scan = numeric(length(bvec))

b_mean = numeric(npar)

brks_h <- c(seq(0,2.5,length.out = nbins), 1000)
brks_kd <- c(seq(0,25,length.out = nbins), 1000)
brks_rc <- c(seq(0,600,length.out = nbins), 1000)

#c(1,7,13,19,25)
for (ib in 1:length(bvec)){
  expt = "hom"
  nsteps = 750000
  N = 512
  RT = 15
  kd = 2 
  h = 0.5
  rI = .02 #bvec[ib] #0.02
  L  = 225
  nx = 450
  b = 0.00589 #bvec[ib] # 0.004 # 0.0022 # as.numeric(sprintf("%.1g",bvec)[ib])
  cd = 0.1
  ch = 0.0681
  kI = 1000
  mu = .01
  irv = irvvec[ib]
  
  b_imit_h = T
  b_imit_kd = T
  b_imit_RT = F
  
  # cat(" ", ib)
  
  if (b_imit_h) h=-1
  if (b_imit_kd) kd=-1
  if (b_imit_RT) RT=-1

  # h
  if (expt == "hom"){
    fname = sprintf("%s/%s/hist_h_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%.3g)_muh(0.02)_murt(1)_mukd(0.2)_twv(20)_irv(%.3g)",
                    homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch  , irv)
  } else{
    fname = sprintf("%s/%s/hist_h_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%.3g)_tmu(%.3g)",
                    homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch  , mu)
  }
  dat <- read.delim(fname, header=F)
  dat <- dat[,-length(dat[1,])]
  hall = as.matrix(dat)
  hall_avg = colMeans(hall[5001:7500,])
  h_scan[ib,] = hall_avg
  #  b_mean[ib]  = mean(hmean[5001:7500])
  
  # kd
  if (expt == "hom"){
    fname = sprintf("%s/%s/hist_kd_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%.3g)_muh(0.02)_murt(1)_mukd(0.2)_twv(20)_irv(%.3g)",
                    homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch   , irv)
  } else{
    fname = sprintf("%s/%s/hist_kd_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%.3g)_tmu(%.3g)",
                    homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch  , mu)
  }
  dat <- read.delim(fname, header=F)
  dat <- dat[,-length(dat[1,])]
  hall = as.matrix(dat)
  hall_avg = colMeans(hall[5001:7500,])
  kd_scan[ib,] = hall_avg
  
  
  # rc
  if (expt == "hom"){
    fname = sprintf("%s/%s/hist_rc_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%.3g)_muh(0.02)_murt(1)_mukd(0.2)_twv(20)_irv(%.3g)",
                    homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch   , irv)
  } else{
    fname = sprintf("%s/%s/hist_rc_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%.3g)_tmu(%.3g)",
                    homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch  , mu)
  }
  dat <- read.delim(fname, header=F)
  dat <- dat[,-length(dat[1,])]
  hall = as.matrix(dat)
  hall_avg = colMeans(hall[5001:7500,])
  rc_scan[ib,] = hall_avg
  
  # rtotal
  if (expt == "hom"){
    fname = sprintf("%s/%s/r_total_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%.3g)_muh(0.02)_murt(1)_mukd(0.2)_twv(20)_irv(%.3g)",
                    homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch   , irv)
  } else{
    fname = sprintf("%s/%s/r_total_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%.3g)_tmu(%.3g)",
                    homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch  , mu)
  }
  dat <- read.delim(fname, header=F)
  r_avg = mean(dat$V1[5001:7500])
  r_scan[ib] = r_avg
  
  
  cat(".")
  
#  png(filename=sprintf("/home/jaideep/austria_project/gpu_codes/figures/imit_rate_ts_long/imit_%.4f.png", bvec[ib]))
#   image(x = times, y= c(brks[1:nbins],brks[nbins]+0.01), (hall+1), col = cols, xlab="time")
#  dev.off()
}
cat("\n")
# image(x = times, y= c(brks[1:nbins],brks[nbins]+0.01), log(hall+1), col = cols)

mids_rc = (brks_rc[1:nbins]+brks_rc[1:nbins+1])/2
weights = matrix(data = rep(mids_rc, npar), nrow=npar, byrow=T)
avg_rc = rowSums(weights*rc_scan)/N

cols = colorRampPalette(colors = c("white", "black"))(100)

vert_cut = 1
op=par(mfrow=c(3,1), mar = c(1,8,1,1), oma=c(4,0.,1,1.5)+0.1, mgp=c(1,1,0), cex.axis=1.8)
image(x = 1:25, y= c(brks_h[1:nbins],brks_h[nbins]+0.001)[1:(nbins*vert_cut)], log(h_scan[,1:(nbins*vert_cut)]+3), col = cols, xaxt="n", xlab="", ylab = "", cex.axis=1.8)
title(ylab="Evolved \nharvesting\n rate", line=3, cex.lab=1.8)

vert_cut=0.5
image(x = 1:25, y= c(brks_kd[1:nbins],brks_kd[nbins]+0.001)[1:(nbins*vert_cut)], log(kd_scan[,1:(nbins*vert_cut)]+3), col = cols, xaxt="n", xlab="", ylab = "", cex.axis=1.8)
title(ylab="Evolved \ndispersal\n distance", line=3, cex.lab=1.8)

vert_cut=0.5
image(x = 1:25, y= c(brks_rc[1:nbins],brks_rc[nbins]+0.001)[1:(nbins*vert_cut)], log(rc_scan[,1:(nbins*vert_cut)]+3), col = cols, xaxt="n", xlab="", ylab = "", cex.axis=1.8)
title(ylab="Evolved \nresource\n consumed", xlab="Benefit of harvesting", line=3, cex.lab=1.8)

# points(r_scan/1.0125e7*600, type="l", lwd=1, col="red")
points(avg_rc, type="l", lwd=1, col="red")

axis(side=1, at=as.integer(seq(1,length(bvec), length.out=5)), labels=sprintf("%4.1g",irvvec)[as.integer(seq(1,length(bvec), length.out=5))], cex.axis=1.8)
mtext("Slope of imitation response function", side = 1, line = 3, outer = F, cex=1.4)

abline(h=avg_rc[1]/100*600, lty=2, col="cyan4")

# find_peaks <- function (x, m = 3){
#   x = c(rep(0,m), x, rep(0,m))
#   shape <- diff(sign(diff(x, na.pad = FALSE)))
#   pks <- sapply(which(shape < 0), FUN = function(i){
#     z <- i - m + 1
#     z <- ifelse(z > 0, z, 1)
#     w <- i + m + 1
#     w <- ifelse(w < length(x), w, length(x))
#     if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
#   })
#   pks <- unlist(pks)
#   pks-m
# }
# 

# hall_avg = colMeans(hall[5001:7500,])
# plot(hall_avg, type="l")


#### R total ####
# for (ib in 1:length(bvec)){
#   expt = "hom"
#   nsteps = 750000
#   N = 512
#   RT = 15
#   kd = 2 
#   h = 0.5
#   rI = .02 #bvec[ib] #0.02
#   L  = 225
#   nx = 450
#   b =  bvec[ib] # 0.0022 # 
#   cd = 0.1
#   ch = 2
#   kI = 300
#   
#   b_imit_h = T
#   b_imit_kd = T
#   b_imit_RT = F
#   
#   cat(" ", ib)
#   
#   if (b_imit_h) h=-1
#   if (b_imit_kd) kd=-1
#   if (b_imit_RT) RT=-1
#   fname = sprintf("%s/%s/r_total_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%g)",
#                   homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch  )
#   
#   dat <- read.delim(fname, header=F)
#   
#   r_avg = mean(dat$V1[5001:7500])
#   
#   r_scan[ib] = r_avg
#   #  b_mean[ib]  = mean(hmean[5001:7500])
#   
#   cat("\n")
#   
# }
# 
# 
# 
# 
# 
# 
# 
# dat=read.delim("/home/jaideep/austria_project/gpu_codes/Consumer-resource-system/src/res200.txt", header=F)
# dat=dat[,-451]
# image(as.matrix(dat))
# which(is.na(dat))
# ########################################################################
# ########################################################################
# 
# ########################################################################
# 
# ########################################################################
# ########################################################################
# 
# 
# #### H scan ####
# 
# homdir = "/home/jaideep/austria_project/gpu_codes/output"
# outdir = "KI_test_ki300"
# 
# bvec = exp(seq(log(0.001), log(1), length.out=50))
# # bvec = exp(seq(log(0.0002), log(0.2), length.out=50))
# npar = length(bvec)
# nbins = 100
# brks <- c(seq(0,2.5,length.out = nbins), 1000)
# b_scan = matrix(data = 0, nrow = npar, ncol=nbins)
# b_mean = numeric(npar)
# 
# for (ib in 1:length(bvec)){
#   expt = "hom"
#   nsteps = 250000
#   N = 512
#   RT = 15
#   kd = 2 
#   h = 0.1
#   rI = .02 #bvec[ib] #0.02
#   kI = 300
#   L  = 225
#   nx = 450
#   b = bvec[ib] # 0.0022 # 
#   cd = 0.1
#   ch = 0.08
#   
#   b_imit_h = T
#   b_imit_kd = T
#   b_imit_RT = F
#   
#   
#   if (b_imit_h) h=-1
#   if (b_imit_kd) kd=-1
#   if (b_imit_RT) RT=-1
#   fname = sprintf("%s/%s/ldc_%s_T(%g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%g)",
#                   homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,  kI,  L,      nx,   b,    cd,    ch  )
#   
#   dat <- read.delim(fname, header=F)
#   dat <- dat[,-length(dat[1,])]
#   
#   nt = nrow(dat)
#   times = dat[,1]
#   hc <- dat[,-1]
#   
#   cat(" ", ib)
#   
#   hall = matrix(data = 0, nrow = nt, ncol=nbins)
#   hmean = numeric(nt)
#   for (i in 1:nt){
#     hist <- hist(as.numeric(hc[i,]), breaks = brks, plot = F)
#     hall[i,] = hist$counts
#     hmean[i] = mean(as.numeric(dat[i,-1]))
#     if (i %% 1000 == 0) cat(".")
#   }
# 
#   hall_avg = colMeans(hall[5001:7500,])
#   
#   b_scan[ib,] = hall_avg
#   b_mean[ib]  = mean(hmean[5001:7500])
#   
#   cat("\n")
# 
# }
# cols = rgb(colorRamp(c("black","red","yellow","white"), space = "rgb", interpolate = "spline")(0:75/75)/255)
# 
# # image(x = times, y= c(brks[1:nbins],brks[nbins]+0.01), log(hall+1), col = cols)
# image(x = 1:50, y= c(brks[1:nbins],brks[nbins]+0.01), log(b_scan+3), col = cols, xaxt="n", xlab="Imitation rate", ylab = "evolved h")
# axis(side=1, at=as.integer(seq(1,50, length.out=5)), labels=format(bvec, digits=2)[as.integer(seq(1,50, length.out=5))])
# points(b_mean~seq(1,50,by=1), col="cyan", type="l", lwd=3)
# 
# write.csv(b_scan, "b_scan_rimit_N512_set2.csv")
# # plot(hmean, type="l")
# 
# 
# 
# #### Rc scan ####
# 
# 
# homdir = "/home/jaideep/austria_project/gpu_codes/output"
# outdir = "rimit_scan_large_b0.02_sync"
# 
# bvec = exp(seq(log(0.001), log(1), length.out=50))
# # bvec = exp(seq(log(0.0002), log(0.2), length.out=50))
# npar = length(bvec)
# nbins = 100
# brks <- c(seq(0,600,length.out = nbins), 1000000)
# rc_scan = matrix(data = 0, nrow = npar, ncol=nbins)
# rc_mean = numeric(npar)
# 
# for (ib in 1:length(bvec)){
#   expt = "hom"
#   nsteps = 750000
#   N = 512
#   RT = 25
#   kd = 8 
#   h = 0.2
#   rI = bvec[ib] # 0.02
#   L  = 225
#   nx = L*2
#   b = 0.002 # 
#   cd = 0.1
#   ch = 0.08
#   
#   b_imit_h = T
#   b_imit_kd = T
#   b_imit_RT = F
#   
#   
#   if (b_imit_h) h=-1
#   if (b_imit_kd) kd=-1
#   if (b_imit_RT) RT=-1
#   fname = sprintf("%s/%s/rc_%s_T(%g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%g)",
#                   homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   L,      nx,   b,    cd,    ch  )
#   
#   dat <- read.delim(fname, header=F)
#   dat <- dat[,-length(dat[1,])]
#   
#   nt = nrow(dat)
#   times = dat[,1]
#   hc <- dat[,-1]
#   
#   cat(" ", ib)
#   
#   hall = matrix(data = 0, nrow = nt, ncol=nbins)
#   hmean = numeric(nt)
#   for (i in 1:nt){
#     hist <- hist(as.numeric(hc[i,]), breaks = brks, plot = F)
#     hall[i,] = hist$counts
#     hmean[i] = mean(as.numeric(dat[i,-1]))
#     if (i %% 1000 == 0) cat(".")
#   }
#   
#   hall_avg = colMeans(hall[5001:7500,])
#   
#   rc_scan[ib,] = hall_avg
#   rc_mean[ib]  = mean(hmean[5001:7500])
#   
#   cat("\n")
#   
# }
# cols = rgb(colorRamp(c("black","red","yellow","white"), space = "rgb", interpolate = "spline")(0:75/75)/255)
# 
# # image(x = times, y= c(brks[1:nbins],brks[nbins]+0.01), log(hall+1))
# image(x = 1:50, y= c(brks[1:nbins],brks[nbins]+0.01), log(rc_scan+1), col = cols, xaxt="n", xlab="Imitation rate", ylab = "Resource consumed")
# axis(side=1, at=as.integer(seq(1,50, length.out=5)), labels=format(bvec, digits=2)[as.integer(seq(1,50, length.out=5))])
# points(rc_mean~seq(1,50,by=1), col="white", type="l", lwd=3)
# 
# write.csv(rc_scan, "rc_scan_rimit_N512.csv")
# # plot(hmean, type="l")
# 
# 
# 
# #### KD scan ####
# 
# 
# homdir = "/home/jaideep/austria_project/gpu_codes/output"
# outdir = "rimit_scan_large_b0.02_sync"
# 
# bvec = exp(seq(log(0.001), log(1), length.out=50))
# # bvec = exp(seq(log(0.0002), log(0.2), length.out=50))
# npar = length(bvec)
# nbins = 100
# brks <- c(seq(0,16,length.out = nbins), 1000000)
# kd_scan = matrix(data = 0, nrow = npar, ncol=nbins)
# kd_mean = numeric(npar)
# 
# for (ib in 1:length(bvec)){
#   expt = "hom"
#   nsteps = 750000
#   N = 512
#   RT = 25
#   kd = 8 
#   h = 0.2
#   rI = bvec[ib] # 0.02
#   L  = 225
#   nx = L*2
#   b =  0.002 # 
#   cd = 0.1
#   ch = 0.08
#   
#   b_imit_h = T
#   b_imit_kd = T
#   b_imit_RT = F
#   
#   
#   if (b_imit_h) h=-1
#   if (b_imit_kd) kd=-1
#   if (b_imit_RT) RT=-1
#   fname = sprintf("%s/%s/kd_%s_T(%g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%g)",
#                   homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   L,      nx,   b,    cd,    ch  )
#   
#   dat <- read.delim(fname, header=F)
#   dat <- dat[,-length(dat[1,])]
#   
#   nt = nrow(dat)
#   times = dat[,1]
#   hc <- dat[,-1]
#   
#   cat(" ", ib)
#   
#   hall = matrix(data = 0, nrow = nt, ncol=nbins)
#   hmean = numeric(nt)
#   for (i in 1:nt){
#     hist <- hist(as.numeric(hc[i,]), breaks = brks, plot = F)
#     hall[i,] = hist$counts
#     hmean[i] = mean(as.numeric(dat[i,-1]))
#     if (i %% 1000 == 0) cat(".")
#   }
#   
#   hall_avg = colMeans(hall[5001:7500,])
#   
#   kd_scan[ib,] = hall_avg
#   kd_mean[ib]  = mean(hmean[5001:7500])
#   
#   cat("\n")
#   
# }
# cols = rgb(colorRamp(c("black","red","yellow","white"), space = "rgb", interpolate = "spline")(0:75/75)/255)
# 
# # image(x = times, y= c(brks[1:nbins],brks[nbins]+0.01), log(hall+1), col=cols)
# image(x = 1:50, y= c(brks[1:nbins],brks[nbins]+0.01), log(kd_scan+2), col = cols, xaxt="n", xlab="Imitation rate", ylab = "Resource consumed")
# axis(side=1, at=as.integer(seq(1,50, length.out=5)), labels=format(bvec, digits=2)[as.integer(seq(1,50, length.out=5))])
# 
# write.csv(kd_scan, "kd_scan_N512.csv")
# # plot(hmean, type="l")
# 
# points(kd_mean~seq(1,50,by=1), type="l", col="white", lwd=3)
# 
# 
# 


