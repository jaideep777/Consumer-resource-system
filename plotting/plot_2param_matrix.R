find_peaks <- function (x, m = 3){
  x = c(rep(0,m), x, rep(0,m))
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks-m
}

find_peaks_top2 = function(x, m){
  pks = find_peaks(x=x, m=m)
  if (length(pks) < 2){
    lo = pks
    hi = pks
  } else if (length(pks) == 2){
    lo = pks[1]
    hi = pks[2]
  } else {
    pksid = order(x[pks], decreasing = T)[1:2]
    lo = pks[pksid[1]]
    hi = pks[pksid[2]]
  }
  c(min(hi,lo), max(hi,lo))
}

homdir = "/home/jaideep/austria_project/gpu_codes/output"
outdir = "map_muxkI_b0.00144_rI0.1" 

bvec = exp(seq(log(0.0002), log(0.2), length.out=50))
bvec = bvec[(1:50)%%2 != 0]
kivec = exp(seq(log(10), log(1000), length.out=25))#[-3]
rivec = exp(seq(log(0.001), log(1), length.out=25))
muvec = exp(seq(log(0.01), log(1), length.out=25))
chvec = exp(seq(log(0.001), log(10), length.out=25))
  
npar = 25
nbins = 100
times = 1:7500

brks_h <- c(seq(0,2.5,length.out = nbins), 1000)
brks_kd <- c(seq(0,25,length.out = nbins), 1000)
brks_rc <- c(seq(0,600,length.out = nbins), 1000)
brks_ldc <- c(seq(0,50,length.out = nbins), 1000)

mids_h = brks_h[1:nbins]
mids_rc = (brks_rc[1:nbins]+brks_rc[1:nbins+1])/2
mids_ldc = brks_ldc[1:nbins] #(brks_ldc[1:nbins]+brks_ldc[1:nbins+1])/2
mids_kd = brks_kd[1:nbins]

weights_h = matrix(data = rep(mids_h, npar), nrow=npar, byrow=T)
weights_rc = matrix(data = rep(mids_rc, npar), nrow=npar, byrow=T)
weights_ldc = matrix(data = rep(mids_ldc, npar), nrow=npar, byrow=T)
weights_kd = matrix(data = rep(mids_kd, npar), nrow=npar, byrow=T)

h_lo_mat = matrix(nrow=npar, ncol=npar)
h_hi_mat = matrix(nrow=npar, ncol=npar)
h_avg_mat = matrix(nrow=npar, ncol=npar)

kd_lo_mat = matrix(nrow=npar, ncol=npar)
kd_hi_mat = matrix(nrow=npar, ncol=npar)
kd_avg_mat = matrix(nrow=npar, ncol=npar)

ld_lo_mat = matrix(nrow=npar, ncol=npar)
ld_hi_mat = matrix(nrow=npar, ncol=npar)
ld_avg_mat = matrix(nrow=npar, ncol=npar)

r_env_mat = matrix(nrow=npar, ncol=npar)
rc_avg_mat = matrix(nrow=npar, ncol=npar)

for (iy in 1:npar){
  
  h_scan = matrix(data = 0, nrow = npar, ncol=nbins)
  kd_scan = matrix(data = 0, nrow = npar, ncol=nbins)
  rc_scan = matrix(data = 0, nrow = npar, ncol=nbins)
  ldc_scan = matrix(data = 0, nrow = npar, ncol=nbins)
  r_scan = numeric(length(bvec))

  ##### x scan #####
  for (ix in 1:length(bvec)){
    expt = "het"
    nsteps = 750000
    N = 512
    RT = 15
    kd = 2 
    h = 0.5
    rI = 0.1 #rivec[iy] #.02 #bvec[ib] #0.02
    L  = 225
    nx = 450
    b = 0.00144 #bvec[ix] # 0.004 # 0.0022 # as.numeric(sprintf("%.1g",bvec)[ib])
    cd = 0.1
    ch = 0.08 #chvec[iy]
    kI = kivec[ix]
    mu = muvec[iy] #.01
    
    b_imit_h = T
    b_imit_kd = T
    b_imit_RT = F
    
    # cat(" ", ix)
    
    if (b_imit_h) h=-1
    if (b_imit_kd) kd=-1
    if (b_imit_RT) RT=-1
    
    expt_params = sprintf("%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%.3g)", #_muh(0.02)_murt(1)_mukd(0.2)_twv(20)", )
                          expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch  )
    
    if (expt == "het") expt_params = paste(expt_params, sprintf("_tmu(%.3g)", mu), sep="")
    
    
    # h
    fname = paste(homedir,"/",outdir,"/","hist_h_", expt_params, sep="")
    dat <- read.delim(fname, header=F)
    dat <- dat[,-length(dat[1,])]
    dist_ts = as.matrix(dat)
    dist_avg = colMeans(dist_ts[5001:7500,])
    h_scan[ix,] = dist_avg
    
    # kd
    fname = paste(homedir,"/",outdir,"/","hist_kd_", expt_params, sep="")
    dat <- read.delim(fname, header=F)
    dat <- dat[,-length(dat[1,])]
    dist_ts = as.matrix(dat)
    dist_avg = colMeans(dist_ts[5001:7500,])
    kd_scan[ix,] = dist_avg
    
    # rc
    fname = paste(homedir,"/",outdir,"/","hist_rc_", expt_params, sep="")
    dat <- read.delim(fname, header=F)
    dat <- dat[,-length(dat[1,])]
    dist_ts = as.matrix(dat)
    dist_avg = colMeans(dist_ts[5001:7500,])
    rc_scan[ix,] = dist_avg
    
    # ldc
    fname = paste(homedir,"/",outdir,"/","hist_ldc_", expt_params, sep="")
    dat <- read.delim(fname, header=F)
    dat <- dat[,-length(dat[1,])]
    dist_ts = as.matrix(dat)
    dist_avg = colMeans(dist_ts[5001:7500,])
    ldc_scan[ix,] = dist_avg
    
    # r_total
    fname = paste(homedir,"/",outdir,"/","r_total_", expt_params, sep="")
    dat <- read.delim(fname, header=F)
    r_avg = mean(dat$V1[5001:7500])
    r_scan[ix] = r_avg
    
    cat(".")
  }
  
  h_hi = numeric(npar)
  h_lo = numeric(npar)
  kd_hi = numeric(npar)
  kd_lo = numeric(npar)
  ldc_hi = numeric(npar)
  ldc_lo = numeric(npar)
  
  for (i in 1:npar){
    pks = find_peaks_top2(x=ldc_scan[i,], m=1)
    ldc_hi[i] = brks_ldc[pks[2]]
    ldc_lo[i] = brks_ldc[pks[1]]
    
    pks = find_peaks_top2(x=kd_scan[i,], m=3)
    kd_hi[i] = brks_kd[pks[2]]
    kd_lo[i] = brks_kd[pks[1]]
    
    pks = find_peaks_top2(x=h_scan[i,], m=5)
    h_hi[i] = brks_h[pks[2]]
    h_lo[i] = brks_h[pks[1]]
    
  }
  
  avg_h = rowSums(weights_h*h_scan)/N
  avg_rc = rowSums(weights_rc*rc_scan)/N
  avg_ldc = rowSums(weights_ldc*ldc_scan)/N
  avg_kd = rowSums(weights_kd*kd_scan)/N
  
  # write to matrices
  
  h_lo_mat[iy,] = h_lo
  h_hi_mat[iy,] = h_hi
  h_avg_mat[iy,] = avg_h
  
  kd_lo_mat[iy,] = kd_lo
  kd_hi_mat[iy,] = kd_hi
  kd_avg_mat[iy,] = avg_kd

  ld_lo_mat[iy,] = ldc_lo
  ld_hi_mat[iy,] = ldc_hi
  ld_avg_mat[iy,] = avg_ldc
  
  r_env_mat[iy,] = r_scan
  rc_avg_mat[iy,] = avg_rc
  
  cat("\n")
  
  ##### #####
  
}

write.csv(x = h_lo_mat, file="h_lo.csv")
write.csv(x = h_hi_mat, file="h_hi.csv")
write.csv(x = h_avg_mat, file="h_avg.csv")

write.csv(x = kd_lo_mat, file="kd_lo.csv")
write.csv(x = kd_hi_mat, file="kd_hi.csv")
write.csv(x = kd_avg_mat, file="kd_avg.csv")

write.csv(x = ld_lo_mat, file="ld_lo.csv")
write.csv(x = ld_hi_mat, file="ld_hi.csv")
write.csv(x = ld_avg_mat, file="ld_avg.csv")

write.csv(x = rc_avg_mat, file="rc_avg.csv")
write.csv(x = r_env_mat, file="r_env.csv")

# op=par(mfrow=c(3,1), mar = c(1,5,1,1) + 0.1, oma=c(4,0.,1,1.5)+0.1, mgp=c(1,1,0))
# # par(mfrow=c(3,1), mgp=c(1,1,0))
# image(x = 1:25, y= c(brks_h[1:nbins],brks_h[nbins]+0.001), log(h_scan+3), col = cols, xaxt="n", xlab="", ylab = "", cex.axis=1.5)
# title(ylab="Evolved h", line=3, cex.lab=1.5)
# 
# image(x = 1:25, y= c(brks_kd[1:nbins],brks_kd[nbins]+0.001), log(kd_scan+3), col = cols, xaxt="n", xlab="", ylab = "", cex.axis=1.5)
# title(ylab="Evolved kd", line=3, cex.lab=1.5)
# 
# image(x = 1:25, y= c(brks_rc[1:nbins],brks_rc[nbins]+0.001), log(rc_scan+3), col = cols, xaxt="n", xlab="", ylab = "", cex.axis=1.5)
# title(ylab="Evolved rc", xlab="Benefit of harvesting", line=3, cex.lab=1.5)
# points(r_scan/1e7*500, type="l", lwd=1, col="white")
# 
# axis(side=1, at=as.integer(seq(1,length(bvec), length.out=5)), labels=format(bvec, digits=1, scientific = F)[as.integer(seq(1,length(bvec), length.out=5))], cex.axis=1.5)
# mtext("Benefit of harvesting", side = 1, line = 3, outer = F, cex=1.)
# 

cols = rgb(colorRamp(c("black","red","yellow","white"), space = "rgb", interpolate = "spline")(0:75/75)/255)
image(t((h_hi_mat-h_lo_mat)), xlab="", ylab="", xaxt="no", yaxt="no", col=cols)
image(t((rc_avg_mat)), xlab="", ylab="", xaxt="no", yaxt="no", col=cols)
image(t((ld_avg_mat)), xlab="", ylab="", xaxt="no", yaxt="no", col=cols)
image(t((ld_hi_mat-ld_lo_mat)), xlab="", ylab="", xaxt="no", yaxt="no", col=cols)
# image(t((h_lo_mat)), xlab="", ylab="", xaxt="no", yaxt="no", col=cols)

# par(mfrow=c(1,1))
# # filled.contour(t(hdiffmat), xlab="", yab="", xaxt="no", yaxt="no", col = rgb(colorRamp(c("black","red","yellow","white"), space = "rgb", interpolate = "spline")(0:25/25)/255))
# image(t(rcmat), xlab="", ylab="", xaxt="no", yaxt="no", col=cols)
# axis(side=1, at=seq(1,length(kivec), length.out=5)/25, labels=format(kivec, digits=1, scientific = F)[as.integer(seq(1,length(kivec), length.out=5))], cex.axis=1.2)
# axis(side=2, at=seq(1,length(rivec), length.out=5)/25, labels=format(rivec, digits=1, scientific = F)[as.integer(seq(1,length(rivec), length.out=5))], cex.axis=1.2)
# mtext("Imitation kernel size", side = 1, line = 3, outer = F, cex=1.2)
# mtext("Imitation rate", side = 2, line = 3, outer = F, cex=1.2)
# 
# image(as.matrix(seq(min(rcmat), max(rcmat), length.out=10)), y=1, x=as.matrix(seq(min(rcmat), max(rcmat), length.out=10)), col=cols)
# image(as.matrix(seq(min(hdiffmat), max(hdiffmat), length.out=10)), y=1, x=as.matrix(seq(min(hdiffmat), max(hdiffmat), length.out=10)), col=cols)
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


