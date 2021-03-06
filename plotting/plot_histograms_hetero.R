homdir = "/home/jaideep/austria_project/gpu_codes/output"
outdir = "LI_lowHet_b_scan_RT15_rI0.1"

bvec = exp(seq(log(0.0002), log(.2), length.out=50))
#  bvec = exp(seq(log(0.001), log(1), length.out=50))
npar = length(bvec)
nbins = 100
b_scan = matrix(data = 0, nrow = npar, ncol=nbins)
b_mean = numeric(npar)

times = 1:7500

cols = rgb(colorRamp(c("black","red","yellow","white"), space = "rgb", interpolate = "spline")(0:75/75)/255)

expt = "het"
nsteps = 750000
N = 512
RT = 15
kd = 8 
h = 0.2
L  = 225
nx = 450
cd = 0.1

brks <- c(seq(0,2.5,length.out = nbins), 1000)
for (ib in 1:length(bvec)){
  rI = 0.1 #bvec[ib] #0.02
  b =  bvec[ib] # 0.0022 # 
  mu = 1 # bvec[ib]
  ch = 0.08
  
  b_imit_h = T
  b_imit_kd = T
  b_imit_RT = F
  
  cat(" ", ib)
  
  if (b_imit_h) h=-1
  if (b_imit_kd) kd=-1
  if (b_imit_RT) RT=-1
  fname = sprintf("%s/%s/hist_h_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%g)_tmu(%.3g)",
                  homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   L,      nx,   b,    cd,    ch  , mu)
  
  dat <- read.delim(fname, header=F)
  dat <- dat[,-length(dat[1,])]
  hall = as.matrix(dat)
  
  hall_avg = colMeans(hall[5001:7500,])
  
  b_scan[ib,] = hall_avg
  #  b_mean[ib]  = mean(hmean[5001:7500])
  
  cat("\n")
  
  #  png(filename=sprintf("/home/jaideep/austria_project/gpu_codes/figures/imit_rate_ts_long/imit_%.4f.png", bvec[ib]))
  #  image(x = times, y= c(brks[1:nbins],brks[nbins]+0.01), log(hall+1), col = cols, main=sprintf("Imit rate = %.4f",bvec[ib]), xlab="time", ylab="Kd")
  #  dev.off()
}
# image(x = times, y= c(brks[1:nbins],brks[nbins]+0.01), log(hall+1), col = cols)
image(x = 1:50, y= c(brks[1:nbins],brks[nbins]+0.01), log(b_scan+3), col = cols, xaxt="n", xlab="Benefit of harvesting", ylab = "evolved h")
axis(side=1, at=as.integer(seq(1,length(bvec), length.out=5)), labels=format(bvec, digits=2)[as.integer(seq(1,length(bvec), length.out=5))])


hall_avg = colMeans(hall[5001:7500,])
plot(hall_avg, type="l")

