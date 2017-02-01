homdir = "/home/jaideep/austria_project/gpu_codes/output"
outdir = "smooth_imit_response" 

expt = "hom"
nsteps = 750000
N = 512
RT = 15
kd = 2 
h = 0.5
rI = .02 #bvec[ib] #0.02
L  = 225
nx = 450
b = 0.00589 # 0.0022 # as.numeric(sprintf("%.1g",bvec)[ib])
cd = 0.1
ch = 0.0681
kI = 1000
mu = .01
irv = 1

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
                  homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch  ,irv)
} else{
  fname = sprintf("%s/%s/hist_h_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%.3g)_tmu(%.3g)",
                  homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch  , mu)
}
dat <- read.delim(fname, header=F)
dat <- dat[,-length(dat[1,])]
hall_h = as.matrix(dat)
#  b_mean[ib]  = mean(hmean[5001:7500])

# kd
if (expt == "hom"){
  fname = sprintf("%s/%s/hist_kd_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%.3g)_muh(0.02)_murt(1)_mukd(0.2)_twv(20)_irv(%.3g)",
                  homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch  ,irv)
} else{
  fname = sprintf("%s/%s/hist_kd_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%.3g)_tmu(%.3g)",
                  homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch  , mu)
}
dat <- read.delim(fname, header=F)
dat <- dat[,-length(dat[1,])]
hall_k = as.matrix(dat)


# rc
if (expt == "hom"){
  fname = sprintf("%s/%s/hist_rc_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%.3g)_muh(0.02)_murt(1)_mukd(0.2)_twv(20)_irv(%.3g)",
                  homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch  ,irv)
} else{
  fname = sprintf("%s/%s/hist_rc_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%.3g)_tmu(%.3g)",
                  homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch  , mu)
}
dat <- read.delim(fname, header=F)
dat <- dat[,-length(dat[1,])]
hall_r = as.matrix(dat)

# rtotal
if (expt == "hom"){
  fname = sprintf("%s/%s/r_total_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%.3g)_muh(0.02)_murt(1)_mukd(0.2)_twv(20)_irv(%.3g)",
                  homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch  ,irv)
} else{
  fname = sprintf("%s/%s/r_total_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%.3g)_tmu(%.3g)",
                  homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch  , mu)
}
dat <- read.delim(fname, header=F)
r_avg = mean(dat$V1[5001:7500])


id = seq(0,nsteps/100, by=10)
id[1]=1
# layout(rbind(c(1,3),c(2,4)))
layout(rbind(c(1,1),c(2,2),c(3,3)))
par(cex.axis=1.8, cex.lab=1.8, mar=c(4,7,1,1), oma=c(1,1,1,1))
cols = colorRampPalette(colors = c("black", "blue", "cyan", "white"))(100)
image(x = id*100/1e5, y= c(brks_h[1:nbins],brks_h[nbins]+0.01), log(hall_h[id,]+1), col = cols, xlab="Time", ylab="Harvesting \n rate")
image(x = id*100/1e5, y= brks_kd[1:(nbins/2.5)], log(hall_k[id,1:(nbins/2.5)]+1), col = cols, xlab="Time", ylab="Dispersal \n distance")
image(x = id*100/1e5, y= brks_rc[1:(nbins/2.5)], log(hall_r[id,1:(nbins/2.5)]+1), col = cols, xlab="Time", ylab="Dispersal \n distance")

