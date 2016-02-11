homdir = "/home/jaideep/austria_project/gpu_codes/output"
outdir = "LI_test"

expt = "hom"
nsteps = 750000
N = 100 #512
RT = 25
kd = 8 
h = 0.2
rI = 0.02 #bvec[ib] #0.02
L  = 100 #225
nx = 200 #450
b =  .002 #bvec[ib] # 0.0022 # 
cd = 0.1
ch = 0.08

nbins = 100
brks <- c(seq(0,2.5,length.out = nbins), 1000)

times = 1:7500
cols = rgb(colorRamp(c("black","red","yellow","white"), space = "rgb", interpolate = "spline")(0:75/75)/255)

b_imit_h = T
b_imit_kd = T
b_imit_RT = F

if (b_imit_h) h=-1
if (b_imit_kd) kd=-1
if (b_imit_RT) RT=-1
fname = sprintf("%s/%s/hist_h_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%g)",
                homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   L,      nx,   b,    cd,    ch  )

dat <- read.delim(fname, header=F)
dat <- dat[,-length(dat[1,])]
hall = as.matrix(dat)

image(x = times, y= c(brks[1:nbins],brks[nbins]+0.01), 
      log(hall+1), col = cols, 
      main=sprintf("Imit rate = %.4f",b), 
      xlab="time", ylab="Kd")

