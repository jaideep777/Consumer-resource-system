homdir = "/home/jaideep/austria_project/gpu_codes/output"
outdir = "coexist_b0.002" 

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


h0_m = .1+(0:10)*(2.5-.1)/10
h0_k = .1+(0:10)*(2.5-.1)/10
pm_mat = matrix(nrow = length(h0_m), ncol = length(h0_m), data=NA)
rc_mat = matrix(nrow = length(h0_m), ncol = length(h0_m), data=NA)

times = 1:7500

cols = rgb(colorRamp(c("black","red","yellow","white"), space = "rgb", interpolate = "spline")(0:75/75)/255)

for (i in 0:10){
for (j in 0:i){
  
  if (j<i){
  expt = "hom"
  nsteps = 750000
  N = 512
  RT = 15
  kd = 2 
  h = 0.5
  rI = .02 #bvec[ib] #0.02
  L  = 225
  nx = 450
  b =  0.002 # bvec[ib] # 0.0022 # 
  cd = 0.1
  ch = 0.08
  kI = 300
  
  b_imit_h = T
  b_imit_kd = T
  b_imit_RT = F
  
  cat(i, " ", j, "\n")
  
  if (b_imit_h) h=-1
  if (b_imit_kd) kd=-1
  if (b_imit_RT) RT=-1
  fname = sprintf("%s/%s/hist_h_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%g)_h0m(%g)_h0k(%g)",
                  homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch  , h0_m[j+1], h0_k[i+1])
  
  dat <- read.delim(fname, header=F)
  dat <- dat[,-length(dat[1,])]
  hall = as.matrix(dat)
  brks <- c(seq(0,2.5,length.out = nbins), 1000)
  
  hall_start = hall[1,]
  id_m = which(hall_start>0)[1]
  id_k = which(hall_start>0)[2]

  hall_avg = colMeans(hall[5001:7500,])
  
  pm = hall_avg[id_m]/sum(hall_avg[c(id_m,id_k)])
  
  pm_mat[i+1,j+1] = pm

  
  fname = sprintf("%s/%s/hist_rc_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%g)_h0m(%g)_h0k(%g)",
                  homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch  , h0_m[j+1], h0_k[i+1])
  
  datrc <- read.delim(fname, header=F)
  datrc <- datrc[,-length(datrc[1,])]
  hallrc = as.matrix(datrc)
  brksrc <- c(seq(0,600,length.out = nbins), 1000)
  dbrk = brksrc[2]-brksrc[1]
  mids = brksrc[1:nbins]+dbrk/2
  
  hall_avgrc = colMeans(hallrc[5001:7500,])
  
  rc_mat[i+1,j+1] = as.numeric(sum(hall_avgrc*mids)/N)
  
  #  b_mean[ib]  = mean(hmean[5001:7500])
  
  cat("\n")
  }
#  png(filename=sprintf("/home/jaideep/austria_project/gpu_codes/figures/imit_rate_ts_long/imit_%.4f.png", bvec[ib]))
#   image(x = times, y= c(brks[1:nbins],brks[nbins]+0.01), (hall+1), col = cols, xlab="time")
#  dev.off()
}
}

hall_avg = colMeans(hall[5001:7500,])
plot(hall_avg, type="l")

cols = rgb(colorRamp(c("black","red","yellow","white"), space = "rgb", interpolate = "spline")(0:75/75)/255)
# image
image(pm_mat[2:11,1:10], x=1:10, y=1:10, xaxt = "n", yaxt="n", xlab="harvesting rate of killers", ylab="harvesting rate of milkers", col=cols)
axis(side = 1,at = 1:10, labels = h0_k[2:11])
axis(side = 2,at = 1:10, labels = h0_m[1:10])
# legend
image(as.matrix(seq(0, 1, length.out=100)), y=1, x=as.matrix(seq(0, 1, length.out=100)), col=cols)


# image
image(rc_mat[2:11,1:10], x=1:10, y=1:10, xaxt = "n", yaxt="n", xlab="harvesting rate of killers", ylab="harvesting rate of milkers", col=cols)
axis(side = 1,at = 1:10, labels = h0_k[2:11])
axis(side = 2,at = 1:10, labels = h0_m[1:10])
# legend
image(as.matrix(seq(min(rc_mat, na.rm = T), max(rc_mat, na.rm = T), length.out=100)), y=1, x=as.matrix(seq(min(rc_mat, na.rm = T), max(rc_mat, na.rm = T), length.out=100)), col=cols)

