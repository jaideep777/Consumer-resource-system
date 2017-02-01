homdir = "/home/jaideep/austria_project/gpu_codes/output"
outdir = "msy_26x26" 

# bvec = exp(seq(log(0.0002), log(0.2), length.out=50))
# bvec = bvec[(1:50)%%2 != 0]
kivec = exp(seq(log(10), log(1000), length.out=25))[-3]
rivec = exp(seq(log(0.001), log(1), length.out=25))
#muvec = exp(seq(log(0.01), log(1), length.out=25))
#bvec = bvec[1:20]
#bvec[1] = bvec[2]
#bvec = exp(seq(log(0.001), log(1), length.out=50))
#bvec = exp(seq(log(1), log(10), length.out=17))
#bvec = bvec[-1]
#bvec = bvec[seq(1,50,by=2)]
# bvec = exp(seq(log(0.0002), log(0.2), length.out=50))

# hvec = seq(0,2.5, length.out=11) #c(0.1,0.15,0.2,0.25,0.3)
# kdvec = seq(0,10, length.out=11)
hvec = seq(0,1.25, length.out=26) #c(0.1,0.15,0.2,0.25,0.3)
kdvec = seq(0,5, length.out=26)
npar = 26

rcmat = matrix(nrow=npar, ncol=npar)
rmat = matrix(nrow=npar, ncol=npar)
for (ir in 1:npar){
  for (ib in 1:npar){
    expt = "hom"
    nsteps = 7500
    N = 512
    RT = 15
    kd = kdvec[ir]
    h =  hvec[ib]
    rI = 0 #rivec[ir] #0.02
    L  = 225
    nx = 450
    b =  0 #bvec[ib] # 0.004 # 0.0022 # 
    cd = 0
    ch = 0 #rivec[ir] #0.08
    kI = 0 #kivec[ib] #1000
    mu = 0 #muvec[ir]
    
    b_imit_h = F
    b_imit_kd = F
    b_imit_RT = F
    
    cat(".")
    
    if (b_imit_h) h=-1
    if (b_imit_kd) kd=-1
    if (b_imit_RT) RT=-1
    
    # rtotal
    #     fname = sprintf("%s/%s/r_total_%s_T(%.3g)_N(%g)_RT(%g)_kd(%g)_h(%g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%.3g)_tmu(%.3g)",
    #                     homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch  , mu)
    fname = sprintf("%s/%s/r_total_%s_T(%.3g)_N(%g)_RT(%g)_kd(%.3g)_h(%.3g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%.3g)",
                    homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch  )
    dat1 <- read.delim(fname, header=F)
    
    fname = sprintf("%s/%s/rc_avg_%s_T(%.3g)_N(%g)_RT(%g)_kd(%.3g)_h(%.3g)_rI(%.3g)_kI(%.3g)_L(%g)_nx(%g)_b(%.3g)_cd(%g)_ch(%.3g)",
                    homdir,outdir,  expt, nsteps/1000,  N, RT, kd, h,     rI,   kI, L,      nx,   b,    cd,    ch  )
    dat2 <- read.delim(fname, header=F)
    
    # plot(dat2$V1, type="l", ylim=c(0,300))
    # points(dat1$V1*50/1.0125e7, type="l", col="red")
    # abline(h=50, lty=2, col="cyan4")
    
    arl = nsteps/5
    start = arl*0.8+1
    end = arl
    rmat[ir,ib] = mean(dat1$V1[start:end])/1.0125e7*50
    rcmat[ir,ib] = mean(dat2$V1[start:end])
    
    cat(".") 
  }  
  cat("\n")
}

cols = colorRampPalette(colors = c("black", "red", "yellow", "white"))(100)
colbr_rc = seq(0,100,length.out=101)
colbr_r = seq(0,50,length.out=101)

layout(rbind(c(1,2),c(1,2),c(1,2),c(1,2),c(3,4)))
par(mar=c(4,6,1,1), oma=c(5,1,1,1), cex.axis=1.8)
image(t(rcmat), xlab="", ylab="", xaxt="no", yaxt="no", col = cols, breaks = colbr_rc)
axis(side=1, at=seq(0,1, length.out=5), labels=sprintf("%.4g", hvec)[as.integer(seq(1,length(hvec), length.out=5))])
axis(side=2, at=seq(0,1, length.out=5), labels=sprintf("%.4g", kdvec)[as.integer(seq(1,length(kdvec), length.out=5))])
mtext("Harvesting rate", side = 1, line = 3, outer = F, cex=1.3)
mtext("Dispersal distance", side = 2, line = 3, outer = F, cex=1.3)
# mtext("Spatial correlation length of \nresource growth rate", side = 2, line = 3, outer = F, cex=1.2)

image(t(rmat), xlab="", ylab="", xaxt="no", yaxt="no", col = cols, breaks = colbr_r)
axis(side=1, at=seq(0,1, length.out=5), labels=sprintf("%.4g", hvec)[as.integer(seq(1,length(hvec), length.out=5))])
axis(side=2, at=seq(0,1, length.out=5), labels=sprintf("%.4g", kdvec)[as.integer(seq(1,length(kdvec), length.out=5))])
mtext("Harvesting rate", side = 1, line = 3, outer = F, cex=1.3)
mtext("Dispersal distance", side = 2, line = 3, outer = F, cex=1.3)
# mtext("Spatial correlation length of \nresource growth rate", side = 2, line = 3, outer = F, cex=1.2)

image(as.matrix(colbr_rc), y=1, x=as.matrix(colbr_rc), col=cols, yaxt="n", xlab="", ylab="")
mtext("Average resource\n consumed", side = 1, line = 5, outer = F, cex=1.3)
image(as.matrix(colbr_r), y=1, x=as.matrix(colbr_r), col=cols, yaxt="n", xlab="", ylab="")
mtext("Average resource\n left", side = 1, line = 5, outer = F, cex=1.3)


