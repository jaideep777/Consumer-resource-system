for (i in 1:npar){
  pks = find_peaks(x=ldc_scan[i,], m=1)
  if (length(pks) == 1) pks = c(pks,pks)
  pks = sort(pks)[1:2]
  ldc_hi[i] = brks_ldc[pks[2]]
  ldc_lo[i] = brks_ldc[pks[1]]
  
  pks = find_peaks(x=h_scan[i,], m=5)
  if (length(pks) == 1) pks = c(pks,pks)
  pks = sort(pks)[1:2]
  h_hi[i] = brks_h[pks[2]]
  h_lo[i] = brks_h[pks[1]]
  
}

# kd_scan = kd_scan*t(matrix(data = rep(1*(ldc_hi>0), nbins), nrow=nbins, byrow=T))
# kd_scan[,1] = kd_scan[,1] + 1*(ldc_hi==0)*N

for (i in 1:npar){
  pks = find_peaks(x=kd_scan[i,], m=3)
  if (length(pks) == 1) pks = c(pks,pks)
  pks = sort(pks)[1:2]
  kd_hi[i] = brks_kd[pks[2]]
  kd_lo[i] = brks_kd[pks[1]]
}



mids_rc = (brks_rc[1:nbins]+brks_rc[1:nbins+1])/2
weights = matrix(data = rep(mids_rc, npar), nrow=npar, byrow=T)
avg_rc = rowSums(weights*rc_scan)/N

mids_ldc = brks_ldc[1:nbins] #(brks_ldc[1:nbins]+brks_ldc[1:nbins+1])/2
weights = matrix(data = rep(mids_ldc, npar), nrow=npar, byrow=T)
avg_ldc = rowSums(weights*ldc_scan)/N


cols = colorRampPalette(colors = c("white", "black"))(100)

vert_cut = 1
op=par(mfrow=c(3,1), mar = c(1,8,1,1), oma=c(4,0.,1,1.5)+0.1, mgp=c(1,1,0), cex.axis=1.8)
image(x = 1:25, y= c(brks_h[1:nbins],brks_h[nbins]+0.001)[1:(nbins*vert_cut)], log(h_scan[,1:(nbins*vert_cut)]+3), col = cols, xaxt="n", xlab="", ylab = "", cex.axis=1.8)
title(ylab="Evolved \nharvesting\n rate", line=3, cex.lab=1.8)
points(h_hi, type="l", col="red")
points(h_lo, type="l", col="green")

vert_cut=0.5
image(x = 1:25, y= c(brks_kd[1:nbins],brks_kd[nbins]+0.001)[1:(nbins*vert_cut)], log(kd_scan[,1:(nbins*vert_cut)]+3), col = cols, xaxt="n", xlab="", ylab = "", cex.axis=1.8)
title(ylab="Evolved \ndispersal\n distance", line=3, cex.lab=1.8)
points(kd_hi, type="l", col="red")
points(kd_lo, type="l", col="green")

vert_cut=0.2
image(x = 1:25, y= c(brks_ldc[1:nbins],brks_ldc[nbins]+0.001)[1:(nbins*vert_cut)], log(ldc_scan[,1:(nbins*vert_cut)]+3), col = cols, xaxt="n", xlab="", ylab = "", cex.axis=1.8)
title(ylab="Evolved \nactual\n dispersal", line=3, cex.lab=1.8)
points(avg_ldc, type="l", lwd=1, col="red")
points(ldc_hi, type="l", col="red")
points(ldc_lo, type="l", col="green")

# vert_cut=0.5
# image(x = 1:25, y= c(brks_rc[1:nbins],brks_rc[nbins]+0.001)[1:(nbins*vert_cut)], log(rc_scan[,1:(nbins*vert_cut)]+3), col = cols, xaxt="n", xlab="", ylab = "", cex.axis=1.8)
# title(ylab="Evolved \nresource\n consumed", xlab="Benefit of harvesting", line=3, cex.lab=1.8)
# 
# points(r_scan/1.0125e7*300, type="l", lwd=1, col="blue")
# points(avg_rc, type="l", lwd=1, col="red")

axis(side=1, at=as.integer(seq(1,length(bvec), length.out=5)), labels=sprintf("%4.1g",irvvec)[as.integer(seq(1,length(bvec), length.out=5))], cex.axis=1.8)
mtext("Slope of imitation response function", side = 1, line = 3, outer = F, cex=1.4)
