dat <- read.delim("/home/jaideep/gpu_codes/milker_killer_v2_graphics_sep/src/ts_h.txt", header=F)
dat <- dat[,-length(dat[1,])]

nt = nrow(dat)
times = dat[,1]
hc <- dat[,-1]

nbins = 100
hall = matrix(data = 0, nrow = nt, ncol=nbins)
hmean = numeric(nt)
for (i in 1:nt){
  hist <- hist(as.numeric(hc[i,]), breaks = seq(0,2,length.out = nbins+1), plot = F)
  hall[i,] = hist$counts
  hmean[i] = mean(as.numeric(dat[i,-1]))
}

image(x = times, y= seq(0,2,length.out = nbins+1), hall)
plot(hmean, type="l")

