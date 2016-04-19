#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

int main(){
	int N = 1000
	int id[N];
	float x[N];
	for (int i=0; i<N; ++i) id[i] = i;
	
	x = sample(x = seq(1,70, length.out=N),size = N, replace=T, prob=c(seq(1,50, length.out=0.7*N),seq(70,1, length.out=0.3*N)) )
	hist(x)
	Ki.sd = 20
	prob = dnorm(x = x, mean=0, sd = Ki.sd)*sqrt(2*pi*Ki.sd^2)
	plot(prob~x)

	n.samples = 10000
	rs = numeric(n.samples)
	iters.req = numeric(n.samples)
	for (i in 1:n.samples){
	  accept= F
	  r = -1
	  niter = 0
	  while (accept == F){
		chosen.one = sample(id,1)
		if (runif(1) < prob[chosen.one]){
		  r = chosen.one
		  accept = T
		}
		niter = niter+1
	  }
	  rs[i] = r
	  iters.req[i] = niter
	  if (i %% 100 == 0) cat(i,"\n")
	}


}
