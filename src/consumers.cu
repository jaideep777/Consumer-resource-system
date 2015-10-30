#include <cmath>
#include <fstream>
#include <cstdlib>
#include <string>
#include <curand_kernel.h>
using namespace std;

#include "../headers/consumers.h"
//#include "../headers/graphics.h"
#include "../utils/cuda_vector_math.cuh"
#include "../utils/cuda_device.h"


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// KERNEL to set up RANDOM GENERATOR STATES
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

__global__ void csRngStateSetup_kernel(int * rng_Seeds, curandState * rngStates, int nc){
	int tid_u = threadIdx.x + blockIdx.x*blockDim.x;	// each block produces different random numbers
	if (tid_u >= nc) return;
	
	curand_init (rng_Seeds[tid_u], 0, 0, &rngStates[tid_u]);
}

#define CS_PP_SEED 777 //time(NULL)

void ConsumerSystem::initRNG(){
	srand(CS_PP_SEED);
	for (int i=0; i<nc; ++i) cs_seeds_h[i] = rand(); 
	cudaMemcpy( cs_seeds_dev, cs_seeds_h, sizeof(int)*nc, cudaMemcpyHostToDevice);

	int nt = min(256, nc), nb = (nc-1)/nt+1;
	csRngStateSetup_kernel <<< nb, nt>>> (cs_seeds_dev, cs_dev_XWstates, nc);
	getLastCudaError("RNG_kernel_launch");
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//          CONSUMER SYSTEM
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void ConsumerSystem::init(Initializer &I){
	nx = I.getScalar("nx"); 
	ny = I.getScalar("ny");
	L  = I.getScalar("L");
	nc = I.getScalar("nc");
	dt = I.getScalar("dt");
	dL = L/nx;
	
	b = I.getScalar("b");
	cd = I.getScalar("cd");
	ch = I.getScalar("ch");
	b_imit_h  = I.getScalar("imitate_h");
	b_imit_rt = I.getScalar("imitate_RT");
	b_imit_kd = I.getScalar("imitate_Kd");
	rImit = I.getScalar("imitation_rate");
	
	ke_lmax = I.getScalar("Ke_cutoff");		// bound for exploitation kernel in length units
	ke_nmax = int(ke_lmax/dL);
	ke_sd = I.getScalar("Ke_sd");
	cout << "ke_nmax = " << ke_nmax << endl;

	
	// calculate and allocate exploitation kernels
	int ke_arrlen = (2*ke_nmax+1)*(2*ke_nmax+1);  // ke goes from [-ke_nmax, ke_nmax]
	ke = new float[ke_arrlen];
	for (int j=-ke_nmax; j<=ke_nmax; ++j){
		for (int i=-ke_nmax; i<=ke_nmax; ++i){
			float ker = 1-(i*i+j*j)*dL*dL/ke_sd/ke_sd;
			if (ker < 0) ker = 0;
			ke[ix2(i+ke_nmax, j+ke_nmax, 2*ke_nmax+1)] = ker;
		}
	}
	cudaMalloc((void**)&ke_dev, sizeof(float)*ke_arrlen);
	cudaMemcpy(ke_dev, ke, ke_arrlen*sizeof(float), cudaMemcpyHostToDevice);
	
	ke_all = new float[nx*ny];
	for (int i=0; i<nx*ny; ++i) ke_all[i] = 0;
	cudaMalloc((void**)&ke_all_dev, sizeof(float)*nx*ny);

	// init consumers
	consumers.resize(nc);
	for (int i=0; i<nc; ++i){
		consumers[i].pos   = runif2(0,L,0,L);
		consumers[i].pos_i = pos2cell(consumers[i].pos, dL);
		consumers[i].RT    = I.getScalar("RT0");
		consumers[i].Kdsd  = I.getScalar("kdsd0");
		consumers[i].h     = I.getScalar("h0");
		consumers[i].vc    = 0;
		consumers[i].vc_avg = 0;
		if (i<nc/2) consumers[i].h = 0.01;
		else 		consumers[i].h = 0.01;
		
//		cout << consumers[i].pos.x << " " << consumers[i].pos.y << ", "   
//			 << consumers[i].pos_i.x << " " << consumers[i].pos_i.y   << endl;

	}

	cudaMalloc((void**)&consumers_dev, nc*sizeof(Consumer));
	cudaMemcpy(consumers_dev, &consumers[0], nc*sizeof(Consumer), cudaMemcpyHostToDevice);

	vc_Tw = I.getScalar("payoff_Tw");
	vc_window = new float[nc*vc_Tw];
	for (int i=0; i<nc*vc_Tw; ++i) vc_window[i] = 0;
	cudaMalloc((void**)&vc_window_dev, sizeof(float)*nc*vc_Tw);
	cudaMemcpy(vc_window_dev, vc_window, nc*vc_Tw*sizeof(float), cudaMemcpyHostToDevice);

	cs_seeds_h = new int[nc];
	cudaMalloc((void**)&cs_seeds_dev, nc*sizeof(int));
	cudaMalloc((void**)&cs_dev_XWstates, nc*sizeof(curandState));
	
	initRNG();
	
	// create Pointset shape to display consumers
	cons_shape = PointSet("res", false, nc, 0, L);	// the res shader is a generic shader for colormaps
	cons_shape.nVertices = nc;
	cons_shape.createShaders();
	float2 tmp[nc]; 
	for (int i=0; i<nc; ++i) {
		tmp[i] = cell2pos(consumers[i].pos_i, dL);
//		cout << consumers[i].pos_i.x << " " << consumers[i].pos_i.y  << ", " << tmp[i].x << " " << tmp[i].y << endl;
	}
	cons_shape.createVBO(tmp, cons_shape.nVertices*sizeof(float2));	
	cons_shape.createColorBuffer();
	cons_shape.setDefaultColor();
	cons_shape.palette = createPalette_ramp(nc, Colour_rgb(0,0.9,0), Colour_rgb(1,0,0));
	printPalette(cons_shape.palette);
}


void ConsumerSystem::graphics_updateArrays(){

	// positions buffer
	cudaMemcpy2D((void*)&consumers[0].pos_i, sizeof(Consumer), (void*)&consumers_dev[0].pos_i, sizeof(Consumer), sizeof(int2),  nc, cudaMemcpyDeviceToHost);	
	float2 tmp[nc]; 
	for (int i=0; i<nc; ++i) {
		tmp[i] = cell2pos(consumers[i].pos_i, dL);
//		cout << consumers[i].pos_i.x << " " << consumers[i].pos_i.y  << ", " << tmp[i].x << " " << tmp[i].y << endl;
	}
	glBindBuffer(GL_ARRAY_BUFFER, cons_shape.vbo_ids[0]); 	// Bring 1st buffer into current openGL context
	glBufferData(GL_ARRAY_BUFFER, nc*sizeof(float2), tmp, GL_DYNAMIC_DRAW); 
	glBindBuffer(GL_ARRAY_BUFFER, 0); 	// Bring 1st buffer into current openGL context
	
	// color buffer
	float h_tmp[nc];
	for (int i=0; i<nc; ++i) {
		h_tmp[i] = consumers[i].h;
	}
	cons_shape.updateColors(h_tmp, nc);
	
}


// blocks run over consumers
// threads run over gridcells in Ke

__global__ void calc_exploitation_kernels_kernel(float* ke_all, Consumer* cons, int nc, float* ke, int rkn, int nx){
	
	int ke_nx = (2*rkn+1);

	if (threadIdx.x >= ke_nx*ke_nx) return;
	if (blockIdx.x >= nc) return;

	// get position and harvesting rate of consumer from block ID
	int ixc = cons[blockIdx.x].pos_i.x;
	int iyc = cons[blockIdx.x].pos_i.y;
	float hc = cons[blockIdx.x].h;

	// i and j are 2D indices on Ke. They range from [-rkn, rkn]
	int i = (threadIdx.x % ke_nx) - rkn;	
	int j = (threadIdx.x / ke_nx) - rkn;
	
	// i and j are also 2D indices on Ke_all. They must be periodic over the grd boundary
	int iK = makePeriodicID(ixc+i, nx);
	int jK = makePeriodicID(iyc+j, nx);
	
	// multiple consumers might be adding to the same grid cell. Hence need atomic add.
	atomicAdd(&ke_all[ix2(iK, jK, nx)], hc * ke[ix2(i+rkn, j+rkn, ke_nx)]);

}


__global__ void resetKeAll_kernel(float * keall, int nx, int ny){
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= nx*ny) return;

	keall[tid] = 0;
}



void::ConsumerSystem::updateExploitationKernels(){
	
//	for (int i=0; i<nx*ny; ++i) ke_all[i] = 0;	// reset exploitation kernels on host
//	cudaMemcpy(ke_all_dev, ke_all, nx*ny*sizeof(float), cudaMemcpyHostToDevice); // reset ke_all_dev to zeros
	int nt = 256; int nb = (nx*ny-1)/nt + 1;
	resetKeAll_kernel <<< nb, nt>>> (ke_all_dev, nx, ny);

	nt = (2*ke_nmax+1)*(2*ke_nmax+1); nb = nc; 
	calc_exploitation_kernels_kernel <<< nb, nt >>> (ke_all_dev, consumers_dev, nc, ke_dev, ke_nmax, nx);
	getLastCudaError("exploitation kernel");

//	cudaMemcpy(ke_all, ke_all_dev, nx*ny*sizeof(float), cudaMemcpyDeviceToHost);
//	ofstream fout("ke_all.txt");
//	for (int j=0; j<ny; ++j){
//		for (int i=0; i<nx; ++i){
//			fout << ke_all[ix2(i,j,nx)] << "\t";
//		}
//		fout << "\n";
//	}
//	fout << "\n";
}


__global__ void resetRc_kernel(Consumer * cons, int nc){
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= nc) return;

	cons[tid].rc = 0;
}


// in this kernel, blocks run over grid cells and threads run over consumers.
// this is done to minimize conflicts in atomicAdd
// this puts limitation on max number of consumers, but thats OK.
__global__ void calc_resource_consumed_kernel(float *res, Consumer* cons, int nc, float* ke, int rkn, int nx, float dt){

	int ke_nx = (2*rkn+1);

	if (blockIdx.x >= ke_nx*ke_nx) return;
	if (threadIdx.x >= nc) return;

	// get position and harvesting rate of consumer from block ID
	int ixc = cons[threadIdx.x].pos_i.x;
	int iyc = cons[threadIdx.x].pos_i.y;
	float hc = cons[threadIdx.x].h;

	// i and j are 2D indices on Ke. They range from [-rkn, rkn]
	int i = (blockIdx.x % ke_nx) - rkn;	
	int j = (blockIdx.x / ke_nx) - rkn;
	
	// i and j are also 2D indices on Ke_all. They must be periodic over the grd boundary
	int iK = makePeriodicID(ixc+i, nx);
	int jK = makePeriodicID(iyc+j, nx);

	atomicAdd(&cons[threadIdx.x].rc, hc * res[ix2(iK, jK, nx)] * ke[ix2(i+rkn, j+rkn, (2*rkn+1))] * dt);

}

void::ConsumerSystem::calcResConsumed(float * resource_dev){
	
	int nt = min(256, nc); int nb = (nc-1)/nt+1; 
	resetRc_kernel <<<nb, nt >>> (consumers_dev, nc);
	
	nb = (2*ke_nmax+1)*(2*ke_nmax+1); nt = nc; 
	calc_resource_consumed_kernel <<< nb, nt >>> (resource_dev, consumers_dev, nc, ke_dev, ke_nmax, nx, dt);
	getLastCudaError("rc kernel");


//	cudaMemcpy2D(&consumers[0].rc, sizeof(Consumer), &consumers_dev[0].rc, sizeof(Consumer), sizeof(float), nc, cudaMemcpyDeviceToHost);
//	for (int i=0; i<nc; ++i){
//		cout << consumers[i].rc << "\n";
//	}
//	cout << "\n\n";
}


__global__ void disperse_kernel(float * res, Consumer* cons, 
								curandState * RNG_states, 
								float L, int nc, int nx ){
	
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= nc) return;
	
	int ixc = cons[tid].pos_i.x;
	int iyc = cons[tid].pos_i.y;

	float p_disperse = 1/(1+exp(10*(res[ix2(ixc, iyc, nx)] - cons[tid].RT)));	
	float b_disperse = curand_uniform(&RNG_states[tid]) < p_disperse;
	
	float len   = fabs(curand_normal(&RNG_states[tid]))*cons[tid].Kdsd;
	float theta = curand_uniform(&RNG_states[tid])*2*3.14159;
	
	float xdisp = b_disperse*len*cos(theta);
	float ydisp = b_disperse*len*sin(theta);
	
	float2 xnew = make_float2(ixc + xdisp, iyc + ydisp );
	makePeriodic(xnew.x, 0, L);
	makePeriodic(xnew.y, 0, L);
	
	cons[tid].pos_i = pos2cell(xnew, L/nx);
	cons[tid].ld  = b_disperse*len;
	cons[tid].nd = b_disperse;
	
}


void ConsumerSystem::disperse(float * resource){
	int nt = min(256, nc); int nb = (nt-1)/nc+1; 
	disperse_kernel <<< nb, nt >>> (resource, consumers_dev, 
									cs_dev_XWstates, 
									L, nc, nx);	
}


// note: vc_window is as follows (because all quantities are in row arrays):
//		t1 t2 t3  .... tw
//	c1   
//	c2
//	c3
//	...
//	cn	
//

__global__ void calc_payoffs_kernel(Consumer * cons, int nc, float b, float cd, float ch){
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= nc) return;
	
	cons[tid].vc = b*cons[tid].rc - cd*cons[tid].ld - ch*cons[tid].h*cons[tid].h;
}

__global__ void avg_payoffs_kernel(Consumer * cons, int nc, int tw){
	
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= nc) return;
	
	cons[tid].vc_avg += (cons[tid].vc - cons[tid].vc_x)/tw;
}


void ConsumerSystem::calcPayoff(int t){
	
	// calculate payoffs for current time step
	int nt = min(256, nc); int nb = (nt-1)/nc+1;
	calc_payoffs_kernel <<<nb, nt>>> (consumers_dev, nc, b, cd, ch); 
	getLastCudaError("payoffs kernel");

	// update the vc_window (discard the oldest column, and replace with current)
	int win_id = t % vc_Tw;
	cudaMemcpy2D(&consumers_dev[0].vc_x, sizeof(Consumer), &vc_window_dev[win_id], sizeof(float)*vc_Tw, sizeof(float), nc, cudaMemcpyDeviceToDevice); 
	cudaMemcpy2D(&vc_window_dev[win_id], sizeof(float)*vc_Tw, &consumers_dev[0].vc, sizeof(Consumer), sizeof(float), nc, cudaMemcpyDeviceToDevice); 

	nt = min(256, nc); nb = (nt-1)/nc+1;
	avg_payoffs_kernel <<< nb, nt >>>  (consumers_dev, nc, vc_Tw);
	getLastCudaError("avg payoffs kernel");

//	cudaMemcpy(vc_window, vc_window_dev, sizeof(float)*nc*vc_Tw, cudaMemcpyDeviceToHost);
//	cudaMemcpy(&consumers[0], consumers_dev, sizeof(Consumer)*nc, cudaMemcpyDeviceToHost);	
//	cout << "\n";
//	for (int i=0; i<nc; ++i){
//		cout << consumers[i].vc << ", " << consumers[i].vc_x << ":\t";
//		for (int t=0; t<vc_Tw; ++t){
//			cout << vc_window[i*vc_Tw + t] << "\t";		
//		}
//		cout << "| " << consumers[i].vc_avg << "\n";
//	}
//	cout << "\n";

}


__global__ void imitate_global_kernel(Consumer* cons, curandState * RNG_states, int nc, float rImit, float dt,
									  bool b_ih, bool b_irt, bool b_ikd){
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= nc) return;
	
	bool b_imit = curand_uniform(&RNG_states[tid]) < rImit*dt;
	if (b_imit){
	
		int id_whom = (1-curand_uniform(&RNG_states[tid]))*(nc-1);	// curand_uniform gives (0,1]. do 1-X to get [0,1)
		float self = id_whom == tid;
		id_whom = self*(nc-1) + (1-self)*id_whom;

		float dv = cons[id_whom].vc_avg - cons[tid].vc_avg;
		float imitation_prob = float(dv > 0);

		if (curand_uniform(&RNG_states[tid]) <= imitation_prob) { 

			float h_new    = cons[id_whom].h    + 0.02*curand_normal(&RNG_states[tid]);	
			float RT_new   = cons[id_whom].RT   + 1.00*curand_normal(&RNG_states[tid]);	
			float Kdsd_new = cons[id_whom].Kdsd + 0.20*curand_normal(&RNG_states[tid]);	

			if (b_ih)  cons[tid].h    = clamp(h_new, 0.f, h_new);	
			if (b_irt) cons[tid].RT   = clamp(RT_new, 0.f, RT_new);	
			if (b_ikd) cons[tid].Kdsd = clamp(Kdsd_new, 0.f, Kdsd_new);	
		}
	}

}


void ConsumerSystem::imitate_global(){
	
	int nt = min(256, nc); int nb = (nt-1)/nc+1;
	imitate_global_kernel <<< nb, nt >>> (consumers_dev, cs_dev_XWstates, nc, rImit, dt,
										  b_imit_h, b_imit_rt, b_imit_kd);
	getLastCudaError("imitate global kernel");
		
}




//void ConsumerSystem::calcPayoffs(int t){
//	int nb = 1; 
//	calc_payoffs_kernel <<<nb, nt >>> (rc_dev, lenDisp_dev, h_dev, nc, 
//									   vc_window_dev, vc_dev, vc_Tw, t,
//									   0.002, 0.1, 0.08);
//									   
//	cudaMemcpy2D(&consumers[0].rc, sizeof(Consumer), rc_dev, sizeof(float), sizeof(float), nc, cudaMemcpyDeviceToHost);
//	cudaMemcpy2D(&consumers[0].h, sizeof(Consumer), h_dev, sizeof(float), sizeof(float), nc, cudaMemcpyDeviceToHost);
//	cudaMemcpy2D(&consumers[0].ld, sizeof(Consumer), lenDisp_dev, sizeof(float), sizeof(float), nc, cudaMemcpyDeviceToHost);
//	cudaMemcpy2D(&consumers[0].nd, sizeof(Consumer), nd_dev, sizeof(float), sizeof(float), nc, cudaMemcpyDeviceToHost);
//	cudaMemcpy2D(&consumers[0].vc, sizeof(Consumer), vc_dev, sizeof(float), sizeof(float), nc, cudaMemcpyDeviceToHost);
//	for (int i=0; i<nc; ++i){
//		
//	}
//									   
//}


