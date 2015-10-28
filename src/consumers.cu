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
	dL = L/nx;
	
	ke_lmax = I.getScalar("Ke_cutoff");		// bound for exploitation kernel in length units
	ke_nmax = int(ke_lmax/(L/nx));
	ke_sd = I.getScalar("Ke_sd");
	cout << "ke_nmax = " << ke_nmax << endl;

	getLastCudaError("resgrid");
	
	// allocate space
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

	consumers.resize(nc);
	for (int i=0; i<nc; ++i){
		consumers[i].pos   = runif2(0,L,0,L);
		consumers[i].pos_i = pos2cell(consumers[i].pos, dL);
		consumers[i].RT    = I.getScalar("RT0");
		consumers[i].Kdsd  = I.getScalar("kdsd0");
		consumers[i].h     = I.getScalar("h0");
		if (i<nc/2) consumers[i].h = 0.1;
		else 		consumers[i].h = 0.3;
		
		cout << consumers[i].pos.x << " " << consumers[i].pos.y << ", "   
			 << consumers[i].pos_i.x << " " << consumers[i].pos_i.y   << endl;

	}
	cudaMalloc((void**)&h_dev, sizeof(float)*nc);
	cudaMalloc((void**)&pos_i_dev, sizeof(int2)*nc);	
	cudaMalloc((void**)&rc_dev, sizeof(float)*nc);
	cudaMalloc((void**)&RT_dev, sizeof(float)*nc);
	cudaMalloc((void**)&kdsd_dev, sizeof(float)*nc);
	
	cudaMalloc((void**)&nd_dev, sizeof(float)*nc);
	cudaMalloc((void**)&lenDisp_dev, sizeof(float)*nc);

	vc_Tw = I.getScalar("payoff_Tw");
	cudaMalloc((void**)&vc_window_dev, sizeof(float)*nc*vc_Tw);
	cudaMalloc((void**)&vc_dev, sizeof(float)*nc);
	float vc
	
	cudaMemcpy2D((void*)h_dev, sizeof(float), (void*)&consumers[0].h, sizeof(Consumer), sizeof(float),  nc, cudaMemcpyHostToDevice);	
	cudaMemcpy2D((void*)RT_dev, sizeof(float), (void*)&consumers[0].RT, sizeof(Consumer), sizeof(float),  nc, cudaMemcpyHostToDevice);	
	cudaMemcpy2D((void*)kdsd_dev, sizeof(float), (void*)&consumers[0].Kdsd, sizeof(Consumer), sizeof(float),  nc, cudaMemcpyHostToDevice);	
	cudaMemcpy2D((void*)pos_i_dev, sizeof(int2), (void*)&consumers[0].pos_i, sizeof(Consumer), sizeof(int2),  nc, cudaMemcpyHostToDevice);	
	getLastCudaError("memcpy2D");

	
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
	cudaMemcpy2D((void*)&consumers[0].pos_i, sizeof(Consumer), (void*)pos_i_dev, sizeof(int2), sizeof(int2),  nc, cudaMemcpyDeviceToHost);	
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


__global__ void calc_exploitation_kernels_kernel(float* ke_all, int2* pos_cell, float* h, int nc, float* ke, int rkn, int nx){
	
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= nc) return;

	int ixc = pos_cell[tid].x;
	int iyc = pos_cell[tid].y;
	float hc = h[tid];

	for (int i=-rkn; i<=rkn; ++i){
		for (int j=-rkn; j<=rkn; ++j){
			int iK = makePeriodicID(ixc+i, nx);
			int jK = makePeriodicID(iyc+j, nx);
			ke_all[ix2(iK, jK, nx)] += hc * ke[ix2(i+rkn, j+rkn, (2*rkn+1))];
		}
	}	
	
}


void::ConsumerSystem::updateExploitationKernels(){
	
	for (int i=0; i<nx*ny; ++i) ke_all[i] = 0;	// reset exploitation kernels on host
	cudaMemcpy(ke_all_dev, ke_all, nx*ny*sizeof(float), cudaMemcpyHostToDevice); // reset ke_all_dev to zeros
	int nt = min(256, nc); int nb = 1; 
	calc_exploitation_kernels_kernel <<< nb, nt >>> (ke_all_dev, pos_i_dev, h_dev, nc, ke_dev, ke_nmax, nx);
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



__global__ void calc_resource_consumed_kernel(float *res, float* rc_vec, int2* pos_cell, float* h, int nc, float* ke, int rkn, int nx){

	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= nc) return;

	int ixc = pos_cell[tid].x;
	int iyc = pos_cell[tid].y;
	float hc = h[tid];

	float R_avail = 0;
	for (int i=-rkn; i<=rkn; ++i){
		for (int j=-rkn; j<=rkn; ++j){
			int iK = makePeriodicID(ixc+i, nx);
			int jK = makePeriodicID(iyc+j, nx);
			R_avail += res[ix2(iK, jK, nx)] * ke[ix2(i+rkn, j+rkn, (2*rkn+1))];
		}
	}	
	rc_vec[tid] = hc * R_avail;

}

void::ConsumerSystem::calcResConsumed(float * resource_dev){
	
	int nt = min(256, nc); int nb = 1; 
	calc_resource_consumed_kernel <<< nb, nt >>> (resource_dev, rc_dev, pos_i_dev, h_dev, nc, ke_dev, ke_nmax, nx);
	getLastCudaError("rc kernel");

//	cudaMemcpy2D(&consumers[0].rc, sizeof(Consumer), rc_dev, sizeof(float), sizeof(float), nc, cudaMemcpyDeviceToHost);
//	for (int i=0; i<nc; ++i){
//		cout << consumers[i].rc << "\n";
//	}
//	cout << "\n";
}


__global__ void disperse_kernel(float * res, int2 * pos_cell, 
								float * kdsd_vec, float * RT_vec, 
								curandState * RNG_states, 
								float L, int nc, int nx, 
								float * lenDisp, float * nd){
	
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= nc) return;
	
	int ixc = pos_cell[tid].x;
	int iyc = pos_cell[tid].y;

	float p_disperse = 1/(1+exp(10*(res[ix2(ixc, iyc, nx)] - RT_vec[tid])));	
	float b_disperse = curand_uniform(&RNG_states[tid]) < p_disperse;
	
	float len   = fabs(curand_normal(&RNG_states[tid]))*kdsd_vec[tid];
	float theta = curand_uniform(&RNG_states[tid])*2*3.14159;
	
	float xdisp = b_disperse*len*cos(theta);
	float ydisp = b_disperse*len*sin(theta);
	
	float2 xnew = make_float2(pos_cell[tid].x + xdisp, 
							  pos_cell[tid].y + ydisp );
	makePeriodic(xnew.x, 0, L);
	makePeriodic(xnew.y, 0, L);
	
	pos_cell[tid] = pos2cell(xnew, L/nx);
	lenDisp[tid]  = b_disperse*len;
	nd[tid] = b_disperse;
	
}


void ConsumerSystem::disperse(float * resource){
	int nt = min(256, nc); int nb = 1; 
	disperse_kernel <<< nb, nt >>> (resource, pos_i_dev, 
									kdsd_dev, RT_dev, 
									cs_dev_XWstates, 
									L, nc, nx,
									lenDisp_dev, nd_dev);	
}


// note: vc_window is as follows (because all quantities are in row arrays):
//		c1 c2 c3 c4 ....
//	t1   
//	t2
//	t3
//	...
//	tw	
//

__global__ void calc_payoffs_kernel(float * rc, float * lend, float * hc, int nc, 
									float * vc_window, float * vc, int tw, int t,
									float b, float cdisp, float charv ){

	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= nc) return;

	float v = b*rc[tid] - cdisp*lend[tid] - charv*hc[tid]*hc[tid]; 
	vc_window[ix2(tid, t%tw, nc)] = v;

	float vavg = 0;
	for (int i=0; i<tw; ++i) vavg = vc_window[ix2(tid, i, nc)];
	vavg = vavg/tw;
	
	vc[tid] = vavg;
									
} 

void ConsumerSystem::calcPayoffs(int t){
	int nt = min(256, nc); int nb = 1; 
	calc_payoffs_kernel <<<nb, nt >>> (rc_dev, lenDisp_dev, h_dev, nc, 
									   vc_window_dev, vc_dev, vc_Tw, t,
									   0.002, 0.1, 0.08);
									   
	cudaMemcpy2D(&consumers[0].rc, sizeof(Consumer), rc_dev, sizeof(float), sizeof(float), nc, cudaMemcpyDeviceToHost);
	cudaMemcpy2D(&consumers[0].h, sizeof(Consumer), h_dev, sizeof(float), sizeof(float), nc, cudaMemcpyDeviceToHost);
	cudaMemcpy2D(&consumers[0].ld, sizeof(Consumer), lenDisp_dev, sizeof(float), sizeof(float), nc, cudaMemcpyDeviceToHost);
	cudaMemcpy2D(&consumers[0].nd, sizeof(Consumer), nd_dev, sizeof(float), sizeof(float), nc, cudaMemcpyDeviceToHost);
	cudaMemcpy2D(&consumers[0].vc, sizeof(Consumer), vc_dev, sizeof(float), sizeof(float), nc, cudaMemcpyDeviceToHost);
	for (int i=0; i<nc; ++i){
		
	}
									   
}


