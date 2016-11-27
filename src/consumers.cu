#include <cmath>
#include <fstream>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <curand_kernel.h>

using namespace std;

#include "../headers/consumers.h"
#include "../headers/resource.h"
//#include "../headers/graphics.h"
#include "../utils/cuda_vector_math.cuh"
#include "../utils/cuda_device.h"
#include "../utils/simple_histogram.h"


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
	
	mu_h = I.getScalar("mu_h");
	mu_RT = I.getScalar("mu_RT");
	mu_kd = I.getScalar("mu_kd");
	
	graphics = I.getScalar("graphicsQual")>0;
	
	ke_lmax = I.getScalar("Ke_cutoff");		// bound for exploitation kernel in length units
	ke_nmax = int(ke_lmax/dL);
	ke_sd = I.getScalar("Ke_sd");
	cout << "ke_nmax = " << ke_nmax << endl;

	ki_sd = I.getScalar("Ki_sd");
	
	// calculate and allocate exploitation kernels
	int ke_arrlen = (2*ke_nmax+1)*(2*ke_nmax+1);  // ke goes from [-ke_nmax, ke_nmax]
	ke = new float[ke_arrlen];
	float kernelSum = 0;
	for (int j=-ke_nmax; j<=ke_nmax; ++j){
		for (int i=-ke_nmax; i<=ke_nmax; ++i){
			float ker = 1-(i*i+j*j)*dL*dL/ke_sd/ke_sd;
			if (ker < 0) ker = 0;
			ke[ix2(i+ke_nmax, j+ke_nmax, 2*ke_nmax+1)] = ker;
			kernelSum += ker;
		}
	}
	cout << "Kernel Sum = " << kernelSum << endl;
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

		consumers[i].ld_cumm = 0;
		consumers[i].nd_cumm = 0;
		
//		if (i<nc/2) consumers[i].h = 0.1;
//		else 		consumers[i].h = 0.3;
		
//		cout << consumers[i].pos.x << " " << consumers[i].pos.y << ", "   
//			 << consumers[i].pos_i.x << " " << consumers[i].pos_i.y   << endl;

	}

	cudaMalloc((void**)&consumers_dev, nc*sizeof(Consumer));
	cudaMalloc((void**)&consumers_child_dev, nc*sizeof(Consumer));
	cudaMemcpy(consumers_dev, &consumers[0], nc*sizeof(Consumer), cudaMemcpyHostToDevice);
	cudaMemcpy(consumers_child_dev, &consumers[0], nc*sizeof(Consumer), cudaMemcpyHostToDevice);

	// payoff window
	vc_Tw = I.getScalar("payoff_Tw");
	vc_window = new float[nc*vc_Tw];
	for (int i=0; i<nc*vc_Tw; ++i) vc_window[i] = 0;
	cudaMalloc((void**)&vc_window_dev, sizeof(float)*nc*vc_Tw);
	cudaMemcpy(vc_window_dev, vc_window, nc*vc_Tw*sizeof(float), cudaMemcpyHostToDevice);

	// output interval
	out_Tw = I.getScalar("out_Tw");

	// RNG
	cs_seeds_h = new int[nc];
	cudaMalloc((void**)&cs_seeds_dev, nc*sizeof(int));
	cudaMalloc((void**)&cs_dev_XWstates, nc*sizeof(curandState));
	
	initRNG();
	
	// imitation stuff
	cudaMalloc((void**)&pd_dev, nc*nc*sizeof(float));
	cudaMalloc((void**)&cumm_prob_dev, nc*(nc+1)*sizeof(float));
	//cudaMalloc((void**)&wsum_dev, nc*sizeof(float));
	
	// array of zeros
	nc_zeros = new float[nc];
	for (int i=0; i<nc; ++i) nc_zeros[i] = 0; // !!
	
	// init histograms
	int nbins = 100;
	brks_h.resize(nbins); 
	for (int i=0; i<nbins; ++i) brks_h[i] = 0 + double(i)/(nbins-1)*2.5;
	brks_h.push_back(1000);	

	brks_rc.resize(nbins); 
	for (int i=0; i<nbins; ++i) brks_rc[i] = 0 + double(i)/(nbins-1)*600.0;
	brks_rc.push_back(1000000);	

	brks_kd.resize(nbins); 
	for (int i=0; i<nbins; ++i) brks_kd[i] = 0 + double(i)/(nbins-1)*25.0;
	brks_kd.push_back(1000000);	

	brks_ldc.resize(nbins); 
	for (int i=0; i<nbins; ++i) brks_ldc[i] = 0 + double(i)/(nbins-1)*50.0;
	brks_ldc.push_back(1000000);	

	brks_ndc.resize(nbins); 
	for (int i=0; i<nbins; ++i) brks_ndc[i] = 0 + double(i)/(nbins-1)*1.0;
	brks_ndc.push_back(1000000);	
	
	if (graphics){
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

	//printConsumers();

}


void ConsumerSystem::initIO(Initializer &I){
       
	// ~~~~~~~~~~~~~~~~~ EXPT DESC ~~~~~~~~~~~~~~~~~~~~      
	stringstream sout;
	sout << setprecision(3);
	sout << I.getString("exptName");
	sout << "_T("   << I.getScalar("nsteps")/1000
		 << ")_N("  << nc
		 << ")_RT(" << ((b_imit_rt)? -1.0:I.getScalar("RT0"))
		 << ")_kd(" << ((b_imit_kd)? -1.0:I.getScalar("kdsd0"))
		 << ")_h("  << ((b_imit_h)?  -1.0:I.getScalar("h0"))
		 << ")_rI(" << rImit
		 << ")_kI(" << ki_sd
		 << ")_L("  << L
		 << ")_nx(" << nx
		 << ")_b(" << b
		 << ")_cd(" << cd
		 << ")_ch(" << ch;
	if (I.getString("exptName") == "het") sout << ")_tmu(" << tmu;
	sout << ")_muh(" << mu_h
		 << ")_murt(" << mu_RT
		 << ")_mukd(" << mu_kd
		 << ")_twv(" << vc_Tw
		 << ")";

	string exptDesc = sout.str(); sout.clear();

	cout << "experiment: " << exptDesc << endl; 

	if (I.getScalar("dataOut")>0) {
		string output_dir = I.getString("homeDir") + "/" + I.getString("outDir");
		system(string("mkdir " + output_dir).c_str());

		fout_h[0].open (string(output_dir + "/h_" + exptDesc).c_str());
		fout_rc[0].open(string(output_dir + "/rc_" + exptDesc).c_str());
		fout_sd[0].open(string(output_dir + "/kd_" + exptDesc).c_str());
		fout_ldc[0].open(string(output_dir + "/ldc_" + exptDesc).c_str());
		fout_ndc[0].open(string(output_dir + "/ndc_" + exptDesc).c_str());

//		if (!fout_h[0].is_open()) cout << "failed to open h file." << endl;
//		if (!fout_rc[0].is_open()) cout << "failed to open rc file." << endl;
//		if (!fout_sd[0].is_open()) cout << "failed to open sd file." << endl;

		fout_h[1].open (string(output_dir + "/hist_h_" + exptDesc).c_str());    
		fout_rc[1].open(string(output_dir + "/hist_rc_" + exptDesc).c_str());   
		fout_sd[1].open(string(output_dir + "/hist_kd_" + exptDesc).c_str());   
		fout_ldc[1].open(string(output_dir + "/hist_ldc_" + exptDesc).c_str());
		fout_ndc[1].open(string(output_dir + "/hist_ndc_" + exptDesc).c_str());

		fout_r.open(string(output_dir + "/r_total_" + exptDesc).c_str());
//		fout_wsum.open(string(output_dir + "/wsum_" + exptDesc).c_str());
	}

}


void ConsumerSystem::closeIO(){
	fout_h[0].close();      
	fout_rc[0].close();     
	fout_sd[0].close();     
	fout_ldc[0].close();     
	fout_ndc[0].close();     

	fout_h[1].close();      
	fout_rc[1].close();     
	fout_sd[1].close();     
	fout_ldc[1].close();     
	fout_ndc[1].close();     

	fout_r.close();
//	fout_wsum.close();
}


void ConsumerSystem::writeState(int istep, ResourceGrid * resgrid){

	cudaMemcpy(&consumers[0], consumers_dev, sizeof(Consumer)*nc, cudaMemcpyDeviceToHost);

	// output indivuduals
	fout_h[0] << istep << "\t";
	fout_rc[0] << istep << "\t";
	fout_sd[0] << istep << "\t";
	fout_ldc[0] << istep << "\t";
	fout_ndc[0] << istep << "\t";
//	fout_wsum << istep << "\t";
	for (int i=0; i<nc; ++i){
		fout_h[0]  << consumers[i].h << "\t";
		fout_rc[0] << consumers[i].rc << "\t";
		fout_sd[0] << consumers[i].Kdsd << "\t";
		fout_ldc[0] << consumers[i].ld_cumm << "\t";
		fout_ndc[0] << consumers[i].nd_cumm << "\t";
//		fout_wsum << consumers[i].wsum << "\t";
	}
	fout_h[0] << endl;
	fout_rc[0] << endl;
	fout_sd[0] << endl;
	fout_ldc[0] << endl;
	fout_ndc[0] << endl;
//	fout_wsum << endl;

	// output histograms
	vector <float> hvec(nc), rcvec(nc), kdvec(nc), ldcvec(nc), ndcvec(nc);
	for (int i=0; i<nc; ++i){
		hvec[i]  = consumers[i].h;
		rcvec[i] = consumers[i].rc;
		kdvec[i] = consumers[i].Kdsd;
		ldcvec[i] = consumers[i].ld_cumm/out_Tw;
		ndcvec[i] = consumers[i].nd_cumm/out_Tw;
	}
	Histogram h_h(hvec,  brks_h);
	Histogram h_rc(rcvec, brks_rc);
	Histogram h_kd(kdvec, brks_kd);
	Histogram h_ldc(ldcvec, brks_ldc);
	Histogram h_ndc(ndcvec, brks_ndc);

	vector <float> counts_h  = h_h.getCounts();
	vector <float> counts_rc = h_rc.getCounts();
	vector <float> counts_kd = h_kd.getCounts();
	vector <float> counts_ldc = h_ldc.getCounts();
	vector <float> counts_ndc = h_ndc.getCounts();
	
	for (int i=0; i<counts_h.size(); ++i){
		fout_h[1]  << counts_h[i] << "\t";
		fout_rc[1] << counts_rc[i] << "\t";
		fout_sd[1] << counts_kd[i] << "\t";
		fout_ldc[1] << counts_ldc[i] << "\t";
		fout_ndc[1] << counts_ndc[i] << "\t";
	}
	fout_h[1] << endl;
	fout_rc[1] << endl;
	fout_sd[1] << endl;
	fout_ldc[1] << endl;
	fout_ndc[1] << endl;
	
	// output total resource remaining
	float totRes = resgrid->sumResource();
	fout_r << totRes << endl;

	// reset cummulative values ld and nd
	cudaMemcpy2D(&consumers_dev[0].ld_cumm, sizeof(Consumer), nc_zeros, sizeof(float), sizeof(float), nc, cudaMemcpyHostToDevice); 
	cudaMemcpy2D(&consumers_dev[0].nd_cumm, sizeof(Consumer), nc_zeros, sizeof(float), sizeof(float), nc, cudaMemcpyHostToDevice); 
	

/*	// DEBUG: Output Ke_all and res at step 200
	if (istep == 20000){
		cudaMemcpy(resgrid->res, resgrid->res_dev, nx*ny*sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(ke_all, ke_all_dev, nx*ny*sizeof(float), cudaMemcpyDeviceToHost);

		ofstream fout;
		
		fout.open("ke200.txt");
		for (int i=0; i<nx; ++i){
			for (int j=0; j<ny; ++j){
				fout << ke_all[i*nx + j] << "\t";
			}
			fout << endl;
		}
		fout << endl;
		fout.close();
		
		fout.open("res200.txt");
		for (int i=0; i<nx; ++i){
			for (int j=0; j<ny; ++j){
				fout << resgrid->res[i*nx + j] << "\t";
			}
			fout << endl;
		}
		fout << endl;

	}
*/

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
	getLastCudaError("reset ke_all kernel");

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
	getLastCudaError("reset_rc kernel");
	
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
								float L, float dL, int nc, int nx ){
	
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
	
	float2 xnew = make_float2(ixc*dL+dL/2 + xdisp, iyc*dL+dL/2 + ydisp);
	makePeriodic(xnew.x, 0, L);
	makePeriodic(xnew.y, 0, L);
	
	cons[tid].pos_i = pos2cell(xnew, dL);
	cons[tid].ld  = b_disperse*len;
	cons[tid].nd = b_disperse;
	
	cons[tid].ld_cumm += b_disperse*len;
	cons[tid].nd_cumm += b_disperse;
}


void ConsumerSystem::disperse(float * resource){
	int nt = min(256, nc); int nb = (nc-1)/nt+1; 
	disperse_kernel <<< nb, nt >>> (resource, consumers_dev, 
									cs_dev_XWstates, 
									L, dL, nc, nx);	
	getLastCudaError("disperse kernel");								
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
	
	cons[tid].vc = b*cons[tid].rc - cd*cons[tid].ld - ch*cons[tid].h*cons[tid].h - 0.001*cons[tid].Kdsd;
}

__global__ void avg_payoffs_kernel(Consumer * cons, int nc, int tw){
	
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= nc) return;
	
	cons[tid].vc_avg += (cons[tid].vc - cons[tid].vc_x)/tw;
}


void ConsumerSystem::calcPayoff(int t){
	
	// calculate payoffs for current time step
	int nt = min(256, nc); int nb = (nc-1)/nt+1;
	calc_payoffs_kernel <<<nb, nt>>> (consumers_dev, nc, b, cd, ch); 
	getLastCudaError("payoffs kernel");

	// update the vc_window (discard the oldest column, and replace with current)
	int win_id = t % vc_Tw;
	cudaMemcpy2D(&consumers_dev[0].vc_x, sizeof(Consumer), &vc_window_dev[win_id], sizeof(float)*vc_Tw, sizeof(float), nc, cudaMemcpyDeviceToDevice); 
	cudaMemcpy2D(&vc_window_dev[win_id], sizeof(float)*vc_Tw, &consumers_dev[0].vc, sizeof(Consumer), sizeof(float), nc, cudaMemcpyDeviceToDevice); 

	nt = min(256, nc); nb = (nc-1)/nt+1;
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

//		cout << consumers[i].rc << " " << consumers[i].ld << "\n";
//	}
//	cout << "\n";

}


/*------------------------------------------------------------------------------

	periodic pairwise distances on GPU
	returns f(x) where x is pairwise distance, and f is gaussian with sd ki

------------------------------------------------------------------------------*/

__global__ void pairwise_distance_kernel(Consumer * cons, float * pd, int nc, float ki, float L, float dL){
	int tx = threadIdx.x, ty = threadIdx.y;
	int bx = blockIdx.x, by = blockIdx.y;
	
	int o = by*16*nc + ty*nc + bx*16 + tx;
	int i = 16*by + ty;
	int j = 16*bx + tx;
//	float x = length(periodicDisplacement(cons[i].pos_i*dL+dL/2, cons[j].pos_i*dL+dL/2, L, L));

	float2 v2other = periodicDisplacement(	make_float2(cons[i].pos_i.x*dL+dL/2, cons[i].pos_i.y*dL+dL/2), 
											make_float2(cons[j].pos_i.x*dL+dL/2, cons[j].pos_i.y*dL+dL/2), 
											L, L );
	float x = length(v2other);

	float prob = expf(-x*x/2/ki/ki);
	pd[o] = prob;
	
	atomicAdd(&cons[by*16+ty].wsum, prob); 
	
}


/*------------------------------------------------------------------------------

	Roulette sampling on GPU

------------------------------------------------------------------------------*/
__global__ void sample_roulette_kernel(int nc, float* prob, float* cumm_prob_all, Consumer * cons, curandState* p_rngStates){
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= nc) return;
	
	float * cumm_prob = &cumm_prob_all[tid*(nc+1)];	// cumm_prob_all row #tid

	// create ranges array 
	cumm_prob[0] = 0;
	for (int i=0; i<nc; ++i) cumm_prob[i+1] = cumm_prob[i] + prob[tid*nc + i]*float(i!=tid);

	// generate range selector
	float a = 1-curand_uniform(&p_rngStates[tid]);  // a is range selector. Must be in [0,1)
	a *= cumm_prob[nc];				// transform a into [0, sum(weights) )

	// get least upper bound
	int r = bin_search_lub(a, cumm_prob, nc+1); 

	// we want r-1 because upper bound is the right edge of the range for the desired element
	cons[tid].imit_whom = r-1;
	
}


/*------------------------------------------------------------------------------

	Rejection sampling on GPU

------------------------------------------------------------------------------*/
__global__ void sample_reject_kernel(int nc, float* pd, Consumer * cons, int *iters, curandState* p_rngStates){
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= nc) return;
	
	bool accept = false;
	int r = -1;
	int niter=0;
	while(!accept){
		// choose a random individual other than self
		int chosen_one = (1-curand_uniform(&p_rngStates[tid]))*(nc-1);
		int self = (chosen_one == tid);
		chosen_one = self*(nc-1) + (1-self)*chosen_one;

		// get probability of choosing for imitation from imitation kernel
		float prob = pd[tid*nc + chosen_one];
		
		// if rnd is < prob, accept.
		if (curand_uniform(&p_rngStates[tid]) < prob){
			r = chosen_one;
			accept=true;
		}
		
		// if iters exceed limit, just choose self.. i.e. no imitation
		++niter;
		if (niter > 5000){
			r = tid;
			break;
		}
	}
	cons[tid].imit_whom = r;	// the id of the individual that consumer tid chooses to imitate
//	iters[tid] = niter; // number of iterations required for choosing
			
}

__global__ void imitate_sync_kernel(Consumer* cons, Consumer* cons_child, curandState * RNG_states, int nc, float rImit, float dt,
									  bool b_ih, bool b_irt, bool b_ikd,
									  float mu_h, float mu_RT, float mu_kd){
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= nc) return;
	
//	float wsum = cons[tid].wsum;
	float imit_feasibility = 1; //wsum*wsum/ ((float(nc)/8)*(float(nc)/8) + wsum*wsum);

	bool b_imit = curand_uniform(&RNG_states[tid]) < rImit*dt;	
	bool b_imit2 = curand_uniform(&RNG_states[tid]) < imit_feasibility;	
	if (b_imit && b_imit2){
	
		int id_whom = cons[tid].imit_whom;
		
		float dv = cons[id_whom].vc_avg - cons[tid].vc_avg;
		float imitation_prob = float(dv > 0);	// no imitation for 0 payoff difference

		if (curand_uniform(&RNG_states[tid]) <= imitation_prob) { 

			float h_new    = cons[id_whom].h    + mu_h*curand_normal(&RNG_states[tid]);	
			float RT_new   = cons[id_whom].RT   + mu_RT*curand_normal(&RNG_states[tid]);	
			float Kdsd_new = cons[id_whom].Kdsd + mu_kd*curand_normal(&RNG_states[tid]);	

			if (b_ih)  cons_child[tid].h    = clamp(h_new, 0.f, h_new);	
			if (b_irt) cons_child[tid].RT   = clamp(RT_new, 0.f, RT_new);	
			if (b_ikd) cons_child[tid].Kdsd = clamp(Kdsd_new, 0.f, Kdsd_new);	
		}
	}

}


void ConsumerSystem::imitate_by_kernel_sync(){
	
	// calculate pairwise distances
	dim3 nt(16,16);
	dim3 nb((nc-1)/16+1, (nc-1)/16+1); 

	cudaMemcpy2D(&consumers_dev[0].wsum, sizeof(Consumer), nc_zeros, sizeof(float), sizeof(float), nc, cudaMemcpyHostToDevice); 
	//cudaMemcpy(wsum_dev, nc_zeros, nc*sizeof(float), cudaMemcpyHostToDevice);	// reset wsum to zero
	pairwise_distance_kernel <<<nb, nt >>> (consumers_dev, pd_dev, nc, ki_sd, L, dL);
	getLastCudaError("pd_kernel");

	// calculate whom to imitate
	if (ki_sd < 15){	// use roulette sampling
		sample_roulette_kernel <<<(nc-1)/64+1, 64 >>> (nc, pd_dev, cumm_prob_dev, consumers_dev, cs_dev_XWstates);
		getLastCudaError("sample_roulette kernel");
	}
	else{				// use rejection sampling
		sample_reject_kernel <<<(nc-1)/64+1, 64 >>> (nc, pd_dev, consumers_dev, NULL, cs_dev_XWstates);
		getLastCudaError("sample_reject kernel");
	}

	// imitate
	cudaMemcpy(consumers_child_dev, consumers_dev, nc*sizeof(Consumer), cudaMemcpyDeviceToDevice);	
	imitate_sync_kernel <<< (nc-1)/64+1, 64 >>> (consumers_dev, consumers_child_dev, cs_dev_XWstates, nc, rImit, dt,
										b_imit_h, b_imit_rt, b_imit_kd,
										mu_h, mu_RT, mu_kd);
	getLastCudaError("imitate sync kernel");
	cudaMemcpy(consumers_dev, consumers_child_dev, nc*sizeof(Consumer), cudaMemcpyDeviceToDevice);	
		
}










/*
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

__global__ void imitate_global_sync_kernel(Consumer* cons, Consumer* cons_child, curandState * RNG_states, int nc, float rImit, float dt,
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

			if (b_ih)  cons_child[tid].h    = clamp(h_new, 0.f, h_new);	
			if (b_irt) cons_child[tid].RT   = clamp(RT_new, 0.f, RT_new);	
			if (b_ikd) cons_child[tid].Kdsd = clamp(Kdsd_new, 0.f, Kdsd_new);	
		}
	}

}


void ConsumerSystem::imitate_global(){
	
	int nt = min(256, nc); int nb = (nc-1)/nt+1;
	imitate_global_kernel <<< nb, nt >>> (consumers_dev, cs_dev_XWstates, nc, rImit, dt,
										  b_imit_h, b_imit_rt, b_imit_kd);
	getLastCudaError("imitate global kernel");
		
}


void ConsumerSystem::imitate_global_sync(){
	
	int nt = min(256, nc); int nb = (nc-1)/nt+1;

	cudaMemcpy(consumers_child_dev, consumers_dev, nc*sizeof(Consumer), cudaMemcpyDeviceToDevice);	

	imitate_global_sync_kernel <<< nb, nt >>> (consumers_dev, consumers_child_dev, cs_dev_XWstates, nc, rImit, dt,
											   b_imit_h, b_imit_rt, b_imit_kd);
	getLastCudaError("imitate global sync kernel");
	
	cudaMemcpy(consumers_dev, consumers_child_dev, nc*sizeof(Consumer), cudaMemcpyDeviceToDevice);	
}



__global__ void find_nn_kernel(Consumer * cons, int nc, float L, float dL){
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= nc) return;

	cons[tid].nn = -1;
	cons[tid].nn_dist = 1e20;
	for (int i=0; i<nc; ++i){

		if (i==tid) continue;	// skip comparison with self
		
		float2 v2other = periodicDisplacement(	make_float2(cons[tid].pos_i.x*dL+dL/2, cons[tid].pos_i.y*dL+dL/2), 
												make_float2(  cons[i].pos_i.x*dL+dL/2,   cons[i].pos_i.y*dL+dL/2), 
												L, L );
		float d2other = length(v2other);

		float Inn = d2other < cons[tid].nn_dist;

		cons[tid].nn_dist = fminf(cons[tid].nn_dist, d2other);
		cons[tid].nn = Inn*i + (1-Inn)*cons[tid].nn;
		
	}		
}



__global__ void imitate_local_sync_kernel(Consumer* cons, Consumer* cons_child, curandState * RNG_states, int nc, float rImit, float dt,
									  bool b_ih, bool b_irt, bool b_ikd){
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= nc) return;
	
	bool b_imit = curand_uniform(&RNG_states[tid]) < rImit*dt;
	if (b_imit){
	
		int id_whom = cons[tid].nn;
		
		float dv = cons[id_whom].vc_avg - cons[tid].vc_avg;
		float imitation_prob = float(dv > 0);

		if (curand_uniform(&RNG_states[tid]) <= imitation_prob) { 

			float h_new    = cons[id_whom].h    + 0.02*curand_normal(&RNG_states[tid]);	
			float RT_new   = cons[id_whom].RT   + 1.00*curand_normal(&RNG_states[tid]);	
			float Kdsd_new = cons[id_whom].Kdsd + 0.20*curand_normal(&RNG_states[tid]);	

			if (b_ih)  cons_child[tid].h    = clamp(h_new, 0.f, h_new);	
			if (b_irt) cons_child[tid].RT   = clamp(RT_new, 0.f, RT_new);	
			if (b_ikd) cons_child[tid].Kdsd = clamp(Kdsd_new, 0.f, Kdsd_new);	
		}
	}

}


void ConsumerSystem::imitate_local_sync(){
	
	int nt = min(256, nc); int nb = (nc-1)/nt+1;

	find_nn_kernel <<< nb, nt >>> (consumers_dev, nc, L, dL);
	getLastCudaError("find nn kernel");

	//printConsumers();

	cudaMemcpy(consumers_child_dev, consumers_dev, nc*sizeof(Consumer), cudaMemcpyDeviceToDevice);	

	imitate_local_sync_kernel <<< nb, nt >>> (consumers_dev, consumers_child_dev, cs_dev_XWstates, nc, rImit, dt,
											  b_imit_h, b_imit_rt, b_imit_kd);
	getLastCudaError("imitate global sync kernel");
	
	cudaMemcpy(consumers_dev, consumers_child_dev, nc*sizeof(Consumer), cudaMemcpyDeviceToDevice);	
}
*/

void ConsumerSystem::freeMemory(){
	cudaFree(consumers_dev);
	cudaFree(consumers_child_dev);
	delete [] vc_window;
	cudaFree(vc_window_dev);
	delete [] ke;
	cudaFree(ke_dev);
	delete [] ke_all;
	cudaFree(ke_all_dev);
	cudaFree(cs_dev_XWstates);
	delete [] cs_seeds_h;
	cudaFree(cs_seeds_dev);
	
	if (graphics){
		cons_shape.deleteShaders();
		cons_shape.deleteVBO();
	}
}



void ConsumerSystem::printConsumers(){

	cudaMemcpy(&consumers[0], consumers_dev, sizeof(Consumer)*nc, cudaMemcpyDeviceToHost);

	for (int i=0; i<nc; ++i){
		cout << consumers[i].pos_i.x*dL+dL/2 << " " << consumers[i].pos_i.y*dL+dL/2 << " " << endl;
	}
	cout << endl;

	for (int tid=0; tid<nc; ++tid){
	cout << tid << ": ";
		for (int i=0; i<nc; ++i){

	//		if (i==tid) continue;	// skip comparison with self
		
			float2 v2other = periodicDisplacement(	make_float2(consumers[tid].pos_i.x*dL+dL/2, consumers[tid].pos_i.y*dL+dL/2), 
													make_float2(  consumers[i].pos_i.x*dL+dL/2,   consumers[i].pos_i.y*dL+dL/2), 
													L, L );
			
			cout << v2other.x << " " << v2other.y << " | ";
			float d2other = length(v2other);

			cout << d2other << " ";
	
//			float Inn = d2other < cons[tid].nn_dist;

//			cons[tid].nn_dist = fminf(cons[tid].nn_dist, d2other);
//			cons[tid].nn = Inn*i + (1-Inn)*cons[tid].nn;
		
		}		
		cout << "\n";
	}	

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


