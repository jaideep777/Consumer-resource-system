#ifndef CONSUMERS_H
#define CONSUMERS_H

#include <string>
#include <curand.h>
#include <curand_kernel.h>

#include "../headers/graphics.h"
#include "../utils/simple_initializer.h"
using namespace std;


class Consumer{
	public:
	float2 pos;		// actual pos (float) 
	int2 pos_i;		// grid cell id of pos

	float h;
	float RT;
	float Kdsd;
	
	float rc;
	
	float ld;
	float nd;
	float vc;
	float vc_x;
	float vc_avg;
	
};

class ConsumerSystem{

	public:
	vector <Consumer> consumers;
	Consumer * consumers_dev, * consumers_child_dev;
	
	bool graphics;
	
	int nx, ny;
	float L, dL;
	int nc;
	float dt;
	
	float ke_lmax;		// bound for exploitation kernel in length units 
	float ke_nmax;		// bound for exploitation kernel in index units (kernel will go from [-ke_n...0...ke_n])
	float ke_sd;
	
	float *ke, *ke_dev;			// exploitation kernels on grid
	float *ke_all, *ke_all_dev;
	
	int vc_Tw;
	float * vc_window, *vc_window_dev; //, * vc_dev;
	
	float b, cd, ch;
	
	bool b_imit_h, b_imit_rt, b_imit_kd;
	float rImit;
	
	curandState * cs_dev_XWstates;
	int *cs_seeds_h, *cs_seeds_dev; 

	string output_dir, expt_desc;
	ofstream fout_h[1], fout_rc[1], fout_sd[1];

	PointSet cons_shape;
	

	
	public:
	void init(Initializer &I);
	void initIO(Initializer &I);
	void closeIO();
	void initRNG();
	
	void updateExploitationKernels();
	void calcResConsumed(float * resource_grid);
	void disperse(float * resource);
	void calcPayoff(int t);
	void calcAvgPayoff();

	void imitate_global();
	void imitate_global_sync();

	void writeState(int istep);

	void graphics_updateArrays();
	
	void freeMemory();
	
};

#endif
