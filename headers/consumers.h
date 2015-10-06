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
	int2 pos_i;	// grid cell id of pos

	float h;
	float RT;
	float Kdsd;
	
	float rc;
	float V;
	
};

class ConsumerSystem{

	public:
	vector <Consumer> consumers;
	
	int nx, ny;
	float L, dL;
	int nc;
	
	float ke_lmax;		// bound for exploitation kernel in length units
	float ke_nmax;		// bound for exploitation kernel in index units
	float ke_sd;
	
	float *ke, *ke_dev;
	float *ke_all, *ke_all_dev;
	
	int2 * pos_i_dev;
	float * h_dev;
	float * rc_dev;
	float * kdsd_dev, *RT_dev;
	
	curandState * cs_dev_XWstates;
	int *cs_seeds_h, *cs_seeds_dev; 

	Shape particles_shape;
	
	public:
	void init(Initializer &I);
	void initRNG();
	
	void updateExploitationKernels();
	void calcResConsumed(float * resource_grid);
	void disperse(float * resource);
	
};

#endif
