#include <iostream>
#include <curand.h>
#include <cstdlib>
#include <unistd.h>
using namespace std;

//#include "../headers/graphics.h"
#include "../headers/consumers.h"
#include "../headers/resource.h"

#include "../utils/cuda_device.h" 
#include "../utils/simple_initializer.h" 
#include "../utils/simple_math.h"

curandGenerator_t generator_host;

//#define SWEEP_LOOP(x) 	for (int i_##x =0; i_##x < x##_sweep.size(); ++i_##x)

int istep = 0;

int launchSim(){
	

	return 0;
}


int main(int argc, char **argv){

	// select device
	int cDevice;
	cDevice = initDevice(argc, argv);

	// read execution parameters
	string config_filename = "../exec_config/execution_config.r";
	if (argc >2) config_filename = argv[2];

	Initializer I(config_filename);
	I.readFile();
	I.printVars();
	

	bool graphics = I.getScalar("graphicsQual")>0;
	if (graphics) init_hyperGL(&argc, argv, I);
	

	ResourceGrid * resGrid = new ResourceGrid;

	ConsumerSystem * csys = new ConsumerSystem;

	vector <float> rimitvec = I.getArray("rimitvec");
	vector <float>     bvec = I.getArray("bvec");
	for (int iri=0; iri<rimitvec.size(); ++iri){
	for (int ib=0; ib<bvec.size(); ++ib){

		resGrid->init(I);
		if (graphics) glRenderer->addShape(&resGrid->res_shape);

		csys->init(I);

		// re-init scan parameter 
		csys->b = bvec[ib];
		csys->rImit = rimitvec[iri];
		
		csys->initIO(I);

		csys->updateExploitationKernels(); 
		if (graphics) glRenderer->addShape(&csys->cons_shape);

		// launch sim
		int nsteps = I.getScalar("nsteps");
		SimpleProgressBar prog(nsteps, &istep, "Diffusion");


		prog.start();
		while(1){       // infinite loop needed to poll anim_on signal.
			if (graphics) glutMainLoopEvent();


//			int i = animate();
			csys->calcResConsumed(resGrid->res_dev);
			resGrid->grow(csys->ke_all_dev);
			csys->disperse(resGrid->res_dev);
			csys->updateExploitationKernels();
			csys->calcPayoff(istep);
			csys->imitate_local_sync();
//			usleep(50e2);   // sleep for 20 ms. This dramatically reduces CPU consumption
			++istep;
//			if (graphics && istep % 1 == 0){
//				  resGrid->graphics_updateArrays();
//				  csys->graphics_updateArrays();
//			}       

			if (istep % 100 == 0){
				csys->writeState(istep); 
			}
		
			prog.update();

			if (istep == nsteps) {
				istep = 0;
				csys->closeIO();
				break;
			}
		}
		// launch sim end.

		resGrid->freeMemory();
		csys->freeMemory();

	}
	}
	
	delete csys;
	delete resGrid;
	
//	if (graphics) delete R;
	
	return 0;
}




	// create output dirs
//	if (dataOut || plotsOut || framesOut) system(string("mkdir " + outDir).c_str());
//	if (dataOut)   system(string("mkdir " + dataDir).c_str());
//	if (plotsOut)  system(string("mkdir " + outDir + "/plots").c_str());
//	if (framesOut) system(string("mkdir " + framesDir).c_str());
	
//	// allocate memory
//	allocArrays();

//	// if graphics are on, initGL
//	if (graphicsQual > 0) initGL(&argc, argv, host_params);

//	// for all the chosen parameter sweeps, init arrays and launch simulations
//	SWEEP_LOOP(mu){ 
//	SWEEP_LOOP(nm){ 
//	SWEEP_LOOP(fb){ 
////	SWEEP_LOOP(as){ 
//	SWEEP_LOOP(rg){ 
//	SWEEP_LOOP(rsb){ 

//		// set parameters
//		mu[0] 			= mu_sweep[i_mu];
//		moveStepsPerGen = nm_sweep[i_nm];
//		fitness_base 	= fb_sweep[i_fb];
////		arenaSize 		= as_sweep[i_as];
//		Rg 				= rg_sweep[i_rg];
//		Rs_base 		= rsb_sweep[i_rsb];

//		// for each parameter set, perform as many ensembles as specified in ens_sweep[]
//		SWEEP_LOOP(ens){ 
//			iEns = ens_sweep[i_ens];
//			initStateArrays();
//			launchExpt();
//		}
//			
//	}}}}}//}

