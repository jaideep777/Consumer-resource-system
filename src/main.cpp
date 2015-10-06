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
	
	// ~~~~~~~~~~~~~~ CPU random generator (MT) ~~~~~~~~~~~~~~~~~~~~~~~~~~~
//	curandCreateGeneratorHost(&generator_host, CURAND_RNG_PSEUDO_MTGP32);
//	curandSetPseudoRandomGeneratorSeed(generator_host, time(NULL));	// seed by time in every expt


	bool graphics = I.getScalar("graphicsQual")>0;
	if (graphics) init_hyperGL(&argc, argv, I);
	

	ResourceGrid * resGrid = new ResourceGrid;
	resGrid->init(I);
	glRenderer->addShape(&resGrid->res_shape);

	ConsumerSystem * csys = new ConsumerSystem;
	csys->init(I);
	csys->updateExploitationKernels();


	// launch sim
	SimpleProgressBar prog(1000, &istep, "Diffusion");

	prog.start();
	while(1){	// infinite loop needed to poll anim_on signal.
		if (graphics) glutMainLoopEvent();
		
//		int i = animate();
		csys->calcResConsumed(resGrid->res_dev);
		resGrid->grow(csys->ke_all_dev);
		csys->disperse(resGrid->res_dev);
		csys->updateExploitationKernels();
		
		++istep;
		resGrid->graphics_updateArrays();
	
		usleep(5e2);	// sleep for 20 ms. This dramatically reduces CPU consumption
		prog.update();

		
		if (istep ==1000) {
			break;
		}
	}
	// launch sim end.


	ofstream fout("res.txt");
	
	int nx = I.getScalar("nx");
	int ny = I.getScalar("ny");
	cudaMemcpy(resGrid->res, resGrid->res_dev, resGrid->nx*resGrid->ny*sizeof(float), cudaMemcpyDeviceToHost);
	for (int j=0; j<ny; ++j){
		for (int i=0; i<nx; ++i){
			fout << resGrid->res[ix2(i,j,resGrid->nx)] << "\t";
		}
		fout << "\n";
	}
	fout << "\n"; 
	
	
//	// disconnect renderer
//	if (graphics) R->disconnect();
//	// --------------------

	resGrid->freeMemory();

//	if (graphics) delete R;
//	delete psys;
	
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

