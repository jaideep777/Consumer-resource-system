## Introduction

This code accompanies the paper "Emergence of Social Inequality in Spatial Ecological Public Goods Games" by Jaideep Joshi, Åke Brännström, and Ulf Dieckmann. 

The full code includes numerous utility functions and classes to enable execution on GPU and to make it generic enough to be able to run multiple parameter configurations. However,  implementations of the core model equations (Eqs. 1-8 in the above paper) are tagged with the word "TAG". These can be viewed in one go by typing 

```
grep -rn "TAG"
``` 


## Installation prerequisites

You will need the NVIDIA CUDA-Toolkit, OpenGL libraries (Glut, GLEW), and libGSL to run this code. 

## Usage

### Compile

In `src/Makefile`, change the CUDA_INSTALL_PATH to your CUDA directory. By default, it is set to `/usr/local/cuda`. Then compile specifying the executable name, e.g.

```
make TARGET=simulate
```

### Set up simulation parameters

All simulation parameters must be listed in an "execution configuration file" in the `exec_config` folder. For a correctly formatted example file, see `exec_config/execution_config_ki10.r`. 

Fore more details on the format of the configuration files, see the documentation of the `Initializer` class in https://github.com/jaideep777/Flare.

### Simulate

After setting the simulation parameters, run the executable as follows:

```
./simulate <device_id> <config_file> 
```

where device_id is the GPU on which to run the executable (if your system has only one NVIDIA GPU, device_id=0), and <config_file> is the relative path the the execution configuration file.

	
### Plot

After a successful simulation run, the output files (histograms of r_H, sigma_D, and Rc) will be stored in the specified output folder. These histograms can be plotted using R. E.g. to plot the timeseries as in Fig. 1, set the parameters as per the output files in `plotting/plot_histograms_singleParam_timeseries.R` and run the script. 

Other R scripts process outputs from multiple simulation runs, such as a parameter scan.



