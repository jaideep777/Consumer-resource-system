## Introduction

This code accompanies the paper "Emergence of Social Inequality in Spatial Ecological Public Goods Games" by Jaideep Joshi, Ake Brannstrom, and Ulf Dieckmann.


## Installation prerequisites

You will need the NVIDIA CUDA-Toolkit, OpenGL libraries (Glut, GLEW), and libGSL to run this code. Once you have installed these libraries.

## Usage

### Compile

In src/Makefile, change the CUDA_INSTALL_PATH to your CUDA directory. By default, it is set to /usr/local/cuda. Then compile speficying the executable name. E.g.

```
make TARGET=simulate
```

### Set up simulation parameters

All simulation parameters must be listed in an "execution configuration file" in the exec_config folder. For a correctly formatted example file, see "execution_config_ki10.r". 

### Simulate

After setting the simulation parameters, run the executable as follows:

```
./simulate <device_id> <config_file> 
```

where device_id is the GPU on which to run the executable (if your system has only one NVIDIA GPU, device_id=0), and <config_file> is the relative path the the execution configuration file.


