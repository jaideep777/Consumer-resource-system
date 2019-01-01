# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Input parameters for Simulations
# Notations used in the paper are given in square brackets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


> STRINGS 
# # > DIR
# Directories for data output
homeDir			/home/jaideep/austria_project/gpu_codes/output		# home dir (no spaces allowed). Directory in which outputs of various experiments eill be stored
outDir  		KI_ki10		    # output dir name (relative to home dir) to store output files from current run.
exptName 	 	hom 			# "hom" for homogeneous resource, OR "het" for heterogeneous resource 
	

> SCALARS
# > GPU_CONFIG
# population
# particles 		512		# total number of consumers
# blockSize	 	256		# threads per block

# > GRAPHICS
# graphics
graphicsQual 	0			# 0 = no graphics, 1 = show graphics
dispInterval  	50 			# display interval in ms, -ve = number of mv steps to run before display
# b_anim_on 		0		  	# turn animation on immediately on start? if false (0), wait for user to press "a" to start sim

# > EXPT
# experiment properties
# b_baseline 		1			# Baseline  - is this a baseline experiment 
# b_constRg  		1			# 

# > RESOURCE GRID
nx 				450			# Resource grid size (x)
ny				450			# Resource grid size (y)
D				0			# Diffusion constant for resource
r				0.2			# [r] Intrinsic growth rate of resource
K				50			# [K] Carrying capacity of resource
L 				225			# System size
 

# > PARTICLES
# movement parameters
dt  			0.1			# [dt] Integration time step 

Ke_sd			4			# [sigma_H] Harvesting radius (exploitation radius)
Ke_cutoff		6			# 
h0				0.1			# [r_H(t=0)] Initial harvesting rate of consumers
kdsd0			2			# [sigma_D(t=0)] Initial dispersal radius of consumers
RT0				15			# [R_T(t=0)] Initial resource threshold for dispersal

nc				512			# [N] Number of consumers


# > SELECTION
# payoff	
b				0.002		# [b] Harvesting benefit
cd				0.1			# [c_D] 	Dispersal cost
ch				0.08		# [c_H] Harvesting cost
payoff_Tw		20			# [T] Window over which payoffs are averaged
out_Tw 			100			# Steps to skip before output 

# > IMIT
imitation_rate	0.02		# [r_I] Imitation rate
Ki_sd			10 			# [sigma_I] Imitation radius

imitate_h		1			# Imitate harvesting rate? (0/1) 
imitate_RT		0			# Imitate dispersal threshold? (0/1)
imitate_Kd		1			# Imitate dispersal radius? (0/1)

mu_h			0.02		# Mutation step size in harveting rate
mu_RT			1			# Mutation step size in dispersal threshold
mu_kd			0.2			# Mutation step size in disperal radius

imresv			100			# [k_I] 1/Smoothness of imitation response function
dresv			10			# [k_D] 1/Smoothness of dispersal response function

# > OUTPUT
# output
# dataOut  		1			# output all values in files?

# > SIM
nsteps			750000		# Number of timesteps to simulate



# Turbulence params
mu_t			0.1			# Controls correlation length of spatial heterogeneity. Higher mu_t -> higher correlation length (Can vary between 0.01-1.00)

nu_t			0.005		# The remaining parameters are for the turbulence engine. See cited reference for details.
xi_t 			0.05
lambda0_t 		400
dt_t 			0.1

> ARRAYS
# > PARAM_SWEEPS
# parameter sets to loop over. These values override any of the parameters set above.
# all of these MUST BE SET AFTER particle-system initialization. 
# Each array must end in -1
# e.g. bvec		0.0002 0.002 0.02 0.2 -1

# Benefit of harvesting
bvec	 
 0.0002	0.0002651423	0.0003515021	0.0004659904	0.0006177687	
 0.000818983	0.0010857351	0.0014393713	 0.001908191	0.0025297104	
 0.0033536659	0.004445993	  0.0058941034	0.0078138799	0.0103589494
 0.0137329769	0.0182059636	0.0241358528	0.0319971744	0.0424190178
 0.056235374	0.0745518744	0.0988342672	0.1310257114	0.1737022748	
	-1

# Imitation rate
rimitvec 0.02 -1 

# Cost of harvesting
chvec 0.08 2 -1

# Spatial correlation length of resource growth rate (for heterogeneous resource simulations)
tmuvec  1 -1

# Imitation Radius
kivec 
10	21.54435	46.41589	100	215.44347	464.15888	1000	316.22777	681.29207
 -1

# Imitation response smoothness
irvvec 100 -1

# Dispersal response smoothness
drvvec 10 -1


