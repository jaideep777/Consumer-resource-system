# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Input parameters for Simulations
# If you change the order of parameters below, you will get what you deserve
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


> STRINGS 
# # > DIR
# Directories for data output
homeDir			/home/jaideep/austria_project/gpu_codes/output		# home dir - no spaces allowed
outDir  		b_scan_large_set2			    # output dir name
exptName 	 	hom 								# expt name
 	

> SCALARS
# > GPU_CONFIG
# population
particles 		512		# total number of particles 4096, 16384, 65536, 262144, 1048576
blockSize	 	256		# threads / block

# > GRAPHICS
# graphics
graphicsQual 	0			# 0 = no graphics, 1 = basic graphics, 2 = good graphics, 3 = fancy graphics, charts etc
dispInterval  	50 			# display interval in ms, -ve = number of mv steps to run before display
b_anim_on 		0		  	# turn animation on immediately on start? if false (0), wait for user to press "a" to start sim

# > EXPT
# experiment properties
b_baseline 		1			# Baseline  - is this a baseline experiment (Rs = 0)
b_constRg  		1			# Grouping  - use constant radius of grouping?

# > RESOURCE GRID
nx 				450
ny				450
D				0
r				0.2
K				50
L 				225			# the size of the entire space (x & y), when body size = 1 --# > determines density
 

# > PARTICLES
# movement parameters
dt  			0.1			# time step (for movement)
h0				0.5
RT0				25
kdsd0			2

nc				512
Ke_sd			4
Ke_cutoff		6

payoff_Tw		5

# > SELECTION
# payoff
b				0.002
cd				0.1
ch				0.08

# > IMIT
imitation_rate	0.02	
imitate_h		1
imitate_RT		0
imitate_Kd		1


# > OUTPUT
# output
dataOut  		1			# output all values in files?

# > SIM
nsteps			750000


# Altruism params

# Turbulence params

> ARRAYS
# > PARAM_SWEEPS
# parameter sets to loop over. These override any of the parameters set above.
# all of these MUST BE SET AFTER particle-system initialization. 
# 	  else, particleSystem will initialize with 1st value in these vectors by default 
# c_full = 0.1 0.12 0.14 0.16 0.19 0.22 0.26 0.3 0.35 0.41 0.48 0.56 0.65 0.77 0.89 1.05 1.22 1.43 1.67 1.96 2.29 2.68 3.13 3.66 4.28 5 5.85 6.84 8 9.36 10.95 12.8  -1
# cS_full = 20 40 60 80 100 120 140 160 180 200 220 240 260 280 300 320 -1
# c = 0.1 0.14 0.19 0.26 0.35 0.48 0.65 0.89 1.22 1.67 2.29 3.13 4.28 5.85 8 10.95 -1
# c_offset = 0.12 0.16 0.22 0.3 0.41 0.56 0.77 1.05 1.43 1.96 2.68 3.66 5 6.84 9.36 12.8 -1
# bvec		0.0002 0.002 0.02 0.2 -1

bvec	 
 0.0067864435 0.0078138799 0.0089968653
 0.0103589494 0.0119272466 0.0137329769 0.0158120864 0.0182059636 0.0209622627 0.0241358528
 0.0277899099 0.0319971744 0.0368413994 0.0424190178 0.0488410619 0.0562353740 0.0647491509
 0.0745518744 0.0858386852 0.0988342672 0.1137973206 0.1310257114 0.1508624013 0.1737022748
 0.2000000000 -1


