# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Input parameters for Simulations
# If you change the order of parameters below, you will get what you deserve
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


> STRINGS 
# # > DIR
# Directories for data output
homeDir			/home/jaideep/austria_project/gpu_codes/output		# home dir - no spaces allowed
outDir  		dispersal_response_smoothness	    # output dir name
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
h0				0.1
RT0				15
kdsd0			2

nc				512
Ke_sd			4
Ke_cutoff		6

Ki_sd			300 

payoff_Tw		20
out_Tw 			100

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

mu_h			0.02
mu_RT			1
mu_kd			0.2
imresv			100

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
 0.0058941034	
-1

rimitvec  .02  -1

chvec 0.0681 -1

tmuvec  1 -1

kivec  1000
 -1

irvvec 
0 -1

drvvec
10.0000000 0.1000000  0.1211528  0.1467799  0.1778279  0.2154435  0.2610157  0.3162278
0.3831187  0.4641589  0.5623413  0.6812921  0.8254042  1.0000000  1.2115277
1.4677993  1.7782794  2.1544347  2.6101572  3.1622777  3.8311868  4.6415888
5.6234133  6.8129207  8.2540419 
-1

# 1.000000    1.333521    1.778279    2.371374    3.162278    4.216965  
# 5.623413    7.498942   10.000000   13.335214   17.782794   23.713737
# -1


# 0.001 0.006 0.03 0.2

# 0.0002	0.0002651423	0.0003515021	0.0004659904	0.0006177687	
# 0.000818983	0.0010857351	0.0014393713	0.001908191	0.0025297104	
# 0.0033536659	0.004445993	  0.0058941034	0.0078138799	0.0103589494
# 0.0137329769	0.0182059636	0.0241358528	0.0319971744	0.0424190178
# 0.056235374	0.0745518744	0.0988342672	0.1310257114	0.1737022748	


# 0.001000000  0.001325711  0.001757511 
# 0.002329952  0.003088844  0.004094915 
# 0.005428675  0.007196857  0.009540955 
# 0.012648552  0.016768329  0.022229965 
# 0.029470517 -1

# 0.001000000 .001151395 0.001325711 0.001526418 0.001757511 0.002023590
# 0.002329952 0.002682696 0.003088844 0.003556480 0.004094915 0.004714866
# 0.005428675 0.006250552 0.007196857 0.008286428 0.009540955 0.010985411
# 0.012648552 0.014563485 0.016768329 0.019306977 0.022229965 0.025595479
# 0.029470517 -1

# 1.154782  1.333521  1.539927  1.778279  2.053525
# 2.371374  2.738420  3.162278  3.651741  4.216965  4.869675
# 5.623413  6.493816  7.498942  8.659643 10.000000	-1

# 0.0002	0.0002651423	0.0003515021	0.0004659904	0.0006177687	
# 0.000818983	0.0010857351	0.0014393713	0.001908191	0.0025297104	
# 0.0033536659	0.004445993	0.0058941034	0.0078138799	0.0103589494	
# 0.0137329769	0.0182059636	0.0241358528	0.0319971744	0.0424190178	-1

# 0.1000000  0.1211528  0.1467799  0.1778279  0.2154435  0.2610157  0.3162278
# 0.3831187  0.4641589  0.5623413  0.6812921  0.8254042  1.0000000  1.2115277
# 1.4677993  1.7782794  2.1544347  2.6101572  3.1622777  3.8311868  4.6415888
# 5.6234133  6.8129207  8.2540419 10.0000000

