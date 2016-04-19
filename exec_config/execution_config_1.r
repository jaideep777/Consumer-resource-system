# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Input parameters for Simulations
# If you change the order of parameters below, you will get what you deserve
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


> STRINGS 
# # > DIR
# Directories for data output
homeDir			/home/jaideep/austria_project/gpu_codes/output		# home dir - no spaces allowed
outDir  		z_for_nvvp			    # output dir name
exptName 	 	het 								# expt name


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
RT0				15
kdsd0			2

nc				512
Ke_sd			4
Ke_cutoff		6

payoff_Tw		20

# > SELECTION
# payoff
b				0.001
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
mu_t			.1		# 0.1
nu_t			.005	# 0.005
xi_t			.05		# 0.05     <- 0.05 - 0.2
lambda0_t		400		# 1.0
dt_t			0.1

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
 0.0002000000 0.0002302791 0.0002651423 0.0003052836 0.0003515021 0.0004047179 0.0004659904
 0.0005365392 0.0006177687 0.0007112961 0.0008189830 0.0009429733 0.0010857351 0.0012501104
 0.0014393713 0.0016572855 0.0019081910 0.0021970823 0.0025297104 0.0029126970 0.0033536659
 0.0038613955 0.0044459930 0.0051190958 0.0058941034
  
 0.0067864435 0.0078138799 0.0089968653
 0.0103589494 0.0119272466 0.0137329769 0.0158120864 0.0182059636 0.0209622627 0.0241358528
 0.0277899099 0.0319971744 0.0368413994 0.0424190178 0.0488410619 0.0562353740 0.0647491509
 0.0745518744 0.0858386852 0.0988342672 0.1137973206 0.1310257114 0.1508624013 0.1737022748
 0.2000000000 -1

rimitvec 0.1 -1

chvec 0.08 -1

tmuvec  1 -1

# 0.01000000 0.01098541 0.01206793 0.01325711 0.01456348 0.01599859 0.01757511 0.01930698
# 0.02120951 0.02329952 0.02559548 0.02811769 0.03088844 0.03393222 0.03727594 0.04094915
# 0.04498433 0.04941713 0.05428675 0.05963623 0.06551286 0.07196857 0.07906043 0.08685114
# 0.09540955 0.10481131 0.11513954 0.12648552 0.13894955 0.15264180 0.16768329 0.18420700
# 0.20235896 0.22229965 0.24420531 0.26826958 0.29470517 0.32374575 0.35564803 0.39069399
# 0.42919343 0.47148664 0.51794747 0.56898660 0.62505519 0.68664885 0.75431201 0.82864277
# 0.91029818 1.00000000 -1

# 0.001000000  0.001325711  0.001757511 
# 0.002329952  0.003088844  0.004094915 
# 0.005428675  0.007196857  0.009540955 
# 0.012648552  0.016768329  0.022229965 
# 0.029470517 -1


# 1.154782  1.333521  1.539927  1.778279  2.053525
# 2.371374  2.738420  3.162278  3.651741  4.216965  4.869675
# 5.623413  6.493816  7.498942  8.659643 10.000000	-1

# 0.001000000 .001151395 0.001325711 0.001526418 0.001757511 0.002023590
# 0.002329952 0.002682696 0.003088844 0.003556480 0.004094915 0.004714866
# 0.005428675 0.006250552 0.007196857 0.008286428 0.009540955 0.010985411
# 0.012648552 0.014563485 0.016768329 0.019306977 0.022229965 0.025595479
# 0.029470517 

# 0.033932218 0.039069399 0.044984327 0.051794747 0.059636233 0.068664885 
# 0.079060432 0.091029818 0.104811313 0.120679264 0.138949549 0.159985872 
# 0.184206997 0.212095089 0.244205309 0.281176870 0.323745754 0.372759372 
# 0.429193426 0.494171336 0.568986603 0.655128557 0.754312006 0.868511374 
# 1.000000000 -1 




