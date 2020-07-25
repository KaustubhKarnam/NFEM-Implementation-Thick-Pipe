### Kaustubh Karnam, 65167 ###
### Computational Material Science, TU Bergakadamie Freiberg ###
### Nonlinear Finite Element Methods; Assignment - SS2020 ###
### Lecturer in charge: Dr. Geralf Hütter ###
################################################################


'''This segment defines the various input parameters that need to be fed to the other segments of the code'''


################################################################
import numpy as np
import math
import matplotlib.pyplot as plt
################################################################

inner_radius = 40		#Inner radius of cylinder
outer_radius = 80 		#Outer radius of cylinder
E = 70000			#Young's Modulus
mu = 0.25			#Poisson ratio
xi = 0				#Value of ξ
weight = 2			#Value of weight in Gauss Quadrature
p_max = 50			#Maximum pressure
ts = 0				#Start time
tl = 2				#Time after which pressure is constant
tf = 10				#Final time
t_step = 0.01			#Time step
T = 1				#Characteristic time scale 
num_ele = 10 			#Number of elements
Q = 35000			#Modulus Q
iteration_val = 20		#Number of iterations of Newton-Raphson scheme