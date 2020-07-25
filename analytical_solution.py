### Kaustubh Karnam, 65167 ###
### Computational Material Science, TU Bergakadamie Freiberg ###
### Nonlinear Finite Element Methods; Assignment - SS2020 ###
### Lecturer in charge: Dr. Geralf Hütter ###
################################################################


'''This segment defines the analytical solution of the problem. 
It is used in the "main" segment for comparing the solutions and generating graphs'''


################################################################
import numpy as np
import matplotlib.pyplot as plt
import math
from input_parameters import *
from mesh_generator import *
################################################################

def analytical_solution(rnodes):
    analytical_disp = []
    for i in range(0 , num_ele+1):
        u = (1+mu)*(p_max/E)*(inner_radius**2/(outer_radius**2-inner_radius**2))*((1-2*mu)*rnodes[i] + outer_radius**2/rnodes[i])
        analytical_disp.append(u)
    return analytical_disp

ana_disp_array = analytical_solution(rnodes)
print(ana_disp_array)