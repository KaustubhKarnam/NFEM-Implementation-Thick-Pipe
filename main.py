### Kaustubh Karnam, 65167 ###
### Computational Material Science, TU Bergakadamie Freiberg ###
### Nonlinear Finite Element Methods; Assignment - SS2020 ###
### Lecturer in charge: Dr. Geralf Hütter ###
################################################################


'''This segment is the main code which calls the segment "Element_Routine" which includes the material routine in it.
This "main" segment has the Newton-Raphson scheme, assembles the Global stiffness matrix, 
assembles the global internal forces vector, updates the global displacement vector,
updates global stresses'''


################################################################
import numpy as np
import math
import matplotlib.pyplot as plt
from mesh_generator import *
from Element_Routine import *
from input_parameters import *
from analytical_solution import *
from material_routine import *
################################################################


## Assigning a zero vector to the following "variables" for the first iteration.
## 1) Overstress
## 2) Stress
## 3) Global displacements of nodes
## 4) Previous delta u which stores values of delta u for current iteration to be passed again
## 5) Previous overstress ~ similar to previous delta u
## 6) A list (later converted to np.array) to contain displacement history for r = outer radius
sig_ov = np.zeros((2,1), dtype = float)
sigma = np.zeros((2,1), dtype = float)
u_glo = np.zeros((nelem+1,1), dtype = float)
prev_delta_u = np.zeros((nelem+1,1), dtype = float)
prev_ov = np.zeros((2,1), dtype = float) 
u_array = [] 


################################################################

## The time condition, the Newton Raphson condition, iteration termination condition, convergence condtion are included   
## Time Condition
for t in np.arange(start = t_step,stop =tf+0.01,step = t_step):
    if t > tl:   
        p = p_max*inner_radius
    else:
        p = p_max*inner_radius*t*1/tl
    it = 0
    ## Newton Raphson Condition
    ## Iteration termination condition
    while it < iteration_val:
        ## Assigning a zero vector/matrix to the following "variables" for the each iteration.
        ## 7) Element stiffness matrix
        ## 8) Element internal force vector
        ## 9) Global stiffness matrix
        ## 10) Global internal forces vector
        ## 11) Global external forces vector
        ## 12) Displacement increment vector
        ## 13) Assigning boundary condition in the global external force vector 
        ## 14) Storing sigma rr and sigma ɸɸ in stress list (later converted to np.array)
        K_ele = np.zeros(2) 
        F_int_ele = np.zeros((2,1)) 
        K_glo = np.zeros((nelem+1,nelem+1))
        F_int_glo = np.zeros((nelem+1,1)) 
        F_ext_glo = np.zeros((nelem+1,1))
        delta_u = np.zeros((nelem+1,1))
        F_ext_glo[0] = p 
        stress_array = [] 


        ## Element loop and assembly
        for i in range(nelem):
            K_ele, F_int_ele, sigma, sig_ov = element_routine(rnodes[i:i+2],u_glo[i:i+2],prev_delta_u[i:i+2],prev_ov)
            stress_array.append(sigma)
            K_glo[i:i+2,i:i+2] = K_glo[i:i+2,i:i+2] + K_ele
            F_int_glo[i:i+2] = F_int_glo[i:i+2] + F_int_ele

        
        ## Solving equation and updating the displacements in global displacement 
        ## Assigning value to prev_delta_u
        delta_u = np.linalg.solve(K_glo,(F_ext_glo-F_int_glo))
        u_glo = u_glo + delta_u
        prev_delta_u = delta_u[:]


        ## Convergence condition
        if (np.linalg.norm(F_ext_glo - F_int_glo,np.inf)) < 0.005*(np.linalg.norm(F_int_glo,np.inf)) and (np.linalg.norm(delta_u,np.inf))< 0.005*(np.linalg.norm(u_glo,np.inf)):
            break
        else:
            it+=1

    u_array.append(u_glo[-1])
    prev_ov = sig_ov 


stress_array = np.array(stress_array) 
u_array = np.array(u_array)

# print(u_array)
# print(u_glo)
# print(stress_array)
# print(analytical_solution(rnodes))