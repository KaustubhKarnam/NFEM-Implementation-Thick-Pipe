### Kaustubh Karnam, 65167 ###
### Computational Material Science, TU Bergakadamie Freiberg ###
### Nonlinear Finite Element Methods; Assignment - SS2020 ###
### Lecturer in charge: Dr. Geralf Hütter ###
################################################################


'''Element Routine which needs to be called from the segment "main". 
This routine takes in the location of nodes in GCS (from the "mesh_generator"), 
global displacements of elements, increment in displacement, overstress value of last step
It also needs the input parameters from segment "input_parameters"
It returns:
1) Element stiffness matrix of current element
2) Internal element force vector
3) Updated sigma (from "material_routine")
4) Updated overstress (from "material_routine")'''


################################################################
import numpy as np
import math
import matplotlib.pyplot as plt
from material_routine import *
from mesh_generator import *
from input_parameters import *
################################################################

def element_routine(rnodes,global_u,delta_u,sig_ov):
    
    # Calculating the
    # 1) Previous overstress = current overstress (pre_sig_ov)
    # 2.1) Derivative of shape function (derivative)
    # 2.2) Jacobian (constant in this case) (J)
    # 3) Shape function (N)
    # 4) Strain displacement (B)
    pre_sig_ov = sig_ov
    derivative = [-1/2,1/2] 
    J = np.matmul(derivative , np.transpose([rnodes[0],rnodes[1]]))
    N = np.array(([0.5*(1-xi),0.5*(1+xi)]))
    B = np.array(([[-1/(2*J),1/(2*J)],[1/(rnodes[0]+rnodes[1]),1/(rnodes[0]+rnodes[1])]]))
    
    
    # Calculating the
    # 5) Strain (strain)
    # 6) Delta Strain (delta_strain)
    # 7) ~ Return part of segment "material_routine" by passing 5,6 ~
    # 8) Internal element force vector (F_int_ele)
    # 9) Element stiffness matrix (K_ele) .
    strain = np.matmul(B,global_u)
    delta_strain = np.matmul(B,delta_u)
    sigma, sig_ov, Ct = material_routine(strain,delta_strain,pre_sig_ov)
    F_int_ele = weight*np.matmul(np.transpose(B),sigma)*np.matmul(np.transpose(N),rnodes)*J
    K_ele = weight*np.matmul(np.transpose(B),np.matmul(Ct,B))*np.matmul(N,rnodes)*J
    
    
    # Returning -  
    # 1) Element stiffness matrix of current element
    # 2) Internal element force vector
    # 3) Updated sigma (from "material_routine")
    # 4) Updated overstress (from "material_routine")
    return K_ele, F_int_ele, sigma, sig_ov