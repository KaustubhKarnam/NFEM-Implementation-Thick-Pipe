### Kaustubh Karnam, 65167 ###
### Computational Material Science, TU Bergakadamie Freiberg ###
### Nonlinear Finite Element Methods; Assignment - SS2020 ###
### Lecturer in charge: Dr. Geralf Hütter ###
################################################################


'''Material Routine which needs to be called from the segment "Element_routine". 
This routine takes in the inputs - strain, partial differentiation of strain and previous overstresses
It also needs the input parameters from segment "input_parameters"
It returns:
1) Updated material tangent stiffness matrix (Ct)
2) Updated stress (stress)
3) Sigma overstress (sig_ov)'''


################################################################
import numpy as np
import math
import matplotlib.pyplot as plt
from input_parameters import *
################################################################


def material_routine(strain,delta_strain,pre_sig_ov):
    
    # 1) Stiffness due to elastic part
    A = E/((1+mu)*(1-(2*mu)))
    C = np.array([[(A*(1-mu)),(A*mu)],[(A*mu),(A*(1-mu))]])  
    # 2) Stiffness due to overstress part  
    t = 1 + t_step/T                           
    C_ov = (Q/t)*(np.array([[2/3,-1/3],[-1/3,2/3]]))

    # 3) Tangential stiffness of material (Ct)
    # 4) Deviator of strain rate (dev_strain_rate)
    # 5) Overstress using Euler backward method (sig_ov)
    # 6) Updated stress (stress)
    Ct = C + C_ov                                           
    dev_strain_rate = np.matmul(np.array([[2/3,-1/3],[-1/3,2/3]]),delta_strain)
    sig_ov = (pre_sig_ov/t)+(Q*dev_strain_rate)/t   
    stress = np.matmul(C,strain)+sig_ov 
    
    # Returning -  
    # 1) Updated material tangent stiffness matrix (Ct)
    # 2) Updated stress
    # 3) Updated overstress
    return stress,sig_ov,Ct 