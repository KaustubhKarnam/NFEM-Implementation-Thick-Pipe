### Kaustubh Karnam, 65167 ###
### Computational Material Science, TU Bergakadamie Freiberg ###
### Nonlinear Finite Element Methods; Assignment - SS2020 ###
### Lecturer in charge: Dr. Geralf Hütter ###
################################################################


'''This segment generates plots for the verification and result requested in the assignment
in section 3.3 and 3.4

The following studies are made
Section 3.3 (Verification)
a) Convergence study with Q = 0, convergence achived in single iteration and towards the exact solution.
b) Convergence study for the visco-elastic case with time increments, convergence towards exact solution.

Section 3.4 (Result)
a) Distributions of ur, sigma rr and sigma ɸɸ at final time (t = tf)
b) Time history of the widening of the pipe ur (r = b,t) for t ∈ [0, tf].
The visco-elastic solution relaxes towards the elastic solution

Miscellaneous
a) Plot of the nodes from code snippet in Section "Appendix" of assignment'''


################################################################
import numpy as np
import math
import matplotlib.pyplot as plt
from main import *
from mesh_generator import *
from Element_Routine import *
from input_parameters import *
from analytical_solution import *
from material_routine import *
################################################################

## Verification
## Section 3.3 
## a) The value of Q was set to zero and the interation value set to 1

fig,ax1 = plt.subplots()
ax1.plot(rnodes,u_glo,color = "red", marker="o", markerfacecolor="black", label="Displacement for Q = 0")
ax1.set(xlabel = 'r (mm)',ylabel = 'Displacement $u_r$ (mm)',title = 'Verification 3.3 - a \n Distribution of $u_r$ with Q = 0')
ax1.legend()
#plt.show()
plt.savefig("Verification 3.3 a.png")

## Section 3.3 
## b) Convergence study for the visco-elastic case with time increments

fig,ax2 = plt.subplots()
ax2.plot(rnodes,u_glo,color = "blue", marker="o", markerfacecolor="black", label="Displacement for visco-elastic case")
ax2.set(xlabel = 'r (mm)',ylabel = 'Displacement $u_r$ (mm)',title = 'Verification 3.3 - b \n Distribution of $u_r$ for visco-elastic case')
ax2.legend()
#plt.show()
plt.savefig("Verification 3.3 b.png")

################################################################

## Result
## Section 3.4 
## a) Distributions of ur, sigma rr and sigma ɸɸ at final time (t = tf)

fig = plt.figure()
ax4 = fig.add_subplot(221)
ax4.plot(rnodes,u_glo,color = "red", marker="o", markerfacecolor="black", label="Distribution of $u_r$")
ax4.set(xlabel = 'r (mm)',ylabel = 'Displacement $u_r$ (mm)',title = 'Distribution of $u_r$ at t = $t_f$')
ax4.legend()
ax5 = fig.add_subplot(222)
ax5.plot(rnodes[1:],stress_array[:,0],color = "green", marker="o", markerfacecolor="black", label="Distribution of $σ_{rr}$")
ax5.set(xlabel = 'r (mm)',ylabel = 'Stress (Mpa)',title = 'Distribution of $σ_{rr}$ at t = $t_f$')
ax5.legend()
ax6 = fig.add_subplot(223)
ax6.plot(rnodes[1:],stress_array[:,1],color = "blue", marker="o", markerfacecolor="black", label="Distribution of $σ_{ɸɸ}$")
ax6.set(xlabel = 'r (mm)',ylabel = 'Stress (Mpa)',title = 'Distribution of $σ_{ɸɸ}$ at t = $t_f$')
ax6.legend()
plt.subplots_adjust(wspace=0.36, hspace=0.45)
#plt.show()
plt.savefig("Result 3.4 a.png")

## Section 3.4 
## b) Time history of the widening of the pipe ur (r = b,t) for t ∈ [0, tf].
duration = np.arange(start = t_step,stop =tf+0.01,step = t_step)
fig,ax7 = plt.subplots()
ax7.plot(duration,u_array,color = "blue", label="Displacement evolution")
ax7.set(xlabel = 'time (s)',ylabel = 'Displacement $u_r$ (mm)',title = 'Result 3.4 - b \n Time history of the widening of the pipe $u_r$ (r = b,t) for t ∈ [0, $t_f$]')
ax7.scatter(duration[-1], ana_disp_array[-1], label = "Analytical Solution", marker = 'o', color ='black')
ax7.legend()
#plt.show()
plt.savefig("Result 3.4 b.png")