#   Generate list of position of nodes according to a geometric series
#   for assignement in "Nonlinear Finite Element Methods" 
#   in summer term 2020
#   lecturer in charge: Dr. Geralf Hütter

################################################################


'''Mesh generator code snippet reference from the appendix of the assignment pdf'''


################################################################
from matplotlib.pyplot import plot,subplots,savefig
from numpy import zeros, array,transpose
import numpy as np
from input_parameters import *
#######################################################################

## Input parameters
a = inner_radius # Inner radius
b = outer_radius # Outer radius
nelem = 10 # number of elements
meshrefinementfactor = 2 
# ratio of element sizes at outer and inner radius
# ratio between element sizes of subsequent elements for a geometric series
q = meshrefinementfactor**(1.0/(nelem-1))

# size of first interval
dr = (b-a)*(1-q)/(1-meshrefinementfactor*q)
rnode = a
rnodes = [a,]

# loop over all elements
for i in range(nelem):
    rnode = rnode + dr
    rnodes.append(rnode)
    dr = dr*q

# visulaizing location of nodes
rnodes = array(rnodes)
array = transpose(rnodes)

"Uncomment the following lines for the plot"

# print(array)
# fig,ax = subplots()
# ax.plot(rnodes,zeros(nelem+1),'x',label="nodes")
# ax.set(xlabel='r')
# ax.legend()
# savefig('Mesh.png')