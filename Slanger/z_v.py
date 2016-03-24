# This program calculates the relative velocity between an object at z1 and an object at z0.
# In this case, z0 is 0.373189359903
# if you do not have the statistics package installed, run:
# $pip install statistics
# in the terminal


import numpy as np
import matplotlib.pyplot as plt
import math
from math import exp 
import statistics as st
from scipy.optimize import curve_fit

def z_to_v(Max_z):

 PI = 3.1415927
 Omega  = 0.3
 Lambda = 0.7
 DH     = 3000.0/0.7  #units of h^-1   DH = 3000/h, h=0.7
 N_z    = 10000
 incr_z = Max_z/float(N_z)
 z_int  = 0.0
 Dc     = 0.0
 for i_z in range (1,N_z):
         z_int = z_int + incr_z
         E_z=math.sqrt((((1.0+z_int)**3)*Omega) + Lambda)
         Dc = Dc + DH*incr_z/E_z

 v=Dc*70

 return v;


v=[0]*4000
z=[0]*4000

for i in range(36000,40000):
 j=i/100000.
 z[i-36000]=j
 v[i-36000]=z_to_v(j)-z_to_v(0.373189359903)
 if (v[i-36000]>1861.9):
  z2=j
  break
print z2

plt.plot(z,v)
plt.ylabel('km/s velocity')
plt.xlabel('redshift')
plt.show()
