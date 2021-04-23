#!/usr/bin/env python
# coding: utf-8

# In[2]:


#Import desired functions
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize
plt.rcParams['axes.facecolor']='white'
plt.rcParams['savefig.facecolor']='white'


# In[58]:


#Initialise values of constants
a = 0.01 #Radius of solenoid
I = 1000 #Magnitude of steady current
mu0 = 4*np.pi*10**-7

#Set the extent of the coordinate axes.
axis_size = 0.025

#Generate a grid of x and y values and an associated grid of radial distance values.
x0 = np.linspace(-axis_size,axis_size,1000)
y0 = np.linspace(-axis_size,axis_size,1000)
x,y = np.meshgrid(x0,y0)
r = np.zeros_like(x)
r = np.sqrt(x**2+y**2)

#Define a function for returning a grid of values of the magnetic field magnitude.
def solenoid_B_field(r,a,n,I,e,mu0):
    
    B_circ = mu0*I*(1/n)/(2*np.pi*r*np.sqrt(4*((np.pi*a)**2)+1/(n**2)))

    B_z = (2*np.pi*mu0*I*a*n)/(np.sqrt(((2*np.pi*a)**2)+1/(n**2)))

    B_in = r<a
    
    B = np.zeros_like(r)

    for i in range(B.shape[0]):
        for j in range(B.shape[1]):
            if B_in[i,j]:
                B[i,j] = B_z
            else:
                B[i,j] = B_circ[i,j]
                
    return B

#Generate plots of the magnetic field for various numbers of turns per unit length.
nvals = np.array((np.arange(250, 5, -5)).tolist() + (np.arange(5,0, -0.25)).tolist())
count = 100
for n in nvals:
    plt.imshow(solenoid_B_field(r,a,n,I,e,mu0), vmin=0, vmax=0.013,extent=[-axis_size, axis_size, -axis_size, axis_size])
    plt.colorbar()
    plt.text(-0.025,-0.035,'$n = {}, I = {}$'.format(str(n), str(I)+' A'))
    plt.xlabel('x [a]')
    plt.ylabel('y [a]')
    plt.savefig('snd'+str(count)+'.png', bbox_inches='tight', dpi=150)
    plt.show()
    count+=1

