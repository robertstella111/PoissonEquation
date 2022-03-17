# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 23:36:58 2021

@author: rober
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import constants
from matplotlib.colors import LogNorm
import sys
import os
import re


xvalue = []
yvalue = []
zvalue = []
n = 0


file_names = os.listdir()
fdata = open("TwoD.txt", 'r')
if fdata.mode == 'r':
    line = fdata.readlines()
    for x in line:
        if n == 0:
            N1 = int(x.split("|")[0])
            N2 = int(x.split("|")[1])
            n = n + 1
        else :
            xvalue.append(float(x.split("|") [0]))
            yvalue.append(float(x.split("|") [1]))
            zvalue.append(float(x.split("|") [2]))
    fdata.close()
            #Function.ausgabe(xvalue, yvalue, n) für maganetica
            
    
    #xvalue = np.multiply(xvalue, a)  
    #yvalue = np.subtract(yvalue, 0.8)
print(N1)
print(N2)
z = np.reshape(zvalue,(N1,N2))
x = np.reshape(xvalue,(N1,N2))
y = np.reshape(yvalue,(N1,N2))



#c = plt.pcolormesh(xvalue, yvalue, zvalue, cmap ='Greens', vmin = -1, vmax = 1)
#plt.colorbar(c)
  
#plt.title('matplotlib.pyplot.pcolormesh() function Example', fontweight ="bold")
#plt.show()
#z = np.cos(xvalue) * np.sin(yvalue) 
#print(z)
cmap = plt.get_cmap('RdYlBu',10)
fig,ax=plt.subplots(1,1)
cp = ax.contourf(x, y, z,  cmap=cmap)
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('Poisson 2D')
ax.set_xlabel('x')
ax.set_ylabel('y')

plt.savefig("Potential2D.pdf")