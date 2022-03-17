import matplotlib.pyplot as plt
import sys
import os
import re
import numpy as np

xvalue = []
yvalue = []
name = 1
zaehler = 0

begin = False
ende = False
typ = False
name = "nF"

fdata = open("Daten.txt", 'r')
if fdata.mode == 'r':
    line = fdata.readlines()
    for x in line:
               
        yvalue.append(float(x.split("|") [1]))
        xvalue.append(float(x.split("|") [0]))
    fdata.close()
            #Function.ausgabe(xvalue, yvalue, n) für maganetica
            
    
    #xvalue = np.multiply(xvalue, a)  
    #yvalue = np.subtract(yvalue, 0.8)

plt.plot(xvalue, yvalue, 'green', label = 'Potential', linewidth = 0.5,)
plt.title('1D Poissongleichung')   


plt.xlabel('Abstand')
plt.grid( linestyle='dotted', linewidth=0.5)
plt.legend()
plt.savefig("OneDPoissongleichung" + ".pdf")
plt.close()

