from tkinter import *
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as anim
import time


     
#x , t = np.meshgrid(np.linspace(-3.14159*10,+3.14159*10,100), np.linspace(1,4,100))

#z = np.cos(x/t)


xvalue = []
yvalue = []
zvalue = []
n = 0


fdata = open("Diff.txt", 'r')
if fdata.mode == 'r':
    line = fdata.readlines()
    for x in line:
        if n == 0:
            T = int(x.split("|")[0])
            N = int(x.split("|")[1])
            n = n + 1
        else :
            xvalue.append(float(x.split("|") [0]))
            yvalue.append(float(x.split("|") [1]))
    fdata.close()
            #Function.ausgabe(xvalue, yvalue, n) für maganetica
            
    
    #xvalue = np.multiply(xvalue, a)  
    #yvalue = np.subtract(yvalue, 0.8)
print(N)
print(T)


x = np.reshape(xvalue,(T,N))
z = np.reshape(yvalue, (T,N))



plt.ion()

figure, ax = plt.subplots(figsize=(8,6))
line1, = ax.plot(x[0],z[0] )

for i in range(T):
   
    line1.set_xdata(x[0])
    line1.set_ydata(z[i])
    
    figure.canvas.draw()
    
    figure.canvas.flush_events()
    time.sleep(float(10/T))


