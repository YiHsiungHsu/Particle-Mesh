# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 14:49:51 2020

@author: jared
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
dt = 0.05
N  = 10
step = 11 

f = open("Particle_position.txt","r")
data = f.readlines()
f.close()
x = []
y = []
z = []

for a in range(step):
    for i in range(a*N,(a+1)*N,1 ):
         sdata = data[i].split()
         print(sdata)
         l = float(sdata[0])
         j = float(sdata[1])
         k = float(sdata[2])
         x.append(l)
         y.append(j)
         z.append(k)
    ax = plt.subplot(projection='3d')  
    ax.scatter(x, y, z, c='r')  
    ax.set_zlabel('Z')  
    ax.set_ylabel('Y')
    ax.set_xlabel('X')
    plt.draw()
    plt.savefig('time = '+ str(a*dt) +'.png')
    plt.close()
    x = []
    y = []
    z = []
