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
N  = 2
T  = 0.05
step = T/dt + 1

f = open("Particle_position.txt","r")
data = f.readlines()
f.close()
x = []
y = []
z = []

for a in range(step):
    for i in range(a*N,(a+1)*N,1 ):
         sdata = data[i].split(" ")
         x.extend(sdata[0])
         y.extend(sdata[1])
         z.extend(sdata[2])
    ax = plt.subplot(projection='3d')  
    ax.scatter(x, y, z, c='r')  
    ax.set_zlabel('Z')  
    ax.set_ylabel('Y')
    ax.set_xlabel('X')
    plt.draw()
    plt.pause(10)
    plt.savefig('time = ',step*dt)
    plt.close()
    x = []
    y = []
    z = []