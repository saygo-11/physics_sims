# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 10:06:24 2026

@author: User
"""

import numpy as np
import matplotlib.pyplot as plt

g = 9.8 # m/s^2
L = 2 # meters
colors = ["green", "blue", "red", "black"]
thetas = [0.4, 1, 2, 3] # rad
dt = 0.1 #sec
t_max = 20 #sec
c = 0

label = lambda t, c : str(thetas[c]) + " rad" if t == 0 else None

for theta in thetas:
   t = 0
   v = 0
   a = 0
   while t <= t_max:
        a = - (g/L) * np.sin(theta) * dt # acc. prev. acc. - the change in acceleration 
        v = v + a * dt 
        theta = theta + v * dt
        plt.scatter(t, theta, s=0.5, color=colors[c], 
                    label = label(t, c))
        t += dt
        
   c += 1
   
plt.title("Simple Pendulum Angular Position over Time")
plt.xlabel("Time(sec)")
plt.ylabel(r"$\theta(t)$")
plt.legend()
plt.show()