# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 15:18:33 2017

@author: e.nguyen-van
"""

import numpy as np
import math
import ATRgeometry


g=ATRgeometry.data(1,0) # arg = Vtsize and inop engines

F1=2.00589106e+04
Fx=10667.774810724921
x=np.array([  4.02189992e-02,   4.53169186e-19,   5.66339329e-19,
        -3.98226012e-19,   8.72664626e-02,   5.84206903e-02,
        -8.90046910e-02,   3.06878749e-02,   8.25018012e-01,
         8.25018012e-01])

fix=np.array([110, -10/180*math.pi, 0, 0]) # omega=Vcos(gamma)/R

V=fix[0]
alpha=x[0]
beta=fix[1]
gamma=fix[2]
omega=fix[-1]
p=x[1]
q=x[2]
r=x[3]
phi=-x[4]
theta=x[5]

sinbank=np.sin(theta)*np.cos(alpha)*np.sin(beta) + np.cos(beta)*np.cos(theta)*np.sin(phi)-np.sin(alpha)*np.sin(beta)*np.cos(theta)*np.cos(phi)
#sinbank=0
A=[0,0]
A[0]=(p*np.sin(alpha) - r*np.cos(alpha))+g.m*9.81*sinbank/(g.m*V) + F1/(g.m*V)-Fx*np.cos(alpha)*np.sin(beta)/(g.m*V)

A[1]=-np.sin(gamma)+np.cos(alpha)*np.cos(beta)*np.sin(theta)-np.sin(beta)*np.sin(phi)*np.cos(theta)-np.sin(alpha)*np.cos(beta)*np.cos(phi)*np.cos(theta)

print("A:")
print(A)