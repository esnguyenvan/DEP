# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 14:05:08 2017

Main file for stability mapping

@author: e.nguyen-van
"""
# import all modules
import ReadCoefAero
import AeroForces
import numpy as np
import math
import scipy.linalg
import matplotlib.pyplot as plt
import ATRgeometry
import equation as e
from scipy.optimize import fsolve, minimize, root

np.set_printoptions(threshold=np.nan)
#List all .stab file from vsp aero and read the coeff

filename=['ATR72_SI_MTOW_Control_flap.stab', 'ATR72_SI_MTOW_Control_flapVmax.stab', 
'ATR72_SI_MTOW_Control_Vmin.stab','ATR72_SI_MTOW_Control_Manoeuver.stab','ATR72_SI_MTOW_Control_Cruise.stab']

# hard coded velocity and corresponding air density
Velocities=(56.1,71,71,90,130)
rho_vec=(1.225,1.225,1.225,1.225*0.7812,1.225*0.4549)

Matrix=ReadCoefAero.ReadCoef(filename)
CoefMatrix=Matrix[:,1:]
BaseAero=Matrix[:,1]

g=ATRgeometry.data(1)
V_base=80
rho_base=AeroForces.InterpolRho(V_base,rho_vec,Velocities)

#plot coefficients evolution

vel=np.linspace(55,135,20)
Cy_beta=np.array([])
CL_alpha=np.array([])
Cy_r=np.array([])
Cy_p=np.array([])
Cn_r=np.array([])
Cl_r=np.array([])

for i in range(len(vel)):
    LocalMatrix=AeroForces.CoefInterpol(vel[i], CoefMatrix, Velocities)
    Cy_beta=np.append(Cy_beta,LocalMatrix[1,1])
    CL_alpha=np.append(CL_alpha,LocalMatrix[2,0])
    Cy_r=np.append(Cy_r,LocalMatrix[1,4])
    Cy_p=np.append(Cy_p,LocalMatrix[1,2])
    Cn_r=np.append(Cn_r,LocalMatrix[-1,4])
    Cl_r=np.append(Cl_r,LocalMatrix[3,4])
    

#plt.figure(1)
#plt.plot(vel,Cy_beta*math.pi/180)
#plt.ylabel('Cy_beta (/°)')
#plt.xlabel('Vel (m/s)')
#plt.grid()
#
#plt.figure(2)
#plt.plot(vel,Cy_r)
#plt.ylabel('Cy_beta (s/rad)')
#plt.xlabel('Vel (m/s)')
#plt.grid()
#
#plt.figure(3)
#plt.plot(vel,Cn_r)
#plt.ylabel('Cy_r (s/rad)')
#plt.xlabel('Vel (m/s)')
#plt.grid()
#
#plt.figure(4)
#plt.plot(vel,CL_alpha/180*math.pi)
#plt.ylabel('CL_alpha (/°)')
#plt.xlabel('Vel (m/s)')
#plt.grid()


# Forces and equations
#x : state [u,v,w,p,q,r,phi,theta,da,de,dr,dx]
xtest=[80,1,18,0,0,0,0,0.01,0,0,0,0.5]

# --- constraints are defined as sets of additional equations to satisfy
# the list "con" is the list of imposed state, it follows the following order
#  con=[V, beta, p, q, r, phi, theta, da, de, dr, dx]
con=np.array([80,30/180*math.pi,0,0,0,0,0,0,0,0,0,0])

# ----- Test equation -----
def X(x):
    # x : state [u,v,w,p,q,r,phi,theta,da,de,dr,dx]
    # need to do some operation before computing the forces
    global CoefMatrix
    global Velocities
    global rho_vec
    global g
    global con
    
#    xreturn=np.empty((0,12))
    
    V=math.sqrt(np.dot(x[0:3],x[0:3]))
    sub_vect=np.array([math.atan(x[2]/x[0]), math.asin(x[1]/V)])
    sub_vect=np.append(sub_vect,x[3:6])
    #(alpha, beta, p, q, r, da, dr, de, dx)
    sub_vect=np.append(sub_vect,[x[-4],x[-2],x[-3],x[-1]])
#    svec_force=np.append(sub_vect,[x[-2],x[-3],x[-1]])
    
    F=AeroForces.CalcForces(V, sub_vect, CoefMatrix, Velocities, rho_vec)
#    print("Printing forces :")
#    print(F)
    #all input taken into account in force computation
    
    vel_vec=x[0:3]
    p=x[3]
    q=x[4]
    r=x[5]
    xdot=-np.cross(np.array([p,q,r]),vel_vec) + 9.81*np.array([-math.sin(x[7]), math.cos(x[7])*math.sin(x[6]), math.cos(x[7])*math.cos(x[6])])+F[0:3]/g.m
#    print(np.cross(np.array([p,q,r]),vel_vec))
    I_inv=np.array([[g.Iz/g.Ix, 0, g.Ixz/g.Ix],
    [0, 1/g.Iy*(g.Iz-g.Ixz**2/g.Ix), 0],
    [g.Ixz/g.Ix, 0, 1]])/(g.Iz-g.Ixz**2/g.Ix)
   
    rot=np.array([p,q,r])
    I=np.array([[g.Ix, 0, -g.Ixz],[0,g.Iy,0],[-g.Ixz,0,g.Iz]]) 
    Mvect=F[3:6]-np.cross(rot,np.array([g.hp,0,0]))-np.cross(rot, np.dot(I,rot))
    
    Mdot=np.dot(I_inv,Mvect);
    xdot=np.append(xdot,Mdot)
    
    # additional constraints
    xadd=np.empty((0,6))
    xadd=V-con[0]
    xadd=np.append(xadd, [math.asin(x[1]/V)-con[1]])
    xadd=np.append(xadd, [x[3:6]-con[2:5]])
    xadd=np.append(xadd, [x[6]-con[5]])
#    print("xadd:")
#    print(xadd)
#    print("xdot:")
#    print(xdot)
#    
    xreturn=np.append([xdot],[xadd])
#    print("return x")
#    print(type(xreturn))
#    print(xreturn)
    
#    return math.sqrt(np.dot(xreturn,xreturn))
    return xreturn

#print(X(xtest))
# eigen values D,V = linalg.eig(a)

# try out a solve run
x0=np.array([95,2,10,0,0,0,0,0.02,0,0,0,0.2])
arg=(CoefMatrix, Velocities, rho_vec,g)
#xslo=fsolve(e.X,x0,arg)
xslo=root(e.X, x0, arg)

bonds=((71.5,110),(-50,50),(-5,20),(-0.5,0.5),(-0.5,0.5),(-0.5,0.5),(-30/180*math.pi,30/180*math.pi),(-30/180*math.pi,30/180*math.pi),(-20,20),(-20,20),(-20,20),(0,2))

#xsol=minimize(X,x0, bounds=bonds)

def postx(x):
    V=math.sqrt(np.dot(x[0:3],x[0:3]))
    alpha=math.atan(x[2]/x[0])/math.pi*180
    beta=math.asin(x[1]/V)/math.pi*180
    phi=x[6]/math.pi*180
    theta=x[7]/math.pi*180
    da=x[8]/math.pi*180
    de=x[9]/math.pi*180
    dr=x[10]/math.pi*180
    
    dx=x[-1]
    
    print("V = {0:0.2f} m/s, alpha = {1:0.1f}°, beta = {2:0.1f}°".format(V,alpha, beta))
    print("phi= {0:0.2f}°, theta={1:0.2f}°".format(phi, theta))
    print("da = {0:0.2f}°, de = {1:0.2f}°, dr = {2:0.2f}°".format(da,de,dr) )
    print(dx)
    
postx(xslo.x)