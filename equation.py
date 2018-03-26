# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 17:01:52 2017

Module defining non linear flight equations

@author: e.nguyen-van
"""

import numpy as np
import math
import AeroForces
from numpy.linalg import inv
# functions for standard earth equation



def X(x, CoefMatrix,Velocities,rho,g):
    # x : state [u,v,w,p,q,r,phi,theta,dfl,da,de,dr,dx]
    # arg : tuple of additional elements, respect orger below:
#    CoefMatrix=arg[0]
#    Velocities = arg[1]
#    rho = arg[2]
#    g = arg[3]
    
    con=np.array([80,0/180*math.pi,0,0,0,0,0,0,0,0,0,0])
    if x[7]<-30/180*math.pi:
        return np.ones(12)*1000-np.ones(12)*x[7]*1000
    elif x[7]>30/180*math.pi:
        return np.ones(12)*1000+np.ones(12)*x[7]*1000
    
    if x[-1]<0:
        return np.ones(12)*1000-np.ones(12)*x[-1]*10000
    elif x[-1]>1:
        return np.ones(12)*1000+np.ones(12)*x[-1]*1000
    
    # need to do some operation before computing the forces
    V=math.sqrt(np.dot(x[0:3],x[0:3]))
    sub_vect=np.array([math.atan(x[2]/x[0]), math.asin(x[1]/V)])
    sub_vect=np.append(sub_vect,x[3:6])
    #(alpha, beta, p, q, r, da, dr, de, dx)
    sub_vect=np.append(sub_vect,[x[-4],x[-2],x[-3],x[-1]])

    
    F=AeroForces.CalcForces(V, sub_vect, CoefMatrix, Velocities, rho)
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
    
    return xreturn

def Constraints_DEP(x, fix, CoefMatrix, rho, g):
    """function defining constraints for power minimization
    inputs:
        -x =[alpha, p, q, r, phi, theta, delta_a, delta_e, delta_r, delta_i]
        x is the state to determine
        length of x except the propulsion levels is 8
        -fix = [V, beta, gamma, omega]
        fix is the vector of parameters whom are fixed by the user
        
    """
    
    # First thing to do is to determine the number of engines on each semi wing
    n_eng=int(g.N_eng/2)

    # --- Now prepare variables for equations ---
    V=fix[0]
    alpha=x[0]
    beta=fix[1]
    gamma=fix[2]
    omega=fix[-1]
    p=x[1]
    q=x[2]
    r=x[3]
    phi=x[4]
    theta=x[5]
    I=np.array([ [g.Ix, 0, -g.Ixz], [0, g.Iy, 0], [-g.Ixz, 0, g.Iz] ])
    
    # --- Compute aerodynamic forces ---
    #here subvector  must be : (alpha, beta, p, q, r, da, dr, de, dx)
    sub_vect=np.array([alpha,beta,p,q,r])
    if g.nofin==False:
        sub_vect=np.append(sub_vect,[x[6],x[8],x[7]]) # rudder is allowed
    else:
        sub_vect=np.append(sub_vect,[x[6],0,x[7]]) # no fin allowed, default case
  
    F=AeroForces.CalcForce_aeroframe_DEP(V, np.copy(CoefMatrix), np.copy(sub_vect), rho, g)
#    print("Forces :")
#    print(F)
    
    #Now compute thrust and moments
    start=8
    if g.nofin==False:
        start=9
    Fx=np.sum(x[start:])*2*g.P_var/(float(g.N_eng)*V)*g.prop_eff
    
    M_y=0 # initialization
    xeng=np.copy(x[start:])
    M_y=-np.dot(xeng,g.PosiEng)
    # for i in range(n_eng):
    #     M_y=M_y+x[start+i]*g.step_y*(n_eng-i)
    # for i in range(n_eng):
    #     M_y=M_y-x[-i-1]*g.step_y*(n_eng-i)
        
    M_y=M_y*2*g.P_var/(float(g.N_eng)*V)*g.prop_eff
#    print(M_y)

#     Now sum up the constraints:
    sinbank=np.sin(theta)*np.cos(alpha)*np.sin(beta) + np.cos(beta)*np.cos(theta)*np.sin(phi)-np.sin(alpha)*np.sin(beta)*np.cos(theta)*np.cos(phi)
    cosbank=np.sin(theta)*np.sin(alpha)+np.cos(beta)*np.cos(theta)*np.cos(phi) 
    
    A=np.zeros(10+g.inop) 
    """
    A0 = x
    A1 = y
    A2 = z
    A3 = l
    A4 = m
    A5 = n
    A6 = phi
    A7 = Omega
    A8 = gamma
    A9 = theta
    """
    A[0]=-9.81*np.sin(gamma)+F[0]/g.m+Fx*np.cos(alpha)*np.cos(beta)/g.m
    A[1]=(p*np.sin(alpha) - r*np.cos(alpha))+g.m*9.81*sinbank/(g.m*V) + F[1]/(g.m*V)-Fx*np.cos(alpha)*np.sin(beta)/(g.m*V)
    A[2]=-(np.sin(beta)*(p*np.cos(alpha)+r*np.sin(alpha))-q*np.cos(beta))/np.cos(beta)+ 9.81*cosbank/(V*np.cos(beta)) + F[2]/(g.m*V*np.cos(beta))-Fx*np.sin(alpha)/(g.m*V*np.cos(beta))
    A[3:6]=np.dot(inv(I), np.array([0,0,M_y])+F[3:6]-np.cross(np.array([p,q,r]),np.dot(I,np.array([p,q,r]))))
    A[6]=p+q*np.sin(phi)*np.tan(theta)+r*np.cos(phi)*np.tan(theta)
    A[7]=-omega + (q*np.sin(phi)+r*np.cos(phi))/np.cos(theta)
    A[8]=-np.sin(gamma)+np.cos(alpha)*np.cos(beta)*np.sin(theta)-np.sin(beta)*np.sin(phi)*np.cos(theta)-np.sin(alpha)*np.cos(beta)*np.cos(phi)*np.cos(theta)
    A[9]=q*math.cos(phi) -r*math.sin(phi)
    for i in range(g.inop):
        A[-1-i]=x[-1-i]
    
    if g.hangar['version']=='original':
        #no DEP with original twin or N engines
        D=np.copy(A)
        for i in range(g.N_eng-g.inop-1):
            AAd=x[start]-x[start+i+1]
            D=np.append(D,[AAd])
        return D
    else:
        return A

def fobjective(x, fix, rho, g):
    
#    Power=np.sum(x[-g.N_eng:])*2*g.P_var/float(g.N_eng)*rho/1.225/1000000
    Power=np.sum(x[-g.N_eng:])*2*g.P_var/float(g.N_eng)/1000000
    
    return Power

def Jac_DEP(x, fix, CoefMatrix, rho, g, h):
    # function to compute the jacobian at a steady state
    # the function is hard coded inside
    # inputs :
    #       -x : steady state vector
    #       -fixtuple : tuple of (fixed param, function param)
    #       -h : step to compute derivative
    
    nfx=7 # number of equations for flight analysis (V, beta, alpha, p, q, r, phi )
    
    step_vec=x*h
    
    for i in range(len(step_vec)):
        # check for zeros
        if step_vec[i]<1e-4:
            step_vec[i]=h
     
#    fx=Constraints_DEP(x, *fixtuple)
    
    dx=np.zeros((nfx,len(x)+2))
    fixtuple=(fix, CoefMatrix, rho, g)
    # compute derivative using centered difference
    #First velocity included in fix
    fix_plus=fix+np.append([fix[0]*h/2.0],np.zeros((len(fix)-1)))
    fix_minus=fix-np.append([fix[0]*h/2.0],np.zeros((len(fix)-1)))
    
    tuple_plus=(fix_plus, CoefMatrix, rho, g)
    tuple_minus=(fix_minus, CoefMatrix, rho, g)
    
    diff=(Constraints_DEP(x,*tuple_plus)-Constraints_DEP(x,*tuple_minus))/fix[0]*h
    dx[:,0]=diff[0:nfx]
    
    #Same for beta
    beta_step=np.zeros((len(fix)))
    beta_step[1]=h/2
    fix_plus=fix+beta_step
    fix_minus=fix-beta_step
    
    tuple_plus=(fix_plus, CoefMatrix, rho, g)
    tuple_minus=(fix_minus, CoefMatrix, rho, g)
    
    diff=(Constraints_DEP(x,*tuple_plus)-Constraints_DEP(x,*tuple_minus))/beta_step[1]
    dx[:,1]=diff[0:nfx]
    
    for j in range(len(x)):
        activex=np.zeros((len(x)))
        activex[j]=1
        dfx=(Constraints_DEP(x+activex*step_vec/2,*fixtuple)-Constraints_DEP(x-activex*step_vec/2,*fixtuple))/np.dot(activex,step_vec)
        dx[:,j+2]=dfx[0:nfx]
    
    # optionally decouple matrix
    return dx