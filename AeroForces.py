# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 10:51:13 2017

Module for coef interpolation and forces calculation

@author: e.nguyen-van
"""
import numpy as np
import math
import ATRgeometry

def Interpol(V, A1, A2, v1, v2):
    # Function to interpol any kind of variables in function of the velocity V
    # input : 
    # V : current velocity
    # A1, A2 : lower matrix, higher matrix
    # v1, v2 : velocities corresponding to matrices
    a=(A2-A1)/(v2-v1)
    b=A1-a*v1
    Areturn=a*V+b
    return Areturn

def InterpolRho(V, rho, v):
    # function to interpol Rho, since altitude is function of Velocity
    if V<v[0] :
        return rho[0]
    
    elif V>v[-1]:
        return rho[-1]
    else:
        exitcondition=1
        length_v=len(v)-1
        i=0
        while exitcondition :
           
            if V==v[i]:
                rhoreturn=rho[i]
                exitcondition=0
            
            elif V>v[i] and V<v[i+1]:
                rhoreturn=Interpol(V, rho[i], rho[i+1], v[i], v[i+1])
                exitcondition=0 #exit
            
            else:
                i=i+1
                
            if i==length_v: #security to exit the while
                print("AeroForces : Error in interpolating rho, returning 0")
                rhoreturn=0
                exitcondition=0
    
    return rhoreturn

def CoefInterpol( V, A, v):
    # A, function takes a numpy array composed of all matrices A [A1; A2; ...], all array types!!
    # v, an array of corresponding velocity
    # V, the velocity at which to compute the coef
    # size of each matrix 6 row, m column
    row=6
    nmatrix=len(A[:,1])/6
    if nmatrix == float:
        #error
        print('Aero forces: general coef matrix not a multiple of 6')
        return 0
    
    if V<v[0] :
        #ill posed problem, the velocity is below smallest velocity for coef
        # use first matrix but print warning
        print("WARNING : velocity, V = {0:0.2f} is below first ref velocity for coef v = {1:0.2f}".format(V,v[0]))
        print("Continue with first matrix")
        return A[0:row,:]
    
    elif V>v[-1]:
        # same idea, the velocity is greater than max used to determine coef. Results in flight faster than cruise
        # use last matrix but print warning
        print("WARNING : velocity, V = {0:0.2f} is higher than last ref velocity for coef v = {1:0.2f}".format(V,v[-1]))
        print("Continue with last matrix")
        return A[-6:]
    
    elif V<0.2:
        # hard coded use the first matrix with flaps
        return A[0:row,:]
        
    else : #otherwise interpolate
        exitcondition=1
        length_v=len(v)-1
        i=1
        while exitcondition :
           
            if V==v[i]:
                Areturn=A[i*row:i*row+row]
                exitcondition=0
            
            elif V>v[i] and V<v[i+1]:
                Areturn=Interpol(V, A[i*row:i*row+row,:], A[(i+1)*row:(i+1)*row+row,:], v[i], v[i+1])
                exitcondition=0 #exit
            
            else:
                i=i+1
                
            if i==length_v+1: #security to exit the while
                print("AeroForces : Error in interpolation, returning 0")
                Areturn=0
                exitcondition=0

    return Areturn

        

""" Deprecated function do not use
def CalcForces(V, x, A, v, rho):
    # 
    # first interpol the coefficients matrix
    CoefMatrix=CoefInterpol(V, A, v)
    g=ATRgeometry.data(1)
    
    #Compute aero forces
    # here x must be of the form (alpha, beta, p, q, r, da, dr, de, dx) (last two punctualy used)
#    F=np.dot(CoefMatrix,x[0:7]) # commented form, modification to account for symmetric drag increase of side slip
    F=np.zeros((3))
    xsym=np.copy(x[0:7])
    xsym[1]=abs(xsym[1]) # make beta always positive since derivatives have already correct sign
    xsym[-4]=abs(xsym[-4]) # make ailerons deflection always positive for drag increase and lift decrease
    xsym[-3]=abs(xsym[-3]) # make rudder deflection always positive for drag increase and lift decrease
    F[0]=np.dot(CoefMatrix[0],xsym)
    F[1]=np.dot(CoefMatrix[1],x[0:7]) #side force
    F[2]=np.dot(CoefMatrix[2],xsym)
    M=np.dot(CoefMatrix[3:6,:],x[0:7]) 
    
    #Compute projection matrix and project
    alpha=x[0]
    beta=x[1]
    H=np.array([[math.cos(alpha)*math.sin(beta), -math.cos(alpha)*math.sin(beta), -math.sin(alpha)],[math.sin(beta), math.cos(beta), 0],[math.sin(alpha)*math.cos(beta), -math.sin(alpha)*math.sin(beta), math.cos(alpha)]])
    
    if V<=v[1] :
        Fbody=np.dot(H,[-F[0]-g.Cd0_fl,F[1],-F[2]-g.CL0_fl]) #project and add alpha=0 coefficients
        Moment=M+np.array([0,x[-2]*g.Cm_de+g.Cm0_fl,0])
    else:
        Fbody=np.dot(H,[-F[0]-g.Cd0,F[1],-F[2]-g.CL0]) #project and add alpha=0 coefficients
        Moment=M+np.array([0,x[-2]*g.Cm_de+g.Cm0,0])
    
    thrust=2*g.basethrust*InterpolRho(V, rho, v)/V*x[-1]
    
    Fbody=0.5*V**2.0*InterpolRho(V, rho, v)*g.S*Fbody+np.array([thrust,0,0])
    Moment=0.5*V**2.0*InterpolRho(V, rho, v)*g.S*g.b*Moment
    
    return np.append(Fbody, Moment)



def CalcForce_aeroframe(V, CoefMatrix, x, rho, g):
    "" Function to compute aerodynamic forces in the velocity frame (aero frame)
    for the original ATR configuration. The propulsion force 'thrust' is computed here
    ""
    
    #Compute aero forces
    # here x must be of the form (alpha, beta, p, q, r, da, dr, de, dx) (last two punctualy used)
   #    F=np.dot(CoefMatrix,x[0:7]) # commented form, modification to account for symmetric drag increase of side slip
    F=np.zeros((3))
    xsym=np.copy(x[0:7])
    xsym[1]=abs(xsym[1]) # make beta always positive since derivatives have already correct sign
    xsym[-4]=abs(xsym[-4]) # make ailerons deflection always positive for drag increase and lift decrease
    xsym[-3]=abs(xsym[-3]) # make rudder deflection always positive for drag increase and lift decrease
    F[0]=np.dot(CoefMatrix[0],xsym)
    F[1]=np.dot(CoefMatrix[1],x[0:7]) #side force
    F[2]=np.dot(CoefMatrix[2],xsym)
    M=np.dot(CoefMatrix[3:6,:],x[0:7]) 
    
    
    #No need to project
#    alpha=x[0]
#    beta=x[1]
#    H=np.array([[math.cos(alpha)*math.sin(beta), -math.cos(alpha)*math.sin(beta), -math.sin(alpha)],[math.sin(beta), math.cos(beta), 0],[math.sin(alpha)*math.cos(beta), -math.sin(alpha)*math.sin(beta), math.cos(alpha)]])
    
    if V<=71 :
        Fbody=[-F[0]-g.Cd0_fl,F[1],-F[2]-g.CL0_fl] # add alpha=0 coefficients
        Moment=M+np.array([0,x[-2]*g.Cm_de+g.Cm0_fl,0])
    else:
        Fbody=[-F[0]-g.Cd0,F[1],-F[2]-g.CL0] # add alpha=0 coefficients
        Moment=M+np.array([0,x[-2]*g.Cm_de+g.Cm0,0])
    
    thrust=2*g.basethrust*rho/V*x[-1]
    
    Fbody=0.5*V**2.0*rho*g.S*Fbody+np.array([thrust,0,0])
    Moment=0.5*V**2.0*rho*g.S*g.b*Moment
    
    return np.append(Fbody, Moment)
"""

def CalcForce_aeroframe_DEP(V, CoefMatrix, x, rho, g):
    """ Function to compute aerodynamic forces in the velocity frame (aero frame)
    for the DEP configuration. The propulsion force and moments are not computed here
    Since V is fixed, Coef Matrix must be calculated before
    Can handle DEP, with and without rudder, 2 or more engines
    """

    #Compute aero forces
    # here x must be of the form (alpha, beta, p, q, r, da, dr, de) (last one punctualy used)
    # set non dim for p,q,r
    nonDim=np.ones(7)
    nonDim[2]=g.b/(2*V)
    nonDim[3]=g.c/(2*V)
    nonDim[4]=g.b/(2*V)
    #    F=np.dot(CoefMatrix,x[0:7]) # commented form, modification to account for symmetric drag increase of side slip
    F=np.zeros((3))
    M=np.zeros((3))
    xsym=np.copy(x[0:-1])
    xsym[1]=abs(xsym[1]) # make beta always positive since derivatives have already correct sign for drag and lift only
    xsym[-3]=abs(xsym[-3]) # make ailerons deflection always positive for drag increase and lift decrease
    xsym[-1]=abs(xsym[-1]) # make rudder deflection always positive for drag increase and lift decrease
    F[0]=np.dot(CoefMatrix[0],xsym)
    F[1]=np.dot(CoefMatrix[1],x[0:-1]) #side force
    F[2]=np.dot(CoefMatrix[2],xsym)
    M=np.dot(CoefMatrix[3:6,:],x[0:-1])
#    print("Printing moment coeff")
#    print(M)

    
    #No need to project
#    alpha=x[0]
#    beta=x[1]
#    H=np.array([[math.cos(alpha)*math.sin(beta), -math.cos(alpha)*math.sin(beta), -math.sin(alpha)],[math.sin(beta), math.cos(beta), 0],[math.sin(alpha)*math.cos(beta), -math.sin(alpha)*math.sin(beta), math.cos(alpha)]])
    if V<=71 :
        Fbody=np.array([-F[0]-g.Cd0_fl,F[1],-F[2]-g.CL0_fl]) # add alpha=0 coefficients
        Moment=M+np.array([0,x[-1]*g.Cm_de+g.Cm0_fl,0])
    else:
        Fbody=np.array([-F[0]-g.Cd0,F[1],-F[2]-g.CL0]) # add alpha=0 coefficients
        Moment=M+np.array([0,x[-1]*g.Cm_de+g.Cm0,0])
        

    Fbody=0.5*V**2.0*rho*g.S*Fbody
    Moment=0.5*V**2.0*rho*g.S*g.b*Moment
    
    return np.append(Fbody, Moment)