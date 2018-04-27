# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 10:56:05 2018

@author: e.nguyen-van
"""

import numpy as np

def BetaSurro(a, Vjet, V, rho):
    """
    This function computes the beta, corrective term of Patterson propeller
    lift model.
    It implements the surrogate model present in the paper "High lift prop
    system for nasa's sceptor"
    Input variables:
        a : aircraft ATR class, with automatic limited propeller distribution
        deltax : vector of throttle level
        V : flight velocity
        rho : actual air density
    Outputs :
        beta : vector of beta value in the order of the deltax given
    """
    
    # first computes thrust
#    Thrust = deltax * 2*a.P_a / (a.N_eng * V)*a.prop_eff
    
    # use mu factor defined by jameson as V/Vj
#    invMu = ( 1 + (Thrust)/(0.5*rho*V**2*a.Sp) )**0.5
    invMu = Vjet/V
    
    # definition of surrogate coefficients
    C0 = np.array([0.378269, 0.748135, -0.179986, -0.056464, -0.146746, -0.015255])
    C1 = np.array([3.071020, -1.769885, 0.436595, 0.148643, -0.989332, 0.197940])
    C2 = np.array([-2.827730, 2.054064, -0.467410, -0.277325, 0.698981, -0.008226])
    C3 = np.array([0.997936, -0.916118, 0.199829, 0.157810, -0.143368, -0.057385])
    C4 = np.array([-0.127645, 0.135543, -0.028919, -0.026546, 0.010470, 0.012221])
    
    #Definition of surrogate vector
    Lratio = a.xp/a.c
    Rratio = a.Dp/(2*a.c)
    print(Lratio)
    print(Rratio)
    
    beta=np.zeros(len(Vjet))
    for i in range(len(Vjet)):
        X = np.array([1, Lratio, Lratio**2, Lratio*invMu[i], invMu[i], invMu[i]**2])
        
        # Compute the whole thing
        beta[i] = np.dot(C0,X) + np.dot(C1,X)*Rratio + np.dot(C2,X)*Rratio**2 + np.dot(C3,X)*Rratio**3 + np.dot(C4,X)*Rratio**4
    
    return beta