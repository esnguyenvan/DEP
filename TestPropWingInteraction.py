# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 13:51:52 2018

@author: e.nguyen-van

testing the propeller-wing interaction
Composite model of Patterson and Jameson

suppose V/Vj=V/Vep and Vep determined from given formula (not Vinf+Vp)

"""

import numpy as np
import ATRgeometryAutomaticPropDia as ATR
import PattersonSurro as Surro
import ReadCoefAero
import AeroForces
import matplotlib.pyplot as plt
import sys
import math

N_eng=20
n_inop=0
Vtail=1.0

g=ATR.data(Vtail,N_eng, n_inop)
g.hangar={'aircraft':'ATR72', 'version':'DEPoriginal'} # new version : depATR72PropWingInterac
#print some info
print("Test Propeller Wing interaction")
print("Composite Patterson Jameson")

print("Configuration: \nEngines : {0:0.0f}\nProp diameter : {1:0.3f}m\n".format(g.N_eng, g.Dp))
print("Engines inoperational : {0:0.0f}".format(g.inop))
print("Engine position :")
strout=""
for i in range(len(g.PosiEng)):
    strout=strout+str(i+1)+" : "+"{:.3}".format(g.PosiEng[i])+"m, "
print(strout)

# initialize a throttle vector
dx=np.linspace(0,1,g.N_eng-g.inop)

if g.inop != 0:
    #complete vector if inop engines
    dx=np.append(dx,np.zeros(g.inop))

print("Throttle vector: ")
strout=""
for i in range(len(dx)):
    strout=strout+str(i+1)+" : "+"{:.2}".format(dx[i])+", "
print(strout)


#flight position
H_base = 0000 # in m the altitude
V_base = 65
beta_base = 0/180*math.pi
gamma = math.atan(3/100)#/180*math.pi # 3% slope gradient
R = 0 # in meters the turn radius
phimax = 5 # in degree the max bank angle authorized
alphamax=10 # in degree to adjust if it is below 71m/s
deltaRmax=30 # in degree

Velocities=(56.1,71,71,90,130)
rho_vec=(1.225,1.225,1.225,1.225*0.7812,1.225*0.4549)
#a_sound=(340.3,340.3,340.3,330.6,310.2)
Mach=[0.0,0.2,0.3,0.4,0.5]

#load coef files
if g.hangar['version']=='original' and g.nofin==True:
    print("WARNING : Using "+g.hangar['version']+" without rudder. Not recommended...")
if g.hangar['version']=='DEPoriginal' and g.inop!=0 and g.nofin==True:
    print("WARNING : Using "+g.hangar['version']+" with inoperative engine and without rudder. Not recommended...")
filename=['ATR72_SI_MTOW_Control_flap.stab', 'ATR72_SI_MTOW_Control_flapVmax.stab',
          'ATR72_SI_MTOW_Control_Vmin.stab','ATR72_SI_MTOW_Control_Manoeuver.stab',
          'ATR72_SI_MTOW_Control_Cruise.stab']
path='ATR72_SI_MTOW_Control_FinLess_DegenGeom_19_2_15h32/'#'D:/e.nguyen-van/Documents/DEP_control/StabilityMapATR72/ATR72_SI_MTOW_Control_FinLess_DegenGeom_19_2_15h32/'
filenameNoFin=[path+'_FinLess_Mach0.stab',path+'_FinLess_Mach02.stab',path+'_FinLess_Mach03.stab',
           path+'_FinLess_Mach04.stab',path+'_FinLess_Mach05.stab']
MatrixNoFin=ReadCoefAero.ReadCoef(filenameNoFin,col=(11,24,37,50,63,76,128),lign=(47,48,49,50,51,52))
MatrixFin=ReadCoefAero.ReadCoef(filename)
# Use coefficient of no fin ATR and add the coefficient for the rudder.
Matrix=np.copy(MatrixNoFin)
Matrix[:,-1]=MatrixFin[:,-1]
#
print("Adjusting Kf, Kh and VT size")
print("New Kf and Kh")
print(g.AdjustVT())
CoefMatrix=g.NicolosiCoef(Matrix[:,1:], Mach)
BaseAero=Matrix[:,1]

if g.hangar['version']=="depATR72PropWingInterac":
    sys.exit("Version not yet implemented")


#interpolated coefficient
atmospher=g.GetAtmo(H_base)
a_sound=atmospher[0]
rho_base=atmospher[1]
M_base=V_base/a_sound
Coef_base=AeroForces.CoefInterpol(M_base, CoefMatrix, Mach)

# check method at mach 0 with flaps and mach 0.3 without
#Compute thrust
Thrust = dx * 2*g.P_a / (g.N_eng * V_base)*g.prop_eff
Tc=Thrust/(0.5*rho_base*g.Sp*V_base**2)

#solve equation for Vp/V
Vp=np.zeros(len(Tc))
for i in range(len(Tc)):
    coef=[1,2*0.1,1,0,-Tc[i]**2]
    roots=np.roots(coef)
    #get the real positive root
    for j in range(len(roots)):
        if np.real(roots[j])>0:
            Vp[i]=np.real(roots[j])*V_base

Vjet = Vp+np.ones(len(Tc))*V_base
