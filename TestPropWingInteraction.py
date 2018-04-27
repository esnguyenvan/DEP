# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 13:51:52 2018

@author: e.nguyen-van

testing the propeller-wing interaction
Composite model of Patterson and Jameson

suppose V/Vj=V/Vep and Vep determined from given formula (not Vinf+Vp)

"""

import numpy as np
#import ATRgeometryAutomaticPropDia as ATR
import PattersonSurro as Surro
import ReadCoefAero
import AeroForces
import matplotlib.pyplot as plt
import sys
import math
import TrialClread as ClRead
import matplotlib.pyplot as plt
import RectangularWingTestData as Plane

N_eng=16
n_inop=0
Vtail=1.0

g=Plane.data(Vtail,N_eng, n_inop)
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
#dx=np.linspace(0,1,g.N_eng-g.inop)
dx = np.ones(N_eng-g.inop)*0.5

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
V_base = 31.4
beta_base = 0/180*math.pi
# define angle of attack
alpha=5.5/180*np.pi #otherwise it is x[0]
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
#Thrust = dx * 2*g.P_a / (g.N_eng * V_base)*g.prop_eff
Thrust = np.ones(len(dx))*750 * g.Sp
Tc=Thrust/(2*rho_base*g.Sp*V_base**2)

#solve equation for Vp/V
Vp=np.zeros(len(Tc))
for i in range(len(Tc)):
    coef=[1,2*np.cos(0),1,0,-Tc[i]**2]
    roots=np.roots(coef)
    #get the real positive root
    for j in range(len(roots)):
        if np.real(roots[j])>0:
            Vp[i]=np.real(roots[j])*V_base
#check thrust
T=2*rho_base*g.Sp*Vp[0]*( (V_base*np.cos(0) + Vp[0])**2 + (V_base*np.sin(0))**2)**0.5

#can compute jet velocity
Vjet = Vp+np.ones(len(Tc))*V_base

# Can compute the surrogate for beta
BetaVec=Surro.BetaSurro(g,Vjet,V_base, rho_base)

# Need to import Cl distribution vsp aero files.
filename='Validationfischer_DegenGeom'#'ATR72_SI_DEP_DegenGeom'
LocalCl = ClRead.ReadSectionCl(filename)

#now for each prop, average the local lift coef
#SectionCl=np.zeros((len(Tc),2))
#for i in range(len(Tc)):
#    
#    #compute the bounds defined by prop position and diameter
#    #suppose wing symetrical
#    startprop = abs(g.PosiEng[i]) - g.Dp/2 # from fuselage
#    stopprop = abs(g.PosiEng[i]) + g.Dp/2 # toward tip
#    
#    indstart = 0
#    indstop = 0
#    for a in range(len(LocalCl[:,1])):
#        if LocalCl[a,0] > startprop:
#            indstart=a-1
#            break
#    for a in range(len(LocalCl[:,1])):
#        if LocalCl[a,0] > stopprop:
#            indstop=a-1
#            break
#
#    SubCl = LocalCl[indstart:indstop,:]
#    SectionCl[i,0] = np.sum(SubCl[:,1]*SubCl[:,-1])/np.sum(SubCl[:,1])
#    SectionCl[i,1] = np.sum(SubCl[:,1])

# define angle of attack
alpha=5.5/180*np.pi #otherwise it is x[0]
# if flaps
alpha_fl = alpha + g.FlEff

# total alpha
alpha_t = alpha - g.alpha_0
alpha_fl_t = alpha_fl - g.alpha_0

# Lift multipier Lm :
V_base_vec = np.ones(len(Tc))*V_base
#Without flaps
Lm = ( 1 - BetaVec*Vp*np.sin(g.ip)/(V_base*np.sin(alpha_t))) * (V_base_vec**2 + 2*V_base_vec*Vp*BetaVec*np.cos(alpha_t + g.ip) + (BetaVec*Vp)**2)**0.5/V_base-1

#with flaps
#determine the number of propellers blowing on flaps
Fltip=g.FusWidth/2+g.bflap/2
PropInFlap=0
for i in range(len(Tc)):
    if g.PosiEng[i]+g.Dp/2 > - Fltip and g.PosiEng[i]-g.Dp/2 > - Fltip:
        PropInFlap=i
        break
    if i==len(Tc):
        sys.exit('Could not determine the number of propeller in flap')
        
# means that N_prop/2-PropInFlap are blowing on flaps
NpBlowFlap = (g.N_eng/2-PropInFlap)*2

#lift multiplier
LmFl = np.zeros(len(Tc))
for i in range(len(Tc)):
    if i >= PropInFlap and i <= PropInFlap+NpBlowFlap:
        LmFl[i] = ( 1 - BetaVec[i]*Vp[i]*np.sin(g.ip)/(V_base*np.sin(alpha_fl_t))) * (V_base_vec[i]**2 + 2*V_base_vec[i]*Vp[i]*BetaVec[i]*np.cos(alpha_fl_t + g.ip) + (BetaVec[i]*Vp[i])**2)**0.5/V_base-1
    else:
        LmFl[i] = ( 1 - BetaVec[i]*Vp[i]*np.sin(g.ip)/(V_base*np.sin(alpha_t))) * (V_base_vec[i]**2 + 2*V_base_vec[i]*Vp[i]*BetaVec[i]*np.cos(alpha_t + g.ip) + (BetaVec[i]*Vp[i])**2)**0.5/V_base-1

##Sum lift multiplier over wing
#Tlm = np.sum(Lm*(g.c*g.Dp)/g.S)
#TlmFl = np.sum(LmFl*(g.c*g.Dp)/g.S)
#
##Local lift coefficient multiplied
#ClLm = np.sum(SectionCl[:,0] * (Lm+1) * SectionCl[:,1]) / g.S
#ClLmFl = np.sum(SectionCl[:,0] * (LmFl+1) * SectionCl[:,1]) / g.S

# Compute total CL with prop
tempCL=0
BlownCl = np.copy(LocalCl)
for i in range(int(len(LocalCl))):
    for a in range(len(Tc)):
        #general case
        if (g.PosiEng[a] - g.Dp/2) < LocalCl[i,0]  < (g.PosiEng[a] + g.Dp/2):
            LiftMulti = Lm[a]
            break
        else:
            LiftMulti = 0

    tempCL=tempCL+LocalCl[i,-1]*LocalCl[i,1]*(LiftMulti+1) #Cl*Slocal
    BlownCl[i,-1] = LocalCl[i,-1]*(LiftMulti+1)
#    print(LocalCl[i,0])
    
CL=tempCL/g.S
print("CL : {0:0.5f}, deltaCL : {1:0.1f}%".format(CL,(CL-0.6054)/0.6054*100))
plt.figure()
plt.plot(BlownCl[:,0],BlownCl[:,-1], color="0.45")
plt.ylabel('Local lift coefficient CL')
plt.xlabel('Span position (m)')
plt.title('Local lift coefficient blown flapped wing')