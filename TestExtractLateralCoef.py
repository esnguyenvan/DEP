# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 17:12:26 2018

@author: e.nguyen-van
"""

import ReadCoefAero
import AeroForces
import numpy as np
import math
import scipy.linalg
import matplotlib.ticker as ticker
import matplotlib as mpl
import matplotlib.pyplot as plt
import ATRgeometry
import equation as e
from scipy.optimize import  minimize
import linecache
#import pdb


font = {'family':'sherif',
        'serif':'Times New Roman',
        'size'   : 22}

plt.rc('font', **font)

# --- Option de sortie console ---
np.set_printoptions(threshold=np.nan)

""" Algorithm set up """
Neng = 12
inop_eng = 0
g = ATRgeometry.data(1.0,Neng,inop_eng) # arg = Vtsize + options(Neng, inop_eng, vertical tail parameters)
# --- Test case and steady parameters
V_base = 90
beta_base = -5/180*math.pi
gamma = 0/180*math.pi

# --- dictionnary for type of aircraft studied. aircraft: ATR72, version : 'original', 'DEPoriginal', 'DEPsmallfin', 'DEPnofin'
g.hangar={'aircraft':'ATR72', 'version':'DEPsmallfin'}

# --- additional parameters (default edited during execution) ---
g.set_nofin(True) # =True means : no rudder used

# --- hard coded velocity, corresponding air density and Sound vel ---
# The velocities corresponds to key points of ATR flights and the ones at which
# VSPaero computes the stability matrices
Velocities=(56.1,71,71,90,130)
rho_vec=(1.225,1.225,1.225,1.225*0.7812,1.225*0.4549)
a_sound=(340.3,340.3,330.6,310.2)
Mach = np.array([0.0,0.2,0.3,0.5])


""" Algorithm start """

#--- List all .stab file from vsp aero and read the coeff ---- 
filename=['ATR72_SI_MTOW_Control_flap.stab', 'ATR72_SI_MTOW_Control_flapVmax.stab',
                  'ATR72_SI_MTOW_Control_Vmin.stab','ATR72_SI_MTOW_Control_Manoeuver.stab',
                  'ATR72_SI_MTOW_Control_Cruise.stab']
#filenameNoFin=['ATRfinless_flap.stab', 'ATRfinless_flapVmax.stab',
#          'ATRfinless_Vmin.stab','ATRfinless_Manoeuver.stab','ATRfinless_Cruise.stab']
path='D:/e.nguyen-van/Documents/DEP_control/StabilityMapATR72/ATR72_SI_MTOW_Control_FinLess_DegenGeom_19_2_15h32/'
filenameNoFin=[path+'_FinLess_Mach0.stab',path+'_FinLess_Mach02.stab',path+'_FinLess_Mach03.stab',
                   path+'_FinLess_Mach04.stab',path+'_FinLess_Mach05.stab']
MatrixNoFin=ReadCoefAero.ReadCoef(filenameNoFin,col=(11,24,37,50,63,76,128),lign=(47,48,49,50,51,52))
MatrixFin=ReadCoefAero.ReadCoef(filename)
Matrix=np.copy(MatrixNoFin)
Matrix[:,-1]=MatrixFin[:,-1]
print("Adjusting Kf, Kh and VT size")
print("New Kf and Kh")
print(g.AdjustVT())
CoefMatrix=g.NicolosiCoef(Matrix[:,1:], Mach)
BaseAero=Matrix[:,1]

cstAR=True
#CoefMatrix-MatrixFin[:,1:]

def ExtractLatCoef(Matrix, limit=0):
    #check input
    if limit==0:
        nmatrix=len(Mach)
    elif limit <= len(Mach):
        nmatrix=limit
    else :
        nmatrix=len(Mach)
        
    VarPosi=(1,2,4)
    EffPosi=(1,3,5)
    DispMat=np.zeros((nmatrix*3,3))
    for i in range(nmatrix):
        for k in range(len(EffPosi)):
            for l in range(len(VarPosi)):
                DispMat[k+i*3,l]=Matrix[EffPosi[k]+i*6,VarPosi[l]]
                
    return DispMat

def PlotMatrixCoef(What,M1,M2=np.zeros((2,1)),M3=np.zeros((2,1)),**keywords):
    # assum matrices are nx3
    #check What to plot
    if What==0:
        limit=0 # plot everything
    elif What<=len(Mach):
        limit=What
    else:
        limit=1
        print("WARNING : out of bound index Mach, using second limit to plot")
        
    #plot coefficients in graph
    xaxis=np.array((1,2,3))
    yLabelList=("$C_Y$","$C_l$","$C_n$")
    LabelList=("M1","M2","M3")
    MarkerStyle=("o","^","v","s","*")
    MarkerColor=[0.8,0.6,0.4,0.2,0]
    SupLim=[0.6,0.4,0.1]
    InfLim=[-0.5,-0.6,-0.3]
    xticks=[r'$\beta$',"p","r"]
    for i in range(3):
        plt.figure()
        plt.ylabel(yLabelList[i])
        plt.grid()
        if limit==0:
            size=len(M1[:,1])//3
            MarkerColor=abs(np.linspace(-1,0.1,size))
            for k in range(size):
                plt.xticks(xaxis,xticks)
                plt.scatter(xaxis, M1[i+k*3,:],color=str(MarkerColor[k]),label=LabelList[0],marker=MarkerStyle[0])#
        else:
            plt.scatter(xaxis, M1[(limit-1)*3+i,:],label=LabelList[0])
            
        if np.sum(M2)!=0:
            if limit==0:
                size=len(M1[:,1])//3
                MarkerColor=abs(np.linspace(-1,0,size))
                for k in range(size):
                    plt.scatter(xaxis, M2[i+k*3,:],color=MarkerColor[k],marker=MarkerStyle[1],label=LabelList[1])
            else:
                plt.scatter(xaxis, M2[(limit-1)*3+i,:],label=LabelList[1])
        
        if np.sum(M3)!=0:
            if limit==0:
                size=len(M1[:,1])//3
                MarkerColor=abs(np.linspace(-1,0,size))
                for k in range(size):
                    plt.scatter(xaxis, M3[i+k*3,:],color=MarkerColor[k],marker=MarkerStyle[2],label=LabelList[2])
            else:
                plt.scatter(xaxis, M3[(limit-1)*3+i,:],label=LabelList[2])
#        plt.legend()
        axe=plt.gca()
        plt.tight_layout()
        axe.set_ylim(InfLim[i],SupLim[i])
        plt.savefig(keywords['FigName'][i]+".pdf")
        
def CoefEvolution(aircraft, BaseMatrix, Sv, Machvec, CsteAR=True):
    # function to build CoefMatrix representing the evolution of VT size
    # Compute damping
    Cybeta=(0,0)
    Cnr=(2,2)
    if CsteAR==True:
        #Use constant Aspect Ratio fin changes
        aircraft.AdjustVT(Sv[0])
    else:
        #use constant span modification
        aircraft.AdjustVTcstSpan(Sv[0])
        
    CoefMatrix=g.NicolosiCoef(BaseMatrix, Machvec)
    ReducedMatrix=ExtractLatCoef(CoefMatrix,1) # extract the first matrix everytime
    Matrix=ReducedMatrix
    ZetaWn=np.array([-0.25*1.225*aircraft.S*65*(aircraft.b**2/aircraft.Iz*ReducedMatrix[Cnr]+ReducedMatrix[Cybeta]/aircraft.m)])
    for i in range(len(Sv)-1):
        if CsteAR==True:
            #Use constant Aspect Ratio fin changes
            aircraft.AdjustVT(Sv[i+1])
        else:
            #use constant span modification
            aircraft.AdjustVTcstSpan(Sv[i+1])
        
        CoefMatrix=g.NicolosiCoef(BaseMatrix, Machvec)
#        print(CoefMatrix)
        ReducedMatrix=ExtractLatCoef(CoefMatrix,1)
        Matrix=np.append(Matrix,ReducedMatrix,axis=0)
        ZetaWn=np.append(ZetaWn,[-0.25*1.225*aircraft.S*65*(aircraft.b**2/aircraft.Iz*ReducedMatrix[Cnr]+ReducedMatrix[Cybeta]/aircraft.m)])
        
    return Matrix, ZetaWn
            
#ExtractLatCoef((CoefMatrix-MatrixFin[:,1:])/MatrixFin[:,1:])
MNoFin=ExtractLatCoef(MatrixNoFin[:,1:])
MFin=ExtractLatCoef(MatrixFin[:,1:])
MVeDSC=ExtractLatCoef(CoefMatrix)

FinSizeRatio=np.linspace(0.3,1,10)
FinEvo, DampingEvo =CoefEvolution(g,Matrix[:,1:],FinSizeRatio,Mach, cstAR)
if cstAR:
    PlotMatrixCoef(0,FinEvo,FigName=('CyCstAR','ClCstAR','CnCstAR'))
else:
    PlotMatrixCoef(0,FinEvo, FigName=('CyCstSpan','ClCstSpan','CnCstSpan'))
#PlotMatrixCoef(0,MFin,MNoFin,MVeDSC)
#plt.figure
#PlotMatrixCoef(0,MVeDSC)

