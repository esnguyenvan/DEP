# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 18:11:08 2017

Geometry class for the ATR72

A few constants are defined here:
    N_eng
    inop_eng
    Vertical Tail correction terms

The rest are terms relative to ATR geometry mass etc...

@author: e.nguyen-van
"""
import math
import sys
import numpy as np

class data:
    # all data from ATR72 go here
    hangar={}
    
    #shared data between class go here:
    # --- Geometry --- 
    S=61.0 #m^2
    b=27.05 #m
    c=2.32 #m
    lv=13 #m approximate, to verify
    zf=2 #3.16 #m, z position of the MAC of the fin, in reality a suitable height
    fswept = 35/180*math.pi
    ftaper=0.55
    fAR=1.57
    
    # --- Mass ---
    x_cg=12.4 # (m)
    m=21500 # Kg
    # Inertia terms are obtained from VSPaero from an homogeneous weight distribution
    Ix=60012 #Kg/m^2
    Iy=860916 #Kg/m^2
    Iz=894527 #Kg/m^2
    Ixz=25120 #Kg/m^2
    
    # --- Power ---
    P_a=2.051*10**6 #per engine
    P_var=P_a
    hp=0 # rotor term
    prop_eff=0.7
    
    
    # ---Unique coeff ---
    basethrust=15000
    Cm_de=-0.035*180/math.pi # per rad
    Cn_beta_f=0#-0.2137
#    Cn_beta_f=-0.15
    
    # alpha=0 coeff
    # with flaps down 30°
    Cd0_fl=0.08128
    CL0_fl=1.29604
    Cm0_fl=0.11806
    
    #without flaps
    Cd0=0.03383
    CL0=0.50264
    Cm0=0.06154
    
    #SmallFin coefficients determined using VeDSC method
    #Edited in 
    Cy_beta = 0
    Cy_p = 0
    Cy_r = 0
    #
    Cl_beta = 0
    Cl_p = 0
    Cl_r = 0
    #
    Cn_beta= 0
    Cn_p = 0
    Cn_n = 0
    # Additional coefficients
    taudr = 0.8 # see nicolosi paper and dela-vecchia thesis
    
    #unique data go here:
    def CalcKf(self, bv,r):
        return 1.4685*(bv/(2*r))**(-0.143)
    
    def CalcKw(self, zw,rf):
        return -0.0131*(zw/rf)**2-0.0459*(zw/rf)+1.0026
        
    def CalcKh(self, zh,bvl, Sh, Sv):
        x=zh/bvl
        Khp=0.906*x**2-0.8403*x+1.1475
        Khs=math.exp(0.0797*math.log(Sh/Sv)-0.0079)
        return 1+Khs*(Khp-1)
    
    def CalcKdr(self, Kf,Av):
        Kdr=(1+((Kf-1)/2.2))*(1.33-0.09*Av) # for T-tail formula
        return Kdr
    
    def set_nofin(self, boolean):
        if type(boolean)==bool:
               self.nofin=boolean # flag to use or not rudder
        else:
            sys.exit("Error type of 'nofin' isn't a boolean. Stopping")
            
    def loadAtmo(self):
        filename='si2py.txt'#'D:/e.nguyen-van/Documents/DEP_control/StabilityMapATR72/si2py.txt'
        sep='\t'
        file=open(filename,'r')
        vecname=file.readline()
        index=0
        VariableList=[]
        condition=True
        while condition:
            VariableList.append(vecname[index:vecname.index(sep,index)])
            if VariableList[-1]=='kvisc':
                condition=False
            index=vecname.index(sep,index)+1
            
        units=file.readline() # skip units
        data=[]
        VariableDic={}
        for j in range(len(VariableList)):
            exec("VariableDic['"+VariableList[j]+"'] = []") #initialize my variables
            
        for line in file:
            mylist=[]
            element=""
            for k in range(len(line)):
                if line[k]!='\t' and line[k]!='\n':
                    element=element+line[k]
                else:
                    mylist.append(element)
                    element=""
            data.append(mylist)
        file.close()
        
        for i in range(len(data)):
            for k in range(len(data[0])-1):
                exec("VariableDic['"+VariableList[k]+"'].append({})".format(float(data[i][k])))
            
        return VariableDic

    
    def __init__(self, VTsize, N_eng=2, inop_eng=0, bv=4.42, r=0.6, zw=1.8, rf=1.3, zh=3.71, bvl=5.91, Sh=11.13, Sv=12.5):
        self.VTsize=VTsize
        self.N_eng=N_eng # number of engines
        self.inop=inop_eng # number of inoperative engines
        
        if N_eng % 2 !=0:
            sys.exit("Number of engine is not even")
        
        if N_eng>2:
            #Adjust engine position
            self.step_y=self.b/self.N_eng
            self.PosiEng=np.append(np.arange(-self.b/2,0,self.step_y),np.arange(self.step_y,(self.b/2+self.step_y),self.step_y))
        else:
            self.step_y=8.1/2.0
            self.PosiEng=np.array([-self.step_y,self.step_y])
                
        self.Sv=Sv
        self.SvBase=Sv
        self.bv=bv
        self.r=r
        self.Av=bv**2/Sv
        self.Sh=Sh
        self.zw=zw
        self.rf=rf
        
        #Nicolosi csts
        self.Kf=self.CalcKf(bv,r)
        self.Kw=self.CalcKw(zw,rf)
        self.Kh=self.CalcKh(zh,bvl,Sh,Sv)
        self.Kdr=self.CalcKdr(self.Kf,self.Av)
        self.taudr=4/Sv*0.9
        #Atmosphere
        self.AtmoDic=self.loadAtmo()

    # functions go here:
    def printVTsize(self):
        print('Asked vertical tail size of current class:')
        print(self.VTsize)
    
    def NicolosiCoef(self, MCoef, Mach):
        # function to compute Cy, Cl and Cn derivatives using VeDSC methods
        # replaces the coefficients in the matrix Mcoef by those computed by VeDSC
        # to call only once, it computes for every velocities/Mach number. Then linear interpol
        # Cy_beta, Cy_p = 0, Cy_r = 0, Cl_beta = 0, Cl_r = 0, Cn_beta= 0, Cn_p = 0, Cn_n = 0
        
        MVeDSC=np.copy(MCoef) # create a copy of the coefficients matrix
        
        K=self.Kf*self.Kh*self.Kw
        for i in range(len(Mach)):
            Av = self.bv**2/self.Sv
#            print(Av)
            cla=2*math.pi # local lift curve slope coefficient of fin (perpendicular to leading edge) per rad
            eta=cla/(2*math.pi)
            av = -(Av * cla * math.cos(self.fswept) * 1/math.sqrt(1-Mach[i]**2*(math.cos(self.fswept))**2))/(Av*math.sqrt(1+4*eta**2/(Av/math.cos(self.fswept))**2)+2*eta*math.cos(self.fswept))# mach number formula
#            print(av)
            VeDSC_Coef=np.array([[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0],[0.0,0.0,0.0,0.0]])
            VeDSC_Coef[0,0]=K*av*self.Sv/self.S
            VeDSC_Coef[0,1]=K*av*self.Sv/self.S*self.zf/self.b*2
            VeDSC_Coef[0,2]=-K*av*self.Sv/self.S*2*self.lv/self.b #apparently x position of fin doesn't change
            VeDSC_Coef[0,3]=-self.Kdr*av*self.taudr*self.Sv/self.S
            VeDSC_Coef[1,0]= K*av*self.Sv/self.S*2*self.zf/self.b
            VeDSC_Coef[1,1]=0
            VeDSC_Coef[1,2] = -K*av*self.Sv/self.S*self.zf/self.b*self.lv/self.b*2.0
            VeDSC_Coef[1,3] = -self.Kdr*av*self.taudr*self.Sv/self.S*2*self.zf/self.b
            VeDSC_Coef[2,0] = -K*av*self.lv/self.b*self.Sv/self.S
            VeDSC_Coef[2,1] = -K*av*self.lv/self.b*self.Sv/self.S*self.zf/self.b*2.0
            VeDSC_Coef[2,2] = K*av*self.lv/self.b*self.Sv/self.S*self.lv/self.b*2.0
            VeDSC_Coef[2,3] = self.Kdr*av*self.taudr*self.lv/self.b*self.Sv/self.S
            # Coefficients are computed now access the right matrix and replace them
            VarPosi=(1,2,4,6)
            EffPosi=(1,3,5)
            NumEff = 6 # number of force equations
            
            for kk in range(len(EffPosi)):
                #Manually change rudder coefficient by simple proportionality
                #MVeDSC[EffPosi[kk]+i*NumEff,-1]=MVeDSC[EffPosi[kk]+i*NumEff,-1]*self.VTsize
                # Replace rudder coefficient
                MVeDSC[EffPosi[kk]+i*NumEff,-1]=VeDSC_Coef[kk,-1]/180.*math.pi
                # Now coefficients are from the finless ATR. Add new coefficients to the matrix
                for jj in range(len(VarPosi)):
                    if VeDSC_Coef[kk,jj]!=0:
                        MVeDSC[EffPosi[kk]+i*NumEff,VarPosi[jj]]=MVeDSC[EffPosi[kk]+i*NumEff,VarPosi[jj]]+VeDSC_Coef[kk,jj]
            print(VeDSC_Coef)
        return MVeDSC
        
    
    def AdjustVT(self, VTNewSize=1):
        # Adjust VT geometry parameters based on VTsize
        # Change only Kf, Kw and Kh, does not modify internal goem param
        
        # local variables
        Sv=self.SvBase
        bv=self.bv
        r=self.r
        Av=self.Av
        
        if VTNewSize!=1:
            Sv=Sv*VTNewSize
        else:
            Sv=Sv*self.VTsize
        
        self.bv=(Av*Sv)**0.5
        self.zh=r+0.77*bv
        self.bvl=bv+r
        
        # Compute new coef if necessary
        self.Sv=Sv
        self.Kf=self.CalcKf(self.bv,r)
        self.Kh=self.CalcKh(self.zh,self.bvl,self.Sh,Sv)
        
        return np.array([self.Kf, self.Kh])
    
    def AdjustVTcstSpan(self, VTNewSize=1):
        # Adjust VT geometry parameters based on VTsize
        # Based on Sv and cste span (increase Av)
        # Change only Kf, Kw and Kh
        
        # local variables
        Sv=self.SvBase
        bv=self.bv
        
        if VTNewSize!=1:
            Sv=Sv*VTNewSize
        else:
            Sv=Sv*self.VTsize
        
        self.Av=bv**2/Sv
        
        # Compute new coef if necessary
        self.Sv=Sv
        self.Kh=self.CalcKh(self.zh,self.bvl,self.Sh,Sv)
        
        return np.array([self.Kh])
    
    def GetAtmo(self,h=0):
        """
        Using the atmosphere model loaded before, it outputs [a_sound, rho] at
        the desired h=altitude. It doesn't perform interpolation.
        """
        Condition=h/500 is int 
        if Condition:
            Indice=h//500+1
            
        else:
            if (h/500)<(h//500+500/2):
                Indice=h//500+1
            else:
                Indice=h//500+2
        
        results=np.array([self.AtmoDic['a'][Indice],self.AtmoDic['dens'][Indice]])
        return results