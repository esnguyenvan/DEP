# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 15:25:07 2018

@author: e.nguyen-van
"""

import numpy as np
import math
import sys

class MapVar:
    Names={'Vel':'Velocity', 'Beta':'Side Slipe angle', 'Gamma':'Climb angle', 'Omega':'Turn rate'}
    Units={'Vel':'m/s','Beta':'\xb0','Gamma':'\xb0','Omega':'\xb0/s'}
    fixposition={'Vel':0, 'Beta':1, 'Gamma':2, 'Omega':3}
    DimUnits="rad"
    
    def __init__(self,Var1,Var2,Min1,Max1,Step1,Min2,Max2,Step2):
        self.FirstVar=Var1
        self.SecondVar=Var2
        self.Min1=Min1
        self.Max1=Max1
        self.Step1=Step1
        self.Min2=Min2
        self.Max2=Max2
        self.Step2=Step2
        self.FirstDim = np.arange(Min1,Max1,Step1)
        self.SecondDim = np.arange(Min2,Max2,Step2)
        
        #default behaviour transform to rad
        if Var1 == "Beta" or Var1 == "Gamma":
            self.FirstDim=self.FirstDim/180*math.pi
        
        if Var2 == "Beta" or Var2 == "Gamma":
            self.SecondDim=self.SecondDim/180*math.pi

    def ReArrangefix(self,i,j, Vel,Beta,gamma,omega):
        fix=np.array([Vel,Beta,gamma,omega])
        
        fix[self.fixposition[self.FirstVar]]=self.FirstDim[i]
        fix[self.fixposition[self.SecondVar]]=self.SecondDim[j]
        
        return fix
    
    def getUnit(self,a):
        if a=='x':
            return self.Units[self.FirstVar]
        elif a=='y':
            return self.Units[self.SecondVar]
        else:
            sys.exit("Error, only 'x' or 'y' entry correct for units")
        
    def getXlocator(self):
        delta=abs(self.Step1/2)
        Locators=np.copy(self.FirstDim[0:-1])
        if self.FirstVar=="Beta" or self.FirstVar=="Gamma":
            if self.DimUnits=="rad":
                Locators=Locators/math.pi*180
        return Locators+delta
        
    def getYlocator(self):
        delta=abs(self.Step2/2)
        Locators=np.copy(self.SecondDim[0:-1])
        if self.SecondVar=="Beta" or self.SecondVar=="Gamma":
            if self.DimUnits=="rad":
                Locators=Locators/math.pi*180

        return Locators+delta
    
    def getName(self,a):
        if a=='x':
            return self.Names[self.FirstVar]
        elif a=='y':
            return self.Names[self.SecondVar]
        else:
            sys.exit("Error, only 'x' or 'y' entry correct for variable name")
            
    def Dim2rad(self):
        if self.FirstVar == "Beta" or self.FirstVar == "Gamma" or self.FirstVar == "Omega":
            self.FirstDim=self.FirstDim/180*math.pi
        
        if self.SecondVar == "Beta" or self.SecondVar == "Gamma" or self.SecondVar == "Omega":
            self.SecondDim=self.SecondDim/180*math.pi
        
        self.DimUnits="rad"
    
    def Dim2deg(self):
        if self.FirstVar == "Beta" or self.FirstVar == "Gamma" or self.FirstVar == "Omega":
            self.FirstDim=self.FirstDim/math.pi*180
        
        if self.SecondVar == "Beta" or self.SecondVar == "Gamma" or self.SecondVar == "Omega":
            self.SecondDim=self.SecondDim/math.pi*180
            
        self.DimUnits="deg"