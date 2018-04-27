# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 11:14:12 2018

File that read 

@author: e.nguyen-van
"""

import linecache
import numpy as np

def ReadSectionCl(filename):
    # Eats the filename without extension
    # Throw a matrix with [YLEposition, Area, Chord, Cl]

#    filename = 'ATR72_SI_DEP_DegenGeom'
    extension = '.fem'
    
    Keywords =('Wing Surface: ', 'SpanStations: ', 'Wing')
    Wingplatform = ('1','2')
    
    
    f=open(filename+extension,'r')
    
    CurrentLine = 1
    nwing = 0 # wing count
    MySec = np.array([]).reshape(0,4)
    
    
    # readcoef function
    def ReadCoef(Nsec, line, filename):
        CLcol=[232,240] # colum where is the Cl
        YLEcol=[20,30] # Y leading edge
        ChCol = [122,130] # Chord col
        Acol = [112,120] # Area
        Dlign = 4 # offset
        Start = Dlign + line
        SecCl = np.zeros((Nsec*1,4)) #will contain YLEposition, Area, Cl
        
        for i in range(Nsec):
            CurrLine = linecache.getline(filename, Start + i)
            SecCl[i,0] = float(CurrLine[YLEcol[0]:YLEcol[1]])
            SecCl[i,1] = float(CurrLine[Acol[0]:Acol[1]])
            SecCl[i,2] = float(CurrLine[ChCol[0]:ChCol[1]])
            SecCl[i,3] = float(CurrLine[CLcol[0]:CLcol[1]])
        
        return SecCl
    
    for line in f:
        # find 'Wing Surface : x'
        if Keywords[0]+Wingplatform[nwing] in line:
            #find the number of section on this wing
            Next = linecache.getline(filename+extension,CurrentLine+1)
            NSection = int(Next[len(Keywords[1]):])
            # call the readcoef function and stack the data in MySec
            MySec = np.vstack([MySec, ReadCoef(NSection,CurrentLine,filename+extension)])
            nwing = nwing+1
            CurrentLine =CurrentLine + 1
        else:
            CurrentLine = CurrentLine+1
        
        if nwing >= 2:
            break
            
    f.close()
    
    return MySec