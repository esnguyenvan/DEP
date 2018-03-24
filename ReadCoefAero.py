# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 14:19:40 2017

Trial to read vsp aero text file and put them into a class

@author: e.nguyen-van
"""
import linecache
import numpy as np


def readstabfile(filename,col, lign):
    file = open(filename, 'r')
    
    # col are formated such that CL=with respect to [alpha, beta, p,q,r, Mach, U, dfl, da, dr]
    # We get rid of Mach, U et dfl
#    col=(11,24,37,50,63,76,128,141) # includes base aero, doesn't take flap effects
#    lign=(49,50,51,52,53,54)
    byte=8
    
    CL=[]
    CD=[]
    CY=[]
    Cl=[]
    Cm=[]
    Cn=[]
    
    # ----- Lift coefficients -------
    Coefstr=linecache.getline(filename, lign[0])
    
    for cursor in range(len(col)):
        if Coefstr[col[cursor]-2] != ' ':
            CL.append(float(Coefstr[col[cursor]-2:col[cursor]+byte]))
        else:        
            CL.append(float(Coefstr[col[cursor]-1:col[cursor]+byte]))
    
    del Coefstr
    
    # ----- drag coefficients -------
    Coefstr=linecache.getline(filename, lign[1])
    for cursor in range(len(col)):
        if Coefstr[col[cursor]-2] != ' ':
            CD.append(float(Coefstr[col[cursor]-2:col[cursor]+byte]))
        else:        
            CD.append(float(Coefstr[col[cursor]-1:col[cursor]+byte]))
    
    del Coefstr
    
    # ----- side force coefficients -------
    Coefstr=linecache.getline(filename, lign[2])
    for cursor in range(len(col)):
        if Coefstr[col[cursor]-2] != ' ':
            CY.append(float(Coefstr[col[cursor]-2:col[cursor]+byte]))
        else:        
            CY.append(float(Coefstr[col[cursor]-1:col[cursor]+byte]))
    
    del Coefstr
        
    # ----- roll moment coefficients -------
    Coefstr=linecache.getline(filename, lign[3])
    
    for cursor in range(len(col)):
        if Coefstr[col[cursor]-2] != ' ':
            Cl.append(float(Coefstr[col[cursor]-2:col[cursor]+byte]))
        else:        
            Cl.append(float(Coefstr[col[cursor]-1:col[cursor]+byte]))
    
    del Coefstr
    
    # ----- pitch moment coefficients -------
    Coefstr=linecache.getline(filename, lign[4])
    
    for cursor in range(len(col)):
        if Coefstr[col[cursor]-2] != ' ':
            Cm.append(float(Coefstr[col[cursor]-3:col[cursor]+byte]))
        else:        
            Cm.append(float(Coefstr[col[cursor]-1:col[cursor]+byte]))
    
    del Coefstr
    
    # ----- yaw moment coefficients -------
    Coefstr=linecache.getline(filename, lign[5])
    
    for cursor in range(len(col)):
        if Coefstr[col[cursor]-2] != ' ':
            Cn.append(float(Coefstr[col[cursor]-2:col[cursor]+byte]))
        else:        
            Cn.append(float(Coefstr[col[cursor]-1:col[cursor]+byte]))
    
    del Coefstr
    
    file.close()
    
    #print(CL)
    #print(CD)
    #print(CY)
    #print(Cl)
    #print(Cm)
    #print(Cn)
    
    
    #derivative=('alpha','beta','p','q','r','a','n') #names of derivatives
    
    # ----- Put into matrices 
    
    Ab=np.array([CD,CY,CL,Cl,Cm,Cn]) # put it back in order such that (u, v, w, p, q, r)
    
    return Ab

def ReadCoef(filename,col=(11,24,37,50,63,76,128,141), lign=(49,50,51,52,53,54)):
    #function takes a list of filenames and read out all the stab files and stack them in an array
    A=readstabfile(filename[0], col, lign)
    for i in range(1,len(filename)):
        A=np.vstack((A,readstabfile(filename[i],col,lign)))
        
    if len(col)!=8:
        dim=len(A[:,1])
        h=np.zeros((dim,1))
        AA=np.hstack((A,h))
        return AA
    else:
        return A




