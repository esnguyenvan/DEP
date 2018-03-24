# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 17:11:31 2018

@author: e.nguyen-van

To read a STL file and plot in mayavi
First created by Junwei Huang @ Toronto, Feb 26, 2013
"""

from numpy import *
from mayavi import mlab

STLfile="test.stl"
f=open(STLfile,'r')

x=[]
y=[]
z=[]

for line in f:
	strarray=line.split()
	if strarray[0]=='vertex':
		x=append(x,double(strarray[1]))
		y=append(y,double(strarray[2]))
		z=append(z,double(strarray[3]))

triangles=[(i, i+1, i+2) for i in range(0, len(x),3)]

mlab.triangular_mesh(x, y, z, triangles)
mlab.show()