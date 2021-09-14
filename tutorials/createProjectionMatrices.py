#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 16:27:13 2021

@author: Vojtěch Kulvait
@license: GNU GPL3

"""

import numpy as np
from denpy import DEN;

#zeroToSource is normally called source  to isocenter
def sourcePosition(omega, zeroToSource):
	return(zeroToSource*np.array([np.cos(np.pi * omega/180.0), np.sin(np.pi * omega/180.0), 0.], dtype=np.float64))

def X1(omega):
	a = np.pi * omega/180.0;
	X1=np.array([[-np.sin(a), np.cos(a), 0., 0.],
			  [np.cos(a), np.sin(a), 0., 0.],
			  [0,0,-1.0,0.],
			  [0,0,0.0,1.0]], dtype=np.float64)
	return X1;

def X2(sourcePosition):
	s = sourcePosition
	X2=np.array([[1, 0, 0., 0.],
		  [0, 1, 0., -(np.sqrt(s[0]*s[0]+s[1]*s[1]))],
		  [0,0,1.0,s[2]],
		  [0,0,0.0,1.0]], dtype=np.float64)
	return X2;

def E(sourceToDetector):
	I=sourceToDetector;
	E=np.array([[1, 0, 0., 0.],
		  [0, 0., 1., 0],
		  [0,-1/I,0.0,0.]], dtype=np.float64)
	return E;

def A1(pixelSizeX, pixelSizeY):
	return np.array([[1.0/pixelSizeX, 0., 0.],[0, 1/pixelSizeY, 0.],[0.,0.,1.0]], dtype=np.float64)

def A2(projectionSizeX, projectionSizeY):
	return np.array([[1.0, 0., (projectionSizeX-1.0)/2.0],[0, 1.0, (projectionSizeY-1.0)/2.0],[0.,0.,1.0]], dtype=np.float64)

#Let's create specified set of projection matrices as np.array
I=749.0
A=1198.0
M=616.0
N=480.0
PX=0.616
PY=0.616

CameraMatrices = np.zeros((3,4,0), dtype=np.float64)
for omega in range(360):
	s = sourcePosition(omega, I) 
	_A2=A2(M, N)
	_A1 = A1(PX, PY)
	_E = E(A)
	_X2 = X2(s) 
	_X1 = X1(omega)
	CM = _A2.dot(_A1).dot(_E).dot(_X2).dot(_X1)
	CameraMatrices = np.dstack((CameraMatrices, CM))
	
DEN.storeNdarrayAsDoubleDEN("CM.den", CameraMatrices, True)