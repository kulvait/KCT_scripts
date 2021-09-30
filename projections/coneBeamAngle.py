#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 2021

Script to output view dependent cone beam angles for given position

Units of lenghth are by convention millimeters

@author: Vojtěch Kulvait
@license: GNU GPL3
"""

import argparse
import numpy as np
from denpy import DEN;

parser = argparse.ArgumentParser()
parser.add_argument("cameraMatrices")
parser.add_argument("x", type=float)
parser.add_argument("y", type=float)
parser.add_argument("z", type=float)

ARG = parser.parse_args()


def normalizeVector(v):
	norm = np.linalg.norm(v)
	return(v/norm)

def sourcePosition(CM):
	CM3x3 = np.delete(CM, 3, 1)
	CM1x3 = CM[:,3]
	ICM3x3 = np.linalg.inv(CM3x3)
	source = ICM3x3.dot(-CM1x3)
	return(source)

def normalToDetector(CM):
	#First I have to orthogonalize third vector towards first two
	source = sourcePosition(CM)
	#Get 0 in a coordinates connected to zero
	zeroLoc = -source
	CM3x3 = np.delete(CM, 3, 1)
	n = normalizeVector(CM3x3[2,:])
	if n.dot(zeroLoc) < 0:
		n = -n
	return(n)

def orthogonalPart(v, w):
	ort = w * v.dot(w)/w.dot(w)
	return v - ort

def VN(CM):
	CM3x3 = np.delete(CM, 3, 1)
	return CM3x3[2,:]

def VX(CM):
	CM3x3 = np.delete(CM, 3, 1)
	vn = VN(CM)
	return orthogonalPart(CM3x3[0,:], vn)

def VY(CM):
	CM3x3 = np.delete(CM, 3, 1)
	vn = VN(CM)
	return orthogonalPart(CM3x3[1,:], vn)

def coneAngle(CM, x, y, z):
	source = sourcePosition(CM)
	point = np.array([x,y,z])
	point_loc = point - source
	point_nrm = normalizeVector(point_loc)
	n=normalToDetector(CM)
	product = n.dot(point_nrm)
	if product > 1.0:
		product = 1.0
	return np.arccos(product)
	

def azimuthalAngle(CM, x, y, z):
	source = sourcePosition(CM)
	point = np.array([x,y,z])
	point_loc = point - source
	vx=VX(CM)
	vy=VY(CM)
	vn=VN(CM)
	PN = point_loc.dot(vn)
	PX = point_loc.dot(vx)/PN
	PY = point_loc.dot(vy)/PN
	if PX*PX+PY*PY<1e-9 :
		return 0
	at = np.arctan2(PY, PX)
	return(at)

#Coordinate on the detector relative to the principal zero
def APX(CM, x, y, z):
	source = sourcePosition(CM)
	point = np.array([x,y,z])
	point_loc = point - source
	vx=VX(CM)
	vn=VN(CM)
	PN = point_loc.dot(vn)
	PX = point_loc.dot(vx)/PN
	return(PX)


#Coordinate on the detector relative to the principal zero
def APY(CM, x, y, z):
	source = sourcePosition(CM)
	point = np.array([x,y,z])
	point_loc = point - source
	vy=VY(CM)
	vn=VN(CM)
	PN = point_loc.dot(vn)
	PY = point_loc.dot(vy)/PN
	return(PY)

#Coordinate on the detector endoded by camera matrix
def PXOFFSET(CM, x, y, z):
	source = sourcePosition(CM)
	point = np.array([x,y,z])
	point_loc = point - source
	CM3x3 = np.delete(CM, 3, 1)
	az=CM3x3[2,:].dot(point_loc)
	ax=CM3x3[0,:].dot(point_loc)
	return(ax/az)

#Coordinate on the detector endoded by camera matrix
def PYOFFSET(CM, x, y, z):
	source = sourcePosition(CM)
	point = np.array([x,y,z])
	point_loc = point - source
	CM3x3 = np.delete(CM, 3, 1)
	az=CM3x3[2,:].dot(point_loc)
	ay=CM3x3[1,:].dot(point_loc)
	return(ay/az)

header = DEN.readHeader(ARG.cameraMatrices)
print("Angle\tConeBeamAngle\tazimuthalAngle\tPX0\tPY0\tPX\tPY")
for i in range(header["zdim"]):
	CM = DEN.getDoubleFrame(ARG.cameraMatrices, i)
	CONEANGLE = coneAngle(CM, ARG.x, ARG.y, ARG.z)
	AZIANGLE = azimuthalAngle(CM, ARG.x, ARG.y, ARG.z)
	source = sourcePosition(CM)
	n = normalToDetector(CM)
	PXO=PXOFFSET(CM, ARG.x, ARG.y, ARG.z)
	PYO=PYOFFSET(CM, ARG.x, ARG.y, ARG.z)
	print("%d\t%f\t%f\t%f\t%f\t%f\t%f"%(i, CONEANGLE*180/np.pi, AZIANGLE*180/np.pi, APX(CM, ARG.x, ARG.y, ARG.z), APY(CM, ARG.x, ARG.y, ARG.z), PXO, PYO))
