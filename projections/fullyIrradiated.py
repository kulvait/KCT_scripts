#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 2021

Script to output maximal dimensions of the cube to be fully irradiated

Units of lenghth are by convention millimeters

@author: Vojtěch Kulvait
@license: GNU GPL3
"""

import argparse
import numpy as np
from denpy import DEN;

parser = argparse.ArgumentParser()
parser.add_argument("cameraMatrices")
parser.add_argument('--voxel-sizex', help="X component of the vector to the corner.",
                    type=float, default=1.0)
parser.add_argument('--voxel-sizey', help="Y component of the vector to the corner.",
                    type=float, default=1.0)
parser.add_argument('--voxel-sizez', help="Z component of the vector to the corner.",
                    type=float, default=1.0)
parser.add_argument('--volume-centerx', help="X coordinate of the volume center in mm, defaults to 0.00.",
                    type=float, default=0.0)
parser.add_argument('--volume-centery', help="Y coordinate of the volume center in mm, defaults to 0.00.",
                    type=float, default=0.0)
parser.add_argument('--volume-centerz', help="Z coordinate of the volume center in mm, defaults to 0.00.",
                    type=float, default=0.0)
parser.add_argument('--projection-sizex', help="X dimension of detector in pixel count, defaults to 616.",
                    type=int, default=616)
parser.add_argument('--projection-sizey', help="Y dimension of detector in pixel count, defaults to 480.",
                    type=int, default=480)
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
	return(np.arctan2(PY, PX))

def APX(CM, x, y, z):
	source = sourcePosition(CM)
	point = np.array([x,y,z])
	point_loc = point - source
	vx=VX(CM)
	vn=VN(CM)
	PN = point_loc.dot(vn)
	PX = point_loc.dot(vx)/PN
	return(PX)


def APY(CM, x, y, z):
	source = sourcePosition(CM)
	point = np.array([x,y,z])
	point_loc = point - source
	vy=VY(CM)
	vn=VN(CM)
	PN = point_loc.dot(vn)
	PY = point_loc.dot(vy)/PN
	return(PY)

def project(CM, x, y, z):
	point = np.array([x,y,z, 1.0])
	det_point = CM.dot(point)
	PX=det_point[0]/det_point[2]
	PY=det_point[1]/det_point[2]
	return(np.array([PX, PY]))

#If all matrices in CMFile project volumeCenter+d*length for d in directions onto the detector
def isConsistent(CMFile, volumeCenter, directions, length):
	header = DEN.readHeader(CMFile)
	for i in range(header["zdim"]):
		CM = DEN.getDoubleFrame(CMFile, i)
		for k in range(directions.shape[1]):
			P = project(CM, float(volumeCenter[0] + length*directions[0,k]), float(volumeCenter[1] + length*directions[1,k]), float(volumeCenter[2] + length*directions[2,k]))
			if P[0] < -0.5 or P[1] < -0.5 or P[0] > ARG.projection_sizex + 0.5 or P[1] > ARG.projection_sizey + 0.5:
				return(False)
	return(True)

volumeCenter=np.array([ARG.volume_centerx, ARG.volume_centery, ARG.volume_centerz])
volumeCenter.shape=(3,1)
v1=np.array([ARG.voxel_sizex, ARG.voxel_sizey, ARG.voxel_sizez])
v2=np.array([-v1[0], v1[1], v1[2]])
v3=np.array([v1[0], -v1[1], v1[2]])
v4=np.array([-v1[0], -v1[1], v1[2]])
v5=np.array([v1[0], v1[1], -v1[2]])
v6=np.array([-v1[0], v1[1], -v1[2]])
v7=np.array([v1[0], -v1[1], -v1[2]])
v8=np.array([-v1[0], -v1[1], -v1[2]])
voxelDirections = np.c_[v1, v2, v3, v4, v5, v6, v7, v8]
length = 1.0
if not isConsistent(ARG.cameraMatrices, volumeCenter, voxelDirections, 0.0):
	print("Volume center does not project onto the detector")
knownFit = 0.0
while isConsistent(ARG.cameraMatrices, volumeCenter, voxelDirections, length):
	knownFit = length
	length = length*2.0
knownUnfit = length

length = knownFit + 0.5*(knownUnfit-knownFit)

for n in range(50):
	consistent = isConsistent(ARG.cameraMatrices, volumeCenter, voxelDirections, length)
	if consistent == True:
		if knownFit < length:
			knownFit = length
	else:
		if knownUnfit > length:
			knownUnfit = length
	length = knownFit + 0.5*(knownUnfit-knownFit)

print("Optimal numbers of voxels are [%f, %f, %f]"%(knownFit*2, knownFit*2, knownFit*2))
directions=voxelDirections
for k in range(directions.shape[1]):
	print("Corner %d has dimensions [%f, %f, %f]."%(k, volumeCenter[0] + knownFit*directions[0,k], volumeCenter[1] + knownFit*directions[1,k], volumeCenter[2] + knownFit*directions[2,k]))
