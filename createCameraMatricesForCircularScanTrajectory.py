#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 16:27:13 2021

Script to create camera matrices related to the circular scan.
Detals described in the posts 
https://kulvait.github.io/KCT_doc/posts/working-with-kct-cbct-2-projective-geometry-and-camera-matrices-to-describe-ct-geometry.html
https://kulvait.github.io/KCT_doc/posts/working-with-kct-cbct-3-python-implementation-of-circular-ct-trajectory.html

Units of lenghth are by convention millimeters

The version prior 21.3.2022 the meaning of the omega was angle between X axis and ccw vector to it that is from principal point to the source.
The version after 21.3.2022 the meaning of the omega is angle between X axis and ccw vector to it that is from source to the principal point, these omegaold = omega+180.
Reason of the change is to make it simmilar to parrallel ray geometry, where direction of the rays is natural vector defining geometry, the vectors VX, VY and (source, detector) in case of the CBCT and VX, VY, ray direction form a right handed coordinate system. This is to be able to easily approximate parallel rays by CBCT and vice versa.

@author: Vojtěch Kulvait
@license: GNU GPL3

"""

import argparse
import numpy as np
import json
import os
import sys
from denpy import DEN;

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("outputMatrixFile")
parser.add_argument("--source-to-isocenter", type=float, default=749.0)
parser.add_argument("--source-to-detector", type=float, default=1198.0)
parser.add_argument("--projection-sizex", type=int, default=616, help="PX dimension of detector in pixel count, defaults to 616.")
parser.add_argument("--projection-sizey", type=int, default=480, help="PY dimension of detector in pixel count, defaults to 480.")
parser.add_argument("--pixel-sizex", type=float, default=0.616, help="X pixel size, defaults to 0.616.")
parser.add_argument("--pixel-sizey", type=float, default=0.616, help="Y pixel size, defaults to 0.616.")
parser.add_argument("--pixel-offsetx", type=float, default=0., help="Principal point of the detector will have PX=0.5*projectionSizeX-0.5+pixelOffsetX")
parser.add_argument("--pixel-offsety", type=float, default=0., help="Principal point of the detector will have PY=0.5*projectionSizeY-0.5+pixelOffsetY.")
parser.add_argument("--number-of-angles", type=int, default=360, help="Number of views of the circular trajectory, .")
parser.add_argument("--omega-zero", type=float, default=0, help="Initial angle omega in degrees.")
parser.add_argument("--omega-angular-range", type=float, default=360, help="This is an angle in degrees, along which possitions are distributed.")
parser.add_argument("--endpoint", action="store_true", default=False, help="If specified include omegaZero+omegaAngularRange as a endpoint of the discretization, if not specified the end point is not included by default as it often coincide with start point."
parser.add_argument("--force", action="store_true")
parser.add_argument("--write-params-file", action="store_true")
parser.add_argument('--_json-message', default="Created using KCT script createCameraMatricesForCircularScanTrajectory.py", help=argparse.SUPPRESS)
ARG = parser.parse_args()

def checkFileExistence(f):
	if os.path.exists(f):
		if os.path.isfile(f):
			if ARG.force:
				os.remove(f)
				return
			else:
				print("File %s exist, add --force to proceed."%f)
				sys.exit(1)
		else:
			print("Path %s exist remove to proceed.", f)
			sys.exit(1)


checkFileExistence(ARG.outputMatrixFile)
paramsFile="%s.params"%ARG.outputMatrixFile
if ARG.write_params_file:
	checkFileExistence(paramsFile)

#just to debug
#import denpy
#print("Denpy version is %s"%denpy.__version__)

#zeroToSource is normally called source  to isocenter
def sourcePosition(omega, zeroToSource):
	return(zeroToSource*np.array([-np.cos(omega), -np.sin(omega), 0.], dtype=np.float64))

def X1(omega):
	X1=np.array([[np.sin(omega), -np.cos(omega), 0., 0.],
			  [np.cos(omega), np.sin(omega), 0., 0.],
			  [0,0,-1.0,0.],
			  [0,0,0.0,1.0]], dtype=np.float64)
	return X1;

def X2(sourcePosition):
	s = sourcePosition
	X2=np.array([[1, 0, 0., 0.],
		  [0, 1, 0., (np.sqrt(s[0]*s[0]+s[1]*s[1]))],
		  [0,0,1.0,s[2]],
		  [0,0,0.0,1.0]], dtype=np.float64)
	return X2;

def E(sourceToDetector):
	I=sourceToDetector;
	E=np.array([[1, 0, 0., 0.],
		  [0., 0., 1., 0],
		  [0., 1./I, 0., 0.]], dtype=np.float64)
	return E;

def A1(pixelSizeX, pixelSizeY):
	return np.array([[1.0/pixelSizeX, 0., 0.],[0, 1/pixelSizeY, 0.],[0.,0.,1.0]], dtype=np.float64)

def A2(projectionSizeX, projectionSizeY):
	return np.array([[1.0, 0., (projectionSizeX-1.0)/2.0],[0, 1.0, (projectionSizeY-1.0)/2.0],[0.,0.,1.0]], dtype=np.float64)

def A3(offsetX, offsetY):
	return np.array([[1.0, 0., offsetX],[0, 1.0, offsetY],[0.,0.,1.0]], dtype=np.float64)

I=ARG.source_to_isocenter
A=ARG.source_to_detector
M=float(ARG.projection_sizey)
N=float(ARG.projection_sizex)
PX=ARG.pixel_sizex
PY=ARG.pixel_sizey
VIEWCOUNT = ARG.number_of_angles
OMEGA = ARG.omega_zero*np.pi/180.0
RANGE = ARG.omega_angular_range*np.pi/180.0
OMEGAINCREMENT = RANGE/VIEWCOUNT

#Let's create specified set of projection matrices as np.array
CameraMatrices = np.zeros((0,4,3), dtype=np.float64)
directionAngles = np.linspace(OMEGA, OMEGA + RANGE, num=VIEWCOUNT, endpoint=ARG.endpoint)

for i in range(len(directionAngles)):
	omega = directionAngles[i]
	s = sourcePosition(omega, I) 
	_A3=A3(ARG.pixel_offsetx, ARG.pixel_offsety)
	_A2=A2(N, M)
	_A1 = A1(PX, PY)
	_E = E(A)
	_X2 = X2(s) 
	_X1 = X1(omega)
	CM = _A3.dot(_A2).dot(_A1).dot(_E).dot(_X2).dot(_X1)
	CameraMatrices = np.dstack((CameraMatrices, CM))

DEN.storeNdarrayAsDEN(ARG.outputMatrixFile, CameraMatrices, ARG.force)
#Solution from https://stackoverflow.com/a/55114771
if ARG.write_params_file:
	with open(paramsFile, 'w') as f:
		json.dump(ARG.__dict__, f, indent=2, sort_keys=True)
