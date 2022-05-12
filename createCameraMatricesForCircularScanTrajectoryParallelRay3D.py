#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 16:27:13 2021

Script to create camera matrices related to the circular scan in Parallel ray geometry
See also
https://kulvait.github.io/KCT_doc/posts/working-with-kct-cbct-5-parallel-beam-geometry.html

Units of lenghth are by convention millimeters

@author: Vojtěch Kulvait
@license: GNU GPL3

"""

import argparse
import numpy as np
import scipy.io
import json
import os
import sys
from denpy import DEN

#Defaults will be taken from https://confluence.desy.de/display/P5I/Detectors, KIT CMOS, no magnification
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("outputMatrixFile")
parser.add_argument("--projection-sizex", type=int, default=5120, help="PX dimension of detector in pixel count, defaults to 5120.")
parser.add_argument("--projection-sizey", type=int, default=3840, help="PY dimension of detector in pixel count, defaults to 3840.")
parser.add_argument("--pixel-sizex", type=float, default=0.064, help="X pixel size, defaults to 0.0064.")
parser.add_argument("--pixel-sizey", type=float, default=0.064, help="Y pixel size, defaults to 0.0064.")
parser.add_argument("--detector-centerx", type=float, default=0., help="Coordinate x of the center of the detector, defaults to 0.0.")
parser.add_argument("--detector-centery", type=float, default=0., help="Coordinate y of the center of the detector, defaults to 0.0.")
parser.add_argument("--detector-centerz", type=float, default=0., help="Coordinate z of the center of the detector, defaults to 0.0.")
parser.add_argument("--number-of-angles", type=int, default=5001, help="Number of views of the circular trajectory, .")
parser.add_argument("--omega-zero", type=float, default=0, help="Initial angle omega in degrees.")
parser.add_argument("--omega-angular-range", type=float, default=360, help="This is an angle in degrees, along which possitions are distributed.")
parser.add_argument("--endpoint", action="store_true", default=False, help="If specified include omegaZero+omegaAngularRange as a endpoint of the discretization, if not specified the end point is not included by default as it often coincide with start point.")
parser.add_argument("--angles-mat", type=str, default=None, help="Sequence of angles stored in mat file.")
parser.add_argument("--force", action="store_true")
parser.add_argument("--write-params-file", action="store_true")
parser.add_argument('--_json-message', default="Created using KCT script createCameraMatricesForCircularScanTrajectoryParallelRay3D.py", help=argparse.SUPPRESS)
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

#Direction so that omega is ccw to the X axis, returns unit vector
def rayDirection(omega):
	return(np.array([np.cos(omega), np.sin(omega), 0.], dtype=np.float64))

M=float(ARG.projection_sizey)
N=float(ARG.projection_sizex)
PX=ARG.pixel_sizex
PY=ARG.pixel_sizey

if ARG.angles_mat is not None:
	matlab_dic = scipy.io.loadmat(ARG.angles_mat)
	directionAngles = matlab_dic["angles"]
else:
	VIEWCOUNT = ARG.number_of_angles
	OMEGA = ARG.omega_zero*np.pi/180.0
	RANGE = ARG.omega_angular_range*np.pi/180.0
	OMEGAINCREMENT = RANGE/VIEWCOUNT
	#Let's create specified set of projection matrices as np.array
	directionAngles = np.linspace(OMEGA, OMEGA + RANGE, num=VIEWCOUNT, endpoint=ARG.endpoint)

CameraMatrices = np.zeros((0,2,4), dtype=np.float64)
for i in range(len(directionAngles)):
	omega = float(directionAngles[i])
	VR = rayDirection(omega)
	VX = np.array([np.sin(omega)*PX, -np.cos(omega)*PX, 0.0], dtype=np.float64)
	a = np.array([np.sin(omega)/PX, -np.cos(omega)/PX, 0.0], dtype=np.float64)
	VY = np.array([0.0, 0.0, -PY], dtype=np.float64)
	b = np.array([0.0, 0.0, -1.0/PY], dtype=np.float64)
	detectorCenter=np.array([ARG.detector_centerx, ARG.detector_centery, ARG.detector_centerz], dtype=np.float64)
	px0 = N * 0.5 - 0.5 - detectorCenter.dot(a)
	py0 = M * 0.5 - 0.5 - detectorCenter.dot(b)
	CM = np.array([np.append(a, px0), np.append(b,py0)])
	CM.shape=(1,2,4)
	CameraMatrices = np.concatenate((CameraMatrices, CM))

DEN.storeNdarrayAsDEN(ARG.outputMatrixFile, CameraMatrices, ARG.force)
#Solution from https://stackoverflow.com/a/55114771
if ARG.write_params_file:
	with open(paramsFile, 'w') as f:
		json.dump(ARG.__dict__, f, indent=2, sort_keys=True)
