#!/usr/bin/env python

"""
@Author: Vojtech Kulvait
@Date: 2020-2022
@License: GPL3 
"""
import argparse
from denpy import DEN
import cv2
import numpy as np

parser = argparse.ArgumentParser(description = 'Morphological operations on DEN files.')
parser.add_argument('input_file', type=str, help = 'Input DEN file')
parser.add_argument('output_file', type=str, help = 'Output DEN file')
parser. add_argument("--erode", action="store_true")
parser. add_argument("--dilate", action="store_true")
parser. add_argument("--open", action="store_true")
parser. add_argument("--close", action="store_true")

ARG = parser.parse_args()


def circular_kernel(sz):
	x,y=np.ogrid[-sz/2:sz/2, -sz/2:sz/2]
	kernel = x*x + y*y <= (sz*sz)/(2*2)
	return kernel.astype(np.uint8)



h = DEN.readHeader(ARG.input_file)
DEN.writeEmptyDEN(ARG.output_file, h["dimspec"], force=True)
for k in range(h["dimspec"][-1]):
	f=DEN.getFrame(ARG.input_file, k)
	if ARG.erode:
		kernel = circular_kernel(10)
		f = cv2.erode(f, kernel, iterations = 1)
	if ARG.dilate:
		kernel = circular_kernel(10)
		f = cv2.dilate(f, kernel, iterations = 1)
	if ARG.open:
		kernel = circular_kernel(5)
		f = cv2.morphologyEx(f, cv2.MORPH_OPEN, kernel)
	if ARG.close:
		kernel = circular_kernel(10)
		f = cv2.morphologyEx(f, cv2.MORPH_CLOSE, kernel)
	DEN.writeFrame(ARG.output_file, k, f, True)
