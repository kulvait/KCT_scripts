#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#pip install git+https://github.com/garrettj403/SciencePlots
from cycler import cycler
import csv 
import argparse
import os
import numpy as np

def convertToGray(rgb):
	h = rgb.lstrip('#')
	col = tuple(int(h[i:i+2], 16) for i in (0, 2, 4))
	print(col)
	gr = 0.2989*col[0]+0.5870*col[1]+0.1140*col[2]
	gr = int(gr)
	return "#%02x%02x%02x"%(gr, gr, gr)

parser = argparse.ArgumentParser()
parser.add_argument("procesCsv")
parser.add_argument("--png", default=argparse.SUPPRESS)
parser.add_argument("--title", default="Errors projecting single $1 \\times 1 \\times 5$ voxel without elevation")
ARG = parser.parse_args()

import matplotlib
import scienceplots
if "png" in ARG:
        matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})
plt.style.use(['science', 'nature'])
col=['#0089f7', '#98e6ff', '#c55e46', '#91b657', '#FF0000'];

default_cycler = (cycler(color=col) +
                  cycler(linestyle=['-', ':', '--', '-.', "-"]))

#if "png" in ARG:
#	matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#plt.style.use(['science', 'nature'])
#https://medialab.github.io/iwanthue/
#col=['#66CCEE', '#228833', convertToGray('#DDAA33'), convertToGray('#004488')];
#col=['#0C5DA5', '#00B945', '#FF2C00', convertToGray('#004488')];
#col=['#0089f7', '#3db42d', '#be37bc', '#ff6e63'];
#col=['#0089f7', '#98e6ff', '#f688f2', '#addc51'];
#col=['#0089f7', '#98e6ff', '#a1c449', '#d44770'];
#col=['#0089f7', '#98e6ff', '#c55e46', '#91b657'];
#col=['#5b2ec3', '#6a8c26', '#c34455', '#6D8288'];
#col=['#0077BB', '#009988', '#EE3311', '#BBBBBB'];

plt.rc('lines', linewidth=1)
plt.rc('axes', prop_cycle=default_cycler)


with open(ARG.procesCsv, newline='') as csvfile:
	spamreader = csv.reader(csvfile, delimiter='\t')
	csv={}
	header=next(spamreader)
	for x in header:
		csv[x] = np.zeros(0);
	for row in spamreader:
		for i in range(len(row)):
			csv[header[i]] = np.append(csv[header[i]], float(row[i]))

#fig = plt.figure(dpi=300, tight_layout=True)
fig = plt.figure(dpi=150, tight_layout=True)
#A4 figure
#fig.set_size_inches(8.27, 11.69, forward=True)

fig.set_size_inches(8.27, 3.5, forward=True)

#for x in header:
#	if x not in ["Angle", "Baseline", "Siddon1", "Siddon2", "Siddon4", "Siddon8", "Siddon16", "Siddon32", "Siddon64", "Siddon256"]:#, "CVP", "TT"]:
#		#plt.plot(csv["Angle"], csv[x]/csv["Baseline"], "o", label=x) 
#		plt.plot(csv["Angle"], 100*csv[x]/csv["Baseline"], label=x) 

line, = plt.plot(csv["Angle"], 100*csv["CVE"]/csv["Baseline"], label="CVP")
line, = plt.plot(csv["Angle"], 100*csv["CVE_relaxed"]/csv["Baseline"], label="CVP relaxed")
line, = plt.plot(csv["Angle"], 100*csv["TT"]/csv["Baseline"], label="TT")
line, = plt.plot(csv["Angle"], 100*csv["Siddon128"]/csv["Baseline"], label="Siddon128")

plt.xlabel("Angle in degrees")
plt.ylabel("Relative error [\%]")
#plt.xlim((0, 360))
plt.xticks(np.arange(0.0,361,45))
plt.title(ARG.title)
plt.legend(frameon=True, loc="upper right")
#plt.xlim((200, 300))
#plt.ylim((0, 0.05))

if "png" in ARG:
	plt.savefig(ARG.png)
else:
	plt.show()


