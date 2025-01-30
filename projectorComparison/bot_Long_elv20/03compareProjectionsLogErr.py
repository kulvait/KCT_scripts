#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
parser.add_argument("--title", default="Comparison of projectors")
ARG = parser.parse_args()

import matplotlib
import scienceplots
if "png" in ARG:
	matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})
plt.style.use(['science', 'nature', 'no-latex'])
#col=['#0089f7', '#98e6ff', '#c55e46', '#91b657', '#FF0000'];

#default_cycler = (cycler(color=col) +
#                  cycler(linestyle=['-', ':', '--', '-.', "-"]))

plt.rc('lines', linewidth=1)
#plt.rc('axes', prop_cycle=default_cycler)

with open(ARG.procesCsv, newline='') as csvfile:
	spamreader = csv.reader(csvfile, delimiter='\t')
	csv={}
	header=next(spamreader)
	for x in header:
		csv[x] = np.zeros(0);
	for row in spamreader:
		for i in range(len(row)):
			csv[header[i]] = np.append(csv[header[i]], float(row[i]))

fig = plt.figure(dpi=150, tight_layout=True)
#A4 figure
#fig.set_size_inches(8.27, 11.69, forward=True)
fig.set_size_inches(8.27, 3.5, forward=True)

for x in header:
	if x not in ["Angle", "Baseline"]:#, "Siddon1", "Siddon2", "Siddon4", "Siddon16", "Siddon32", "Siddon64", "Siddon128", "Siddon256"]:#, "CVP", "TT"]:
		#plt.plot(csv["Angle"], csv[x]/csv["Baseline"], "o", label=x) 
		if x.startswith("Siddon"):
			line, = plt.plot(csv["Angle"], 100*csv[x]/csv["Baseline"], linestyle=":", color="black", label=None, linewidth="0.5") 
			plt.text(csv["Angle"][-1], (100 * csv[x] / csv["Baseline"])[-1], x, fontsize=14, verticalalignment='center')
		else:
			line, = plt.plot(csv["Angle"], 100*csv[x]/csv["Baseline"], label=x, linestyle="-") 

plt.yscale('log')
# Format y-axis ticks as 1, 10, 100 instead of 10^0, 10^1, 10^2
plt.yticks([10**i for i in range(-5, 3)], [str(10**i) for i in range(-5, 3)])
plt.xlabel("Angle in degrees")
plt.ylabel("Relative error [%]")
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


#print(csv)
