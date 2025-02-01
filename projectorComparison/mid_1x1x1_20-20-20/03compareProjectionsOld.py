#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import csv 
import argparse
import os
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument("procesCsv")
parser.add_argument("--png", default=argparse.SUPPRESS)
parser.add_argument("--title", default="Comparison of projectors")
ARG = parser.parse_args()

import matplotlib
if "png" in ARG:
	matplotlib.use('Agg')
import matplotlib.pyplot as plt

with open(ARG.procesCsv, newline='') as csvfile:
	spamreader = csv.reader(csvfile, delimiter='\t')
	csv={}
	header=next(spamreader)
	for x in header:
		csv[x] = np.zeros(0);
	for row in spamreader:
		for i in range(len(row)):
			csv[header[i]] = np.append(csv[header[i]], float(row[i]))

fig = plt.figure(dpi=300, tight_layout=True)
#A4 figure
#fig.set_size_inches(8.27, 11.69, forward=True)
fig.set_size_inches(8.27, 3.5, forward=True)

for x in header:
	if x not in ["Angle", "Baseline", "Siddon1", "Siddon2", "Siddon4", "Siddon32", "Siddon64", "Siddon128", "Siddon256", "CVP_relaxed", "CVP"]:#, "CVP", "TT"]:
		#plt.plot(csv["Angle"], csv[x]/csv["Baseline"], "o", label=x) 
		plt.plot(csv["Angle"], 100*csv[x]/csv["Baseline"], "-", label=x) 

plt.xlabel("Angle in degrees")
plt.ylabel("Relative error [%]")
plt.title(ARG.title)
plt.legend()
plt.show()






if "png" in ARG:
	plt.savefig(ARG.png)
else:
	plt.show()
