#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 16:27:13 2021

Script to parse h5 file and produce dark field, flat field and scan field DEN files

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
parser.add_argument("hd5")
parser.add_argument("outputDir")
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

if not os.path.isfile(ARG.hd5):
    
