#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created 2021

@author: Vojtech Kulvait

This script processes data to get raw DEN volumes

Input dir shall be the directory with the dicom files from the single scan or single series of scan with the same SeriesInstanceUID or alternativelly IrradiationEventUID

This script can process perfusion data as well producing series of DEN files

Normal scan output will be INFO file with basics characteristics, Series_00.den with the raw data and Series_00.tick with the timing information that contain times of acquisitions of individual slices in seconds relative to the first slice acquired  based on the DICOM

By default, output will be in Hounsfield units

To directly convert these to the raw values, with installed dentk package its possible using --convert-hu-to-raw option
--water-value shall be specified
--base2 as Siemens and maybe other scanners use base2 encoding, minimal HU value is -1024, that might have produced negative values in raw conversions if expected that the minimum value is 1000, this setting will put -1024 as a minimal HU value for conversion

License: GNU GPL3

"""
import argparse
import os
import pydicom
import numpy as np
from denpy import DICOM
from denpy import DEN
import copy
import sys
import subprocess as sh


parser = argparse.ArgumentParser()
parser.add_argument("inputDir")
parser.add_argument("outputDir")
parser.add_argument("--force", action="store_true")
parser.add_argument("--dicom-file-suffix", type=str, default="ima")
parser.add_argument("--convert-hu-to-raw", action="store_true")
parser.add_argument("--water-value", default="0.027")
parser.add_argument("--base2", action="store_true")
ARG = parser.parse_args()

dicomFiles = DICOM.getDicomFiles(ARG.inputDir, [ARG.dicom_file_suffix])
dicoms = [pydicom.dcmread(x) for x in dicomFiles]
if len(dicoms) == 0:
	print("No DICOM files found in the %s directory!"%ARG.inputDir)
	sys.exit(2)
if not hasattr(dicoms[0], "SeriesDescription") or dicoms[0].SeriesDescription.startswith('DynMulti4D  1.5  Br36') or dicoms[0].SeriesDescription.startswith('Perfusion  1.5  Br36'):
	dicomUIDs = list(set([x.IrradiationEventUID for x in dicoms]))
else:
	dicomUIDs = list(set([x.SeriesInstanceUID for x in dicoms]))
if len(dicomUIDs) != 1:
	raise ValueError(
		"The number of different UIDs in folder should be 1 but is %d" % len(dicomUIDs))




#Originally the script was working with AcquisitionTime but Philips scanners seem to put it to the beginning of the acquisition
#For Siemens scanners AcquisitionTime is the same as ContentTime and for Philips ContentTime describes scan better
start_time = min([DICOM.dateAndTimeToSeconds(x.ContentDate, x.ContentTime) for x in dicoms])
end_time = max([DICOM.dateAndTimeToSeconds(x.ContentDate, x.ContentTime) for x in dicoms])
print("Duration of the perfusion scan was %0.3fs."%(end_time-start_time))
if end_time - start_time > 100:
	raise ValueError("Length of acquisition should not exceed 100 seconds but it is %f" % (
		end_time - start_time))

INF = []
for i in range(len(dicoms)):
	time = DICOM.dateAndTimeToSeconds(dicoms[i].ContentDate, dicoms[i].ContentTime)
	t = (dicomFiles[i], time - start_time, float(dicoms[i].SliceLocation), dicoms[i])
	INF.append(t)

INX = copy.deepcopy(INF)

if not os.path.exists(ARG.outputDir):
	os.makedirs(ARG.outputDir)
elif not ARG.force:
	sys.exit("Directory %s exists, use --force to force create."%ARG.outputDir)

INF = sorted(INF, key=lambda x: (x[2], x[1]))

z_directions = [x[2] for x in INF]
z_directions = list(set(z_directions))
z_directions = np.sort(z_directions)



numberOfVolumes = int(len(INF) / len(z_directions))

volumeStart = [val for ind, val in enumerate(INF) if ind % numberOfVolumes == 0]
volumeEnd = [val for ind, val in enumerate(INF) if (ind+1) % numberOfVolumes == 0]
timesStart = [x[1] for x in volumeStart]
timesEnd = [x[1] for x in volumeEnd]
acquisition_start = np.min(timesStart)
acquisition_end = np.max(timesEnd)
covered_start = np.max(timesStart)
covered_end = np.min(timesEnd)
total_time=acquisition_end-acquisition_start
covered_time=covered_end-covered_start
z_spacing = np.abs(volumeStart[0][2]-volumeStart[1][2])#Based on volumeStart 0 and 1 elements positions difference
z_spacing = np.mean([j-i for i, j in zip(z_directions[:-1], z_directions[1:])])

#Estimate optimal pause size to fit particular scan


dimensions = dicoms[3].pixel_array.shape
pixelSpacing = dicoms[3].PixelSpacing
sliceThickness = dicoms[3].SliceThickness

info = open(os.path.join(ARG.outputDir,'INFO'), 'w')
info.write("ACQUISITIONSTART=%f\n"%acquisition_start)
info.write("ACQUISITIONEND=%f\n"%acquisition_end)
info.write("TOTALTIME=%f\n"%(acquisition_end-acquisition_start))
info.write("PIXELSPACINGX=%f\n"%pixelSpacing[0])
info.write("PIXELSPACINGY=%f\n"%pixelSpacing[1])
#There is no implementation of non aligned voxels yet in KCT so VOXELSIZEX and VOXELSIZEY being PIXELSPACINGX and PIXELSPACINGY is a good idea
info.write("VOXELSIZEX=%f\n"%pixelSpacing[0])
info.write("VOXELSIZEY=%f\n"%pixelSpacing[1])
#Here I use z_spacing as a equivalent for VOXELSIZEZ as there is not yet implementation of overlapping voxels in a z direction in KCT
info.write("VOXELSIZEZ=%f\n"%z_spacing)
#Volume dimenxions
info.write("VOLUMESIZEX=%d\n"%dimensions[0])
info.write("VOLUMESIZEY=%d\n"%dimensions[1])
info.write("VOLUMESIZEZ=%d\n"%len(z_directions))
info.write("SLICETHICKNESS=%f\n"%sliceThickness)
info.write("ZSPACING=%f\n"%z_spacing)
info.close()


for i in range(numberOfVolumes):
	print("Preparing volume Series_%02d." % (i))
	out_im = np.zeros([dimensions[0], dimensions[1],
                    len(z_directions)], dtype="<f4")
	dimfr = np.zeros([1, 5, len(z_directions)], dtype="<f4")
	print("Processing start file %s with time %f."%(INF[i][0], INF[i][1]))
	for k in range(len(z_directions)):
		index = i + k * numberOfVolumes
		intercept = INF[index][3].RescaleIntercept
		slope = INF[index][3].RescaleSlope
		out_im[:, :, k] = INF[index][3].pixel_array * slope + intercept
		dimfr[:, :, k] = np.reshape(np.array(
			(pixelSpacing[0], pixelSpacing[1], z_spacing, INF[index][2], INF[index][1]), dtype="<f4"), (1, 5))
	out_im = out_im.astype(np.float32)
	denfile = os.path.join(ARG.outputDir, "Series_%02d.den" % (i))
	DEN.storeNdarrayAsDEN(denfile, out_im, force=True)
	DEN.storeNdarrayAsDEN(os.path.join(
		ARG.outputDir, "Series_%02d.tick" % (i)), dimfr, force=True)
	if ARG.convert_hu_to_raw:
		CMD=["dentk-fromhu"]
		if ARG.base2:
			CMD.append("--base2")
		CMD.append("--water-value")
		CMD.append("%s"%(ARG.water_value))
		CMD.append(denfile)
		CMD.append(denfile)
		print(" ".join(CMD))
		proc = sh.call(CMD)
