#!/usr/bin/env python

"""
@Author: Hana Haseljic, Vojtech Kulvait
@Date: 2020-2022
@License: GPL3 
"""

#Example
#perfusion-mosaic.py CBF.den CBV.den MTT.den TTP.den --mask-file ALPHA_INSIDE.den -lp 2 --name pat1 --frames 2 14 21 --lut Rainbow.csv
#perfusion-mosaic.py TTP.den TTPA.den TTPB.den --ranges 0-50 0-50 0-50 --frames 15

#
#Examples for Erlangen data
#

#python perfusion-mosaic-compare.py /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_3.50_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/TTPB.den /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_2.75_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/TTPB.den /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06FourierGauss_3.50Viz_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/TTPB.den --ttp-range-files /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_3.50_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/TTP.den /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_2.75_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/TTP.den /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06FourierGauss_3.50Viz_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/TTP.den --aif 274 237 2 --frames 14 --labels "Gauss 3.5" "Gauss 2.75" "FourierGauss 3.5"

#perfusion-mosaic-compare.py /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_3.50_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/TTPB.den /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_2.75_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/TTPB.den /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06FourierGauss_3.50Viz_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/TTPB.den --aif 274 237 2 --frames 14 --ranges 0-50 0-50 0-50

#python perfusion-mosaic-compare.py /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_3.50_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/CBFA.den /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_2.75_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/CBFA.den /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06FourierGauss_3.50Viz_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/CBFA.den --mask-file /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_3.50_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/BETA.den -lp 1.5 --aif 274 237 2 --frames 14 --labels "Gauss 3.5" "Gauss 2.75" "FourierGauss 3.5"

#python perfusion-mosaic-compare.py /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_3.50_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/CBFA.den /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_2.75_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/CBFA.den /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06FourierGauss_3.50Viz_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/CBFA.den -lp 1.5 --aif 274 237 2 --frames 14 --labels "Gauss 3.5" "Gauss 2.75" "FourierGauss 3.5"

import argparse
import glob
import os
import sys
import subprocess
from denpy import DEN
import numpy as np

from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw
from PIL import ImageEnhance

import matplotlib as plt
plt.use('Agg')
from matplotlib import pyplot
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib import colors

from skimage import exposure
from skimage.exposure import equalize_hist

import math

import csv


# Compute window based on the volume
def computeRange(files, lp, mask = False):
    ranges=[]
    for file in files:
        f = DEN.getNumpyArray(file)
        p0=0
        p1=0
        if (mask):
            m = DEN.getNumpyArray(mask)
            f = np.multiply(f, m)
        # if os.path.splitext(os.path.basename(file))[0][-1] == 'A':
        #     p0=np.percentile(f, 1.5) #based on brain data
        #     p1=np.percentile(f, 98.5)
        # elif os.path.splitext(os.path.basename(file))[0][-1] == 'B':
        #     p0 = np.percentile(f, 0.5) #based on brain and liver data
        #     p1 = np.percentile(f, 99.5)
        # else:
        #     p0 = np.percentile(f, 0) #based on brain and liver data
        #     p1 = np.percentile(f, 100)
        p0 = np.percentile(f, lp)
        p1 = np.percentile(f, 100-lp)

        ranges.append((p0, p1))
    return ranges

# Help function for LUT
def rgb2hex(r,g,b):
    return "#{:02x}{:02x}{:02x}".format(r,g,b)



parser = argparse.ArgumentParser()
parser.add_argument("inputFiles", nARG='+', help="DEN files of the same ")
parser.add_argument("--labels", nARG='+', help="Title labels for perfusion maps, same order as inputFiles")

parser.add_argument("--name", help="Will be saved as PNG", default="output")
parser.add_argument("--dpi", type=int, default=300, help="Set DPI for output image.")

parser.add_argument("--mask-file", help="File to segment the organ")
parser.add_argument("--lower-percentile", '-lp', type=float, help="Recommended values for alpha 1.5 and beta are 0.5.")

parser.add_argument("--lut", help="LUT in csv file. Default used is rainbow.")

group = parser.add_mutually_exclusive_group()
group.add_argument("--ranges", nARG='+', help="Ranges for each file in format min-max, same order as inputFiles. Overwrites lower-percentile amd ttp-range-files.")
group.add_argument("--ttp-range-files", nARG='+', help="TTP DEN files without masking to determine TTP range, same order as inputFiles. Overwrites lower-percentile ranges.")

parser.add_argument("--aif", nargs='+', type=int, help='x y z')

parser.add_argument("--frames", nargs='+', type=int, help="All if not specified")
#parser.add_argument("--units", '-u', action='store_true')

ARG = parser.parse_args()


# Rainbow colormap
custom_cm=['#000000', '#4101a7', '#6504f7', '#6f08ff', '#6a03ff', '#6401ff', '#5e00ff', '#5b01ff', '#5800ff', '#5600ff', '#5200ff', '#4d00ff', '#4800ff', '#4200ff', '#3d00ff', '#3801ff', '#3402ff', '#2c01ff', '#2601ff', '#2002ff', '#1801fe', '#1301fe', '#0f01fe', '#0d02ff', '#0a05ff', '#0a08ff', '#0c0bff', '#090fff', '#0514ff', '#001bff', '#0020ff', '#0024ff', '#0028ff', '#002bff', '#002eff', '#0032ff', '#0036ff', '#003bff', '#003fff', '#0044ff', '#004cff', '#0053ff', '#0057ff', '#005cff', '#0061ff', '#0065ff', '#0069ff', '#006eff', '#0072ff', '#0075ff', '#007aff', '#007ffe', '#0087fa', '#008cf6', '#0091f4', '#0095f5', '#0098f6', '#009bf7', '#009ef6', '#00a3f6', '#00a8f6', '#00adf5', '#00b1f4', '#00b7f5', '#00bef6', '#00c2f4', '#00c7f4', '#00ccf4', '#00d1f5', '#00d4f4', '#00d8f5', '#00dcf6', '#00e0f8', '#00e4f9', '#00e8fa', '#00ecfa', '#00f0fa', '#00f4fc', '#00f6fe', '#00faff', '#00fcff', '#02feff', '#01feff', '#01feff', '#00fdff', '#00fdfc', '#00fdf9', '#00fdf5', '#00fef0', '#00ffeb', '#00ffe5', '#00fedf', '#00feda', '#00fed6', '#00ffd3', '#00ffcd', '#00ffc7', '#00ffc0', '#00ffbb', '#00ffb7', '#00ffb6', '#00ffb2', '#00ffaf', '#00ffaa', '#00ffa6', '#00ffa0', '#00ff9b', '#00ff98', '#00ff94', '#00ff8f', '#00ff8b', '#00ff86', '#00ff82', '#00ff7d', '#00ff79', '#00ff75', '#00ff72', '#00ff6f', '#00ff6c', '#00ff68', '#00ff64', '#00ff5f', '#00ff59', '#00ff54', '#00ff4d', '#02fe46', '#06fc3f', '#07fb3a', '#05fc35', '#02fd31', '#01fe2e', '#00ff2c', '#00ff2a', '#00ff27', '#00ff23', '#02ff20', '#0aff1f', '#12ff1d', '#14fe17', '#0bfc0b', '#04fe04', '#00ff01', '#00ff04', '#0bfe02', '#13fb01', '#19f700', '#1efc02', '#23fd03', '#26fd04', '#24f701', '#28f800', '#2efa01', '#36fe05', '#3bfc04', '#3ffb03', '#44fb03', '#49fb03', '#4efb04', '#53fb05', '#57fc05', '#5bfc05', '#5efc04', '#61fc04', '#65fd06', '#69fe06', '#6dff04', '#70fd02', '#74fc00', '#78fd00', '#7afe00', '#7dfe00', '#80ff00', '#85ff00', '#8aff00', '#90ff00', '#94ff00', '#98ff00', '#9cff00', '#a2ff00', '#a7fe00', '#abfe00', '#affd00', '#b3fe00', '#b8fe00', '#bcfe00', '#c0ff00', '#c3ff00', '#c4ff00', '#c9ff00', '#d2ff00', '#dbff00', '#e1ff00', '#e5ff00', '#e9ff00', '#eeff00', '#f2fb00', '#f5f700', '#f8f200', '#fced00', '#fee800', '#ffe300', '#ffdb00', '#ffd600', '#ffd200', '#ffcd00', '#ffc800', '#ffc400', '#ffc000', '#ffbd00', '#ffba00', '#ffb700', '#feb400', '#fcaf00', '#f8aa00', '#f7a602', '#f7a302', '#f69f00', '#f59b01', '#f59801', '#f59402', '#f69002', '#f78b03', '#f78704', '#f88306', '#f67b03', '#f67402', '#f77004', '#f66b03', '#f66602', '#f76303', '#f85f03', '#f85a02', '#f75601', '#f55300', '#f64e00', '#f84900', '#fd4200', '#ff3d00', '#ff3800', '#ff3400', '#ff2f00', '#ff2900', '#ff2200', '#ff1f00', '#ff1e00', '#ff1d00', '#ff1a00', '#ff1800', '#fe1300', '#f90900', '#f40f02', '#ee1202', '#e40700', '#ed0202', '#f00607', '#e81210']
#['#000000', '#020004', '#040008', '#06000c', '#080010', '#0a0014', '#0c0018', '#0e001c', '#100020', '#120024', '#140028', '#16002c', '#180030', '#1a0034', '#1c0038', '#1e003c', '#200040', '#220044', '#240048', '#26004c', '#280050', '#2a0054', '#2c0058', '#2e005c', '#300060', '#320064', '#340068', '#36006c', '#380070', '#3a0074', '#3c0078', '#3e007c', '#400080', '#3e0084', '#3c0088', '#3a008c', '#380090', '#360094', '#340098', '#32009c', '#3000a0', '#2e00a4', '#2c00a8', '#2a00ac', '#2800b0', '#2600b4', '#2400b8', '#2200bc', '#2000c0', '#1e00c4', '#1c00c8', '#1a00cc', '#1800d0', '#1600d4', '#1400d8', '#1200dc', '#1000e0', '#0e00e4', '#0c00e8', '#0a00ec', '#0800f0', '#0600f4', '#0400f8', '#0200fc', '#0000ff', '#0008f7', '#0010ef', '#0018e7', '#0020df', '#0028d7', '#0030cf', '#0038c7', '#0040bf', '#0048b7', '#0050af', '#0058a7', '#00609f', '#006897', '#00708f', '#007887', '#00807f', '#008877', '#00906f', '#009867', '#00a05f', '#00a857', '#00b04f', '#00b847', '#00c03f', '#00c837', '#00d02f', '#00d827', '#00e01f', '#00e817', '#00f00f', '#00f807', '#00ff00', '#04ff00', '#08ff00', '#0cff00', '#10ff00', '#14ff00', '#18ff00', '#1cff00', '#20ff00', '#24ff00', '#28ff00', '#2cff00', '#30ff00', '#34ff00', '#38ff00', '#3cff00', '#40ff00', '#44ff00', '#48ff00', '#4cff00', '#50ff00', '#54ff00', '#58ff00', '#5cff00', '#60ff00', '#64ff00', '#68ff00', '#6cff00', '#70ff00', '#74ff00', '#78ff00', '#7cff00', '#80ff00', '#84ff00', '#88ff00', '#8cff00', '#90ff00', '#94ff00', '#98ff00', '#9cff00', '#a0ff00', '#a4ff00', '#a8ff00', '#acff00', '#b0ff00', '#b4ff00', '#b8ff00', '#bcff00', '#c0ff00', '#c4ff00', '#c8ff00', '#ccff00', '#d0ff00', '#d4ff00', '#d8ff00', '#dcff00', '#e0ff00', '#e4ff00', '#e8ff00', '#ecff00', '#f0ff00', '#f4ff00', '#f8ff00', '#fcff00', '#ffff00', '#fffd00', '#fffb00', '#fff900', '#fff700', '#fff500', '#fff300', '#fff100', '#ffef00', '#ffed00', '#ffeb00', '#ffe900', '#ffe700', '#ffe500', '#ffe300', '#ffe100', '#ffdf00', '#ffdd00', '#ffdb00', '#ffd900', '#ffd700', '#ffd500', '#ffd300', '#ffd100', '#ffcf00', '#ffcd00', '#ffcb00', '#ffc900', '#ffc700', '#ffc500', '#ffc300', '#ffc100', '#ffc000', '#ffbd00', '#ffba00', '#ffb700', '#ffb400', '#ffb100', '#ffae00', '#ffab00', '#ffa800', '#ffa500', '#ffa200', '#ff9f00', '#ff9c00', '#ff9900', '#ff9600', '#ff9300', '#ff9000', '#ff8d00', '#ff8a00', '#ff8700', '#ff8400', '#ff8100', '#ff7e00', '#ff7b00', '#ff7800', '#ff7500', '#ff7200', '#ff6f00', '#ff6c00', '#ff6900', '#ff6600', '#ff6300', '#ff6000', '#ff5d00', '#ff5a00', '#ff5700', '#ff5400', '#ff5100', '#ff4e00', '#ff4b00', '#ff4800', '#ff4500', '#ff4200', '#ff3f00', '#ff3c00', '#ff3900', '#ff3600', '#ff3300', '#ff3000', '#ff2d00', '#ff2a00', '#ff2700', '#ff2400', '#ff2100', '#ff1e00', '#ff1b00', '#ff1800', '#ff1500', '#ff1200', '#ff0f00', '#ff0c00', '#ff0900', '#ff0600', '#ff0300']

AIF=ARG.aif
# What frames will be processed
inf = DEN.readHeader(ARG.inputFiles[0])
dimspec = inf["dimspec"]
dimx = dimspec[0]
dimy = dimspec[1]
dimz = dimspec[2]

frames=[]
if (ARG.frames):
    frames = ARG.frames
else:
    frames = range(0, zdim)


# Set labels
labels = []
for i in ARG.inputFiles:
    labels.append('')
if (ARG.labels):
    if len(ARG.labels)!=len(ARG.inputFiles):
        print("Not all input files have labels")
    else:
        labels=[]
        for l in ARG.labels:
            labels.append((l))


# Set ranges
ranges = computeRange(ARG.inputFiles, ARG.lower_percentile)
if (ARG.mask_file):
    ranges = computeRange(ARG.inputFiles, ARG.lower_percentile, ARG.mask_file)
if (ARG.ttp_range_files):
    if (not ARG.aif):
        print('AIF is required')
        sys.exit()
    else:
        if AIF[0]>dimx or AIF[1]>dimy or AIF[2]>dimz:
            print('The coordinates of the AIF (%d, %d, %d) are greater than the dimensions of the volume (%d, %d, %d)'%(AIF[0], AIF[1], AIF[2], dimx, dimy, dimz))
            sys.exit()
        print("Using original volumes to determine TTP ranges")
        ranges=[]
        for r in ARG.ttp_range_files:
            a = subprocess.check_output(['dentk-value', str(ARG.aif[0]), str(ARG.aif[1]), str(ARG.aif[2]), r])
            ranges.append((float(a)-10, float(a)+20))
if (ARG.ranges):
    ranges=[]
    for r in ARG.ranges:
        ranges.append(list(map(float, r.split('-'))))

#Set name prefix
name= ARG.name

#Create Rainbow LUT
if (ARG.lut):
    custom_cm=[]
    with open(ARG.lut) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        line_count = 0
        for row in csv_reader:
            custom_cm.append(rgb2hex(int(row[0]), int(row[1]), int(row[2])))

cmap = colors.ListedColormap(custom_cm)

# mask = np.ones((ydim, xdim))

#Set label
units_label = ''
im = ARG.inputFiles[0]
if 'BF' in im:
    units_label = r'$\frac{mL}{100g \cdot\dot min}$'
if 'BV' in im:
    units_label = r'$\frac{mL}{100\dotg}$'
if 'MTT' in im:
    units_label = r'$s$'
if 'TTP' in im:
    units_label = r'$s$'

#Create images
for f in frames:

    fig = pyplot.figure()
    grid = ImageGrid(fig, 111,
                    nrows_ncols=(1, len(ARG.inputFiles)))

    for ax, im, l, r in zip(grid, ARG.inputFiles, labels, ranges):

        img = DEN.getFrame(im, f)

        if (ARG.mask_file):
            mask = np.ones((ydim, xdim))
            mask = DEN.getFrame(ARG.mask_file, f)
            img = np.multiply(img, mask)

        img=ax.imshow(img, cmap='Greys_r', vmin=r[0], vmax=r[1])

        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        #bname=''
        #if os.path.splitext(os.path.basename(im))[0][-1] == 'A' or os.path.splitext(os.path.basename(im))[0][-1]=='B':
        bname = os.path.splitext(os.path.basename(im))[0]
        #else:
        #    bname = os.path.splitext(os.path.basename(im))[0]

        ax.text(1,1,  bname + '\n'+ l, verticalalignment='top', horizontalalignment='left', color='white', fontsize=6)

        cbaxes = inset_axes(ax, width="3%", height="30%", loc=1)
        cb = fig.colorbar(img, cax=cbaxes, orientation='vertical', ticks=[math.ceil(r[0]), math.floor(r[1])])
        cb.ax.tick_params(labelsize=6, size=0, pad=0.1, labelcolor='white', labelleft='on', labelright='off')
        cb.ax.set_xlabel(units_label, fontsize=4)
        cb.ax.xaxis.label.set_color('white')
    pyplot.savefig("%s%03d.png"%(name, f), format="png", dpi=ARG.dpi, bbox_inches='tight', pad_inches = 0)
    pyplot.close(fig)
