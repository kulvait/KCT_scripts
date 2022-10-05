#!/usr/bin/env python

"""
@Author: Hana Haseljic, Vojtech Kulvait
@Date: 2020-2022
@License: GPL3 
"""

#Example
#perfusion-mosaic.py CBF.den CBV.den MTT.den TTP.den --mask-file ALPHA_INSIDE.den -lp 2 --name output_ --frames 2 14 21 --lut Rainbow.csv
#perfusion-mosaic.py CBFB.den CBVB.den MTTB.den TTPB.den --ttp-range-file TTP.den -mtt 0-60 --name output --frames 15

#
#Examples for Erlangen data
#

#python perfusion-mosaic-volumewise.py /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_3.50_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/CBFB.den /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_3.50_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/CBVB.den /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_3.50_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/MTTB.den /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_3.50_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/TTPB.den  --aif 274 237 2 --frames 13

#perfusion-mosaic-volumewise.py /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_3.50_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/CBFB.den /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_3.50_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/CBVB.den /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_3.50_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/MTTB.den /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_3.50_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/TTPB.den   --aif 274 237 2 --frames 13 -lp 1.5

#python perfusion-mosaic-volumewise.py /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_3.50_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/CBFA.den /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_3.50_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/CBVA.den /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_3.50_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/MTTA.den /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_3.50_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/TTPA.den --mask-file /mnt/md1337/kulvait/HanaViz2/stroke_play_midpo/Pat1A/06GaussViz_3.50_1.3.12.2.1107.5.8.15.101021.30000020080108270715200001298AGG/BETA.den  --aif 274 237 2 --frames 13 -lp 1

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


sys.stderr.write("perfusion-mosaic-volumewise.py called from %s\n"%(os.path.abspath(os.getcwd())))

# hc: Rainbow colormap
custom_cm=['#000000', '#020004', '#040008', '#06000c', '#080010', '#0a0014', '#0c0018', '#0e001c', '#100020', '#120024', '#140028', '#16002c', '#180030', '#1a0034', '#1c0038', '#1e003c', '#200040', '#220044', '#240048', '#26004c', '#280050', '#2a0054', '#2c0058', '#2e005c', '#300060', '#320064', '#340068', '#36006c', '#380070', '#3a0074', '#3c0078', '#3e007c', '#400080', '#3e0084', '#3c0088', '#3a008c', '#380090', '#360094', '#340098', '#32009c', '#3000a0', '#2e00a4', '#2c00a8', '#2a00ac', '#2800b0', '#2600b4', '#2400b8', '#2200bc', '#2000c0', '#1e00c4', '#1c00c8', '#1a00cc', '#1800d0', '#1600d4', '#1400d8', '#1200dc', '#1000e0', '#0e00e4', '#0c00e8', '#0a00ec', '#0800f0', '#0600f4', '#0400f8', '#0200fc', '#0000ff', '#0008f7', '#0010ef', '#0018e7', '#0020df', '#0028d7', '#0030cf', '#0038c7', '#0040bf', '#0048b7', '#0050af', '#0058a7', '#00609f', '#006897', '#00708f', '#007887', '#00807f', '#008877', '#00906f', '#009867', '#00a05f', '#00a857', '#00b04f', '#00b847', '#00c03f', '#00c837', '#00d02f', '#00d827', '#00e01f', '#00e817', '#00f00f', '#00f807', '#00ff00', '#04ff00', '#08ff00', '#0cff00', '#10ff00', '#14ff00', '#18ff00', '#1cff00', '#20ff00', '#24ff00', '#28ff00', '#2cff00', '#30ff00', '#34ff00', '#38ff00', '#3cff00', '#40ff00', '#44ff00', '#48ff00', '#4cff00', '#50ff00', '#54ff00', '#58ff00', '#5cff00', '#60ff00', '#64ff00', '#68ff00', '#6cff00', '#70ff00', '#74ff00', '#78ff00', '#7cff00', '#80ff00', '#84ff00', '#88ff00', '#8cff00', '#90ff00', '#94ff00', '#98ff00', '#9cff00', '#a0ff00', '#a4ff00', '#a8ff00', '#acff00', '#b0ff00', '#b4ff00', '#b8ff00', '#bcff00', '#c0ff00', '#c4ff00', '#c8ff00', '#ccff00', '#d0ff00', '#d4ff00', '#d8ff00', '#dcff00', '#e0ff00', '#e4ff00', '#e8ff00', '#ecff00', '#f0ff00', '#f4ff00', '#f8ff00', '#fcff00', '#ffff00', '#fffd00', '#fffb00', '#fff900', '#fff700', '#fff500', '#fff300', '#fff100', '#ffef00', '#ffed00', '#ffeb00', '#ffe900', '#ffe700', '#ffe500', '#ffe300', '#ffe100', '#ffdf00', '#ffdd00', '#ffdb00', '#ffd900', '#ffd700', '#ffd500', '#ffd300', '#ffd100', '#ffcf00', '#ffcd00', '#ffcb00', '#ffc900', '#ffc700', '#ffc500', '#ffc300', '#ffc100', '#ffc000', '#ffbd00', '#ffba00', '#ffb700', '#ffb400', '#ffb100', '#ffae00', '#ffab00', '#ffa800', '#ffa500', '#ffa200', '#ff9f00', '#ff9c00', '#ff9900', '#ff9600', '#ff9300', '#ff9000', '#ff8d00', '#ff8a00', '#ff8700', '#ff8400', '#ff8100', '#ff7e00', '#ff7b00', '#ff7800', '#ff7500', '#ff7200', '#ff6f00', '#ff6c00', '#ff6900', '#ff6600', '#ff6300', '#ff6000', '#ff5d00', '#ff5a00', '#ff5700', '#ff5400', '#ff5100', '#ff4e00', '#ff4b00', '#ff4800', '#ff4500', '#ff4200', '#ff3f00', '#ff3c00', '#ff3900', '#ff3600', '#ff3300', '#ff3000', '#ff2d00', '#ff2a00', '#ff2700', '#ff2400', '#ff2100', '#ff1e00', '#ff1b00', '#ff1800', '#ff1500', '#ff1200', '#ff0f00', '#ff0c00', '#ff0900', '#ff0600', '#ff0300']

siemens_cb_cm = ['#000000', '#0b000b', '#160016', '#210021', '#2b002b', '#2e002e', '#320032', '#360036', '#3a003a', '#3d003d', '#410041', '#450045', '#470049', '#47004d', '#470050', '#470054', '#440055', '#400055', '#3c0055', '#380055', '#380459', '#38075c', '#380b60', '#370e63', '#340e63', '#300e63', '#2c0e63', '#2a0f63', '#2a1363', '#2a1763', '#2a1a63', '#2a1e65', '#2a2269', '#2a256c', '#2a2970', '#2a2a74', '#2a2a78', '#2a2a7c', '#2a2a80', '#2a2a83', '#2a2a87', '#2a2a8b', '#2a2a8e', '#2a2a92', '#2a2a96', '#2a2a9a', '#2a2a9d', '#2a2aa1', '#2a2aa5', '#2a2aa8', '#2a2aac', '#2a2ab0', '#2a2ab3', '#2a2ab7', '#2a2abb', '#2a2abf', '#2a2ac3', '#2a2ac7', '#2a2aca', '#2a2ace', '#2a2ad2', '#2a2ad5', '#2a2ad9', '#2a2add', '#2a2ae0', '#2a2ae4', '#2a2ae8', '#2a2aeb', '#2a2aef', '#2a2af3', '#2a2af6', '#2a2afa', '#2a2afe', '#2a2dfa', '#2a30f3', '#2a34eb', '#2a38e4', '#2a38e0', '#2a38dc', '#2a38d8', '#2a38d5', '#2a3cd5', '#2a40d5', '#2a44d5', '#2a47d4', '#2a47d0', '#2a47cd', '#2a47c9', '#2c49c7', '#2f4cc7', '#3350c7', '#3754c7', '#3855c4', '#3855c0', '#3855bd', '#3855b9', '#3858b5', '#385cb1', '#385fae', '#3863aa', '#3867a6', '#386ba2', '#386e9f', '#38719b', '#387197', '#387194', '#387190', '#3a738e', '#3e778e', '#427b8e', '#467f8e', '#47808c', '#478088', '#478084', '#478081', '#47837d', '#478779', '#478a75', '#478e71', '#47926d', '#47956a', '#479966', '#489c62', '#4b9c5f', '#4f9c5b', '#539c57', '#559d55', '#55a155', '#55a555', '#55a955', '#55aa53', '#55aa4f', '#55aa4b', '#55aa48', '#55ad47', '#55b147', '#55b447', '#55b847', '#55b843', '#55b83f', '#55b83b', '#55b937', '#55bd34', '#55c130', '#55c52c', '#56c729', '#5ac725', '#5ec721', '#61c71e', '#63c91c', '#63cd1c', '#63d01c', '#63d41c', '#63d519', '#63d516', '#63d512', '#63d50e', '#63d90e', '#63dc0e', '#63e00e', '#63e30d', '#63e30a', '#63e306', '#63e302', '#63e400', '#63e800', '#63ec00', '#63ef00', '#65f300', '#69f700', '#6cfa00', '#70fe00', '#74ff00', '#78ff00', '#7cff00', '#80ff00', '#83ff00', '#87ff00', '#8bff00', '#8eff00', '#92ff00', '#96ff00', '#99ff00', '#9dff00', '#a1ff00', '#a5ff00', '#a8ff00', '#acff00', '#b0ff00', '#b3ff00', '#b7ff00', '#bbff00', '#bfff00', '#c3ff00', '#c7ff00', '#caff00', '#ceff00', '#d2ff00', '#d5ff00', '#d9ff00', '#ddff00', '#e0ff00', '#e4ff00', '#e8ff00', '#ebff00', '#efff00', '#f3ff00', '#f6ff00', '#faff00', '#feff00', '#fffd00', '#fff900', '#fff500', '#fff100', '#ffee00', '#ffea00', '#ffe600', '#ffe300', '#fbdb00', '#f7d400', '#f4cc00', '#f1c600', '#f1c200', '#f1be00', '#f1ba00', '#f1b500', '#f1ad00', '#f1a600', '#f19f00', '#f19a00', '#f19600', '#f19200', '#f18f00', '#f18800', '#f18000', '#f17800', '#f17100', '#ed6900', '#ea6200', '#e65b00', '#e35400', '#e35000', '#e34d00', '#e34900', '#e34400', '#e33c00', '#e33500', '#e32d00', '#e32800', '#e32400', '#e32000', '#e31d00', '#e31900', '#e31500', '#e31200', '#e30e00']

siemens_tt_cm = ['#000000', '#380038', '#710071', '#a900a9', '#b500b8', '#b100b8', '#ac00b8', '#aa00ba', '#aa00bf', '#aa00c3', '#a900c7', '#a500c7', '#a100c7', '#9c00c7', '#9800c7', '#9400c7', '#8f00c7', '#8b00ca', '#8700ce', '#8300d2', '#7e00d5', '#7a00d5', '#7500d5', '#7100d5', '#6c00da', '#6800de', '#6400e2', '#6304e3', '#6308e3', '#630ce3', '#610ee3', '#5c0ee3', '#580ee3', '#540ee4', '#4f0ee9', '#4b0eed', '#470ef1', '#4712f1', '#4717f1', '#471bf1', '#431cf1', '#3f1cf1', '#3a1cf1', '#361cf3', '#321cf7', '#2d1cfc', '#2a1dff', '#2a21ff', '#2a26ff', '#2a2aff', '#262afb', '#222af7', '#1d2af2', '#1c2df1', '#1c31f1', '#1c36f1', '#1c38ef', '#1c38eb', '#1c38e7', '#1c38e2', '#1c38de', '#1c38da', '#1c38d5', '#1c3cd5', '#1c41d5', '#1c45d5', '#1c47d2', '#1c47ce', '#1c47ca', '#1c48c7', '#1c4dc7', '#1c51c7', '#1c55c7', '#1755c2', '#1355bd', '#0f55b9', '#0e58b5', '#0e5db0', '#0e61ac', '#0e65a8', '#0e6aa3', '#0e6e9f', '#0e719b', '#0e7197', '#0e7192', '#0e718e', '#0e768e', '#0e7a8e', '#0e7f8e', '#0b808b', '#078087', '#028082', '#00827e', '#008679', '#008b75', '#008f70', '#00936c', '#009768', '#009c63', '#009c5f', '#009c5b', '#009c56', '#009f55', '#00a355', '#00a755', '#00aa53', '#00aa4f', '#00aa4b', '#00aa47', '#00af47', '#00b347', '#00b747', '#00b843', '#00b83f', '#00b83a', '#00bb36', '#00bf31', '#00c42d', '#00c829', '#00cd24', '#00d120', '#00d51c', '#00d518', '#00d513', '#00d50f', '#00d80e', '#00dd0e', '#00e10e', '#00e30c', '#00e308', '#00e303', '#00e400', '#00e800', '#00ed00', '#00f100', '#00f500', '#00f900', '#00fe00', '#03ff00', '#07ff00', '#0cff00', '#10ff00', '#14ff00', '#18ff00', '#1dff00', '#21ff00', '#25ff00', '#2aff00', '#2eff00', '#32ff00', '#36ff00', '#3bff00', '#3fff00', '#44ff00', '#48ff00', '#4dff00', '#51ff00', '#55ff00', '#5aff00', '#5eff00', '#62ff00', '#66ff00', '#6bff00', '#6fff00', '#73ff00', '#78ff00', '#7dff00', '#81ff00', '#85ff00', '#8aff00', '#8eff00', '#92ff00', '#97ff00', '#9bff00', '#9fff00', '#a3ff00', '#a8ff00', '#acff00', '#b0ff00', '#b5ff00', '#b9ff00', '#bdff00', '#c2ff00', '#c7ff00', '#cbff00', '#cfff00', '#d4ff00', '#d8ff00', '#dcff00', '#e0ff00', '#e5ff00', '#e9ff00', '#edff00', '#f1ff00', '#f6ff00', '#faff00', '#feff00', '#fffb00', '#fff700', '#fff300', '#ffef00', '#ffea00', '#ffe600', '#ffe200', '#ffdd00', '#ffd900', '#ffd500', '#ffd100', '#ffcc00', '#ffc800', '#ffc300', '#ffbf00', '#ffba00', '#ffb600', '#ffb200', '#ffad00', '#ffa900', '#ffa500', '#ffa000', '#ff9c00', '#ff9800', '#ff9400', '#ff8f00', '#ff8b00', '#ff8700', '#ff8200', '#ff7e00', '#ff7900', '#ff7500', '#ff7000', '#ff6c00', '#ff6800', '#ff6300', '#ff5f00', '#ff5b00', '#ff5700', '#ff5200', '#ff4e00', '#ff4a00', '#ff4500', '#ff4100', '#ff3c00', '#ff3800', '#ff3300', '#ff2f00', '#ff2b00', '#ff2700', '#ff2200', '#ff1e00', '#ff1a00', '#ff1500', '#ff1100', '#ff0d00', '#ff0900', '#ff0400', '#ff0000']

asist_cm = ['#000000', '#040000', '#080003', '#0c0006', '#100009', '#14000c', '#18000f', '#1c0012', '#200015', '#240018', '#28001b', '#2c001e', '#300021', '#340024', '#380027', '#3c002a', '#40002d', '#440030', '#480033', '#4c0036', '#500039', '#4f003c', '#4e003f', '#4d0042', '#4c0045', '#4b0048', '#4a004b', '#49004e', '#480051', '#470054', '#460057', '#45005a', '#44005d', '#430060', '#420063', '#410066', '#400069', '#3f006c', '#3e006f', '#3d0072', '#3c0075', '#3b0078', '#3a007b', '#39007e', '#380081', '#370084', '#360387', '#35068a', '#34098d', '#330c90', '#320f93', '#311296', '#301599', '#2f189c', '#2e1b9f', '#2d1ea2', '#2c21a5', '#2b24a8', '#2a27ab', '#292aae', '#282db1', '#2730b4', '#2633b7', '#2536ba', '#2439bd', '#233cc0', '#223fc3', '#2142c6', '#2045c9', '#1f48cc', '#1e4bcf', '#1d4ed2', '#1c51d5', '#1b54d8', '#1a57db', '#195ade', '#185de1', '#1760e4', '#1663e7', '#1566ea', '#1469ed', '#136cf0', '#126ff3', '#1172f6', '#1075f9', '#0f78fc', '#0e7bff', '#0d7efa', '#0c81f5', '#0b84f0', '#0a87eb', '#098ae6', '#088de1', '#0790dc', '#0693d7', '#0596d2', '#0499cd', '#039cc8', '#029fc3', '#01a2be', '#00a5b9', '#00a8b4', '#00abaf', '#00aeaa', '#00b1a5', '#00b4a0', '#00b79b', '#00ba96', '#00bd91', '#00c08c', '#00c387', '#00c682', '#00c97d', '#00cc78', '#00cf73', '#00d26e', '#00d569', '#00d864', '#00db5f', '#00de5a', '#00e155', '#00e450', '#00e74b', '#00ea46', '#00ed41', '#00f03c', '#00f337', '#00f632', '#00f92d', '#04fc28', '#08ff23', '#0cff1e', '#10ff19', '#14ff14', '#18ff0f', '#1cff0a', '#20ff05', '#24ff00', '#28ff00', '#2cff00', '#30ff00', '#34ff00', '#38ff00', '#3cff00', '#40ff00', '#44ff00', '#48ff00', '#4cff00', '#50ff00', '#54ff00', '#58ff00', '#5cff00', '#60ff00', '#64ff00', '#68ff00', '#6cff00', '#70ff00', '#74ff00', '#78ff00', '#7cff00', '#80ff00', '#84ff00', '#88ff00', '#8cff00', '#90ff00', '#94ff00', '#98ff00', '#9cff00', '#a0ff00', '#a4ff00', '#a8ff00', '#acff00', '#b0ff00', '#b4ff00', '#b8ff00', '#bcff00', '#c0ff00', '#c4ff00', '#c8ff00', '#ccff00', '#d0ff00', '#d4ff00', '#d8ff00', '#dcff00', '#e0ff00', '#e4ff00', '#e8ff00', '#ecff00', '#f0ff00', '#f4ff00', '#f8ff00', '#fcff00', '#fffc00', '#fff800', '#fff400', '#fff000', '#ffec00', '#ffe800', '#ffe400', '#ffe000', '#ffdc00', '#ffd800', '#ffd400', '#ffd000', '#ffcc00', '#ffc800', '#ffc400', '#ffc000', '#ffbc00', '#ffb800', '#ffb400', '#ffb000', '#ffac00', '#ffa800', '#ffa400', '#ffa000', '#ff9c00', '#ff9800', '#ff9400', '#ff9000', '#ff8c00', '#ff8800', '#ff8400', '#ff8000', '#ff7c00', '#ff7800', '#ff7400', '#ff7000', '#ff6c00', '#ff6800', '#ff6400', '#ff6000', '#ff5c00', '#ff5800', '#ff5400', '#ff5000', '#ff4c00', '#ff4800', '#ff4400', '#ff4000', '#ff3c00', '#ff3800', '#ff3400', '#ff3000', '#ff2c00', '#ff2800', '#ff2400', '#ff2000', '#ff1c00', '#ff1800', '#ff1400', '#ff1000', '#ff0c00', '#ff0800', '#ff0400', '#ff0000']




parser = argparse.ArgumentParser()
parser.add_argument("inputFiles", nargs='+', help="BF, BV, MTT, TTP (in this order)")
parser.add_argument("--name", help="Will be saved as PNG")
parser.add_argument("--dpi", type=int, default=300, help="Set DPI for output image.")
parser.add_argument("--format", default='png', help="Set format for the output image.")

parser.add_argument("--mask-file", help="File to segment the organ", default=argparse.SUPPRESS)
parser.add_argument("--lower-percentile", '-lp', type=float, default=0.0, help="Recommended values for alpha 1.5 and beta 0.5.")

parser.add_argument("--lut", help="LUT in csv file. Default used is rainbow.")
parser.add_argument("--siemens", action='store_true', help="Use Siemens colormaps. Default is rainbow.")
parser.add_argument("--asist", action='store_true', help="Use Kohsuke Kudo colormap. Default is rainbow.")

parser.add_argument("--cbf-range", '-cbf', help="CBF range in format min-max. Overwrites lower-percentile ranges.")
parser.add_argument("--cbv-range", '-cbv', help="CBV range in format min-max. Overwrites lower-percentile ranges.")
parser.add_argument("--mtt-range", '-mtt', help="MTT range in format min-max. Overwrites lower-percentile ranges.")
parser.add_argument("--ttp-range", '-ttp', help="TTP range in format min-max. Overwrites lower-percentile and ttp-range-file")

parser.add_argument("--ttp-range-file", help="TTP file without masking to determine TTP range. Overwrites lower-percentile ranges.")

parser.add_argument("--layout", help="Image grid layout in format rows,column. Options are (1,4), (4,1) and (2,2). Default is (2,2)." )

parser.add_argument("--aif", nargs='+', type=int, help='x y z')

parser.add_argument("--frames", nargs='+', type=int, help="All if not specified")
#parser.add_argument("--units", '-u', action='store_true')

ARG = parser.parse_args()

# Compute window based on the volume
def computeRange(inputFilesArrays, lp, mask_file = None):
    ranges=[]
    for array in inputFilesArrays:
        p0=0
        p1=0
        if mask_file is not None:
            array = np.multiply(array, mask_file)
        # if os.path.splitext(os.path.basename(file))[0][-1] == 'A':
        #     p0=np.percentile(f, 1.5) #based on brain data
        #     p1=np.percentile(f, 98.5)
        # elif os.path.splitext(os.path.basename(file))[0][-1] == 'B':
        #     p0 = np.percentile(f, 0.5) #based on brain and liver data
        #     p1 = np.percentile(f, 99.5)
        # else:
        #     p0 = np.percentile(f, 0) #based on brain and liver data
        #     p1 = np.percentile(f, 100)
        p0 = np.percentile(array, lp)
        p1 = np.percentile(array, 100-lp)
        ranges.append((p0, p1))
    return ranges

# Help function for LUT
def rgb2hex(r,g,b):
    return "#{:02x}{:02x}{:02x}".format(r,g,b)

AIF=ARG.aif
# What frames will be processed
inf = DEN.readHeader(ARG.inputFiles[0])
dimspec = inf["dimspec"]
dimx = dimspec[0]
dimy = dimspec[1]
dimz = dimspec[2]

inputFiles = []
mask_file = None
if "mask_file" in ARG:
	mask_file = DEN.getNumpyArray(ARG.mask_file)
for f in ARG.inputFiles:
	inputFiles.append(DEN.getNumpyArray(f))


frames=[]
if (ARG.frames):
    frames = ARG.frames
else:
    frames = range(0, dimz)

format = ARG.format
formats = ['eps', 'jpeg', 'jpg', 'pdf', 'pgf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff']
if format not in formats:
    print("Format is not supported. Saving PNG.")
    format = 'png'

#Set ranges
ranges = computeRange(inputFiles, ARG.lower_percentile, mask_file)

if (ARG.cbf_range):
    ranges[0]=tuple(map(float, ARG.cbf_range.split('-')))
if (ARG.cbv_range):
    ranges[1]=tuple(map(float, ARG.cbv_range.split('-')))
if (ARG.mtt_range):
    ranges[2]=tuple(map(float, ARG.mtt_range.split('-')))
if (ARG.ttp_range_file):
    if (not ARG.aif):
        print('AIF is required')
        sys.exit()
    else:
        if AIF[0]>dimx or AIF[1]>dimy or AIF[2]>dimz:
            print('The coordinates of the AIF (%d, %d, %d) are greater than the dimensions of the volume (%d, %d, %d)'%(AIF[0], AIF[1], AIF[2], dimx, dimy, dimz))
            sys.exit()
        print("Using original volume to determine TTP range")
        a = subprocess.check_output(['dentk-value', str(ARG.aif[0]), str(ARG.aif[1]), str(ARG.aif[2]), ARG.ttp_range_file])
        ranges[3]=(float(a)-10, float(a)+20)
if (ARG.ttp_range):
    ranges[3]=tuple(map(float, ARG.ttp_range.split('-')))
#Set name prefix
name=""
if (ARG.name):
    name=ARG.name
else:
    name="output"

#Create Rainbow LUT
if (ARG.lut):
    custom_cm=[]
    with open(ARG.lut) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        line_count = 0
        for row in csv_reader:
            custom_cm.append(rgb2hex(int(row[0]), int(row[1]), int(row[2])))

cmaps=[]
if ARG.siemens:
    cmap_cb = colors.ListedColormap(siemens_cb_cm)
    cmap_tt = colors.ListedColormap(siemens_tt_cm)
    cmaps = [cmap_cb, cmap_cb, cmap_tt, cmap_tt]
elif ARG.asist:
    cmap = colors.ListedColormap(asist_cm)
    cmaps = [cmap, cmap, cmap, cmap]
else:
    cmap = colors.ListedColormap(custom_cm)
    cmaps = [cmap, cmap, cmap, cmap]

# if ARG.mask_file and (not ARG.lower_percentile):
#    print('Window has to be setted.')
#    sys.exit()
layout=[2, 2]


if (ARG.layout):
    tmp=list(map(int, ARG.layout.split(',')))
    if tmp[0]*tmp[1]!=4:
        print("Layout (2,2) will be used.")
    else:
        layout=tmp

#Set units_label
units_label = [(r'$\frac{mL}{100g \cdot min}$'), (r'$\frac{mL}{100g}$'), (r'$s$'), (r'$s$')]


#Create images
for f in frames:

    fig = pyplot.figure()

    if layout[0]!=2:
        fig.set_size_inches(layout[0]*15/2.54, layout[1]*15/2.54)

    grid = ImageGrid(fig, 111,
                    nrows_ncols=(layout[0], layout[1]))


    for ax, im, r, l, c, stringFile in zip(grid, inputFiles, ranges, units_label, cmaps, ARG.inputFiles):

        img = im[f]

        if mask_file is not None:
            mask = mask_file[f]
            img = np.multiply(img, mask)
        img=ax.imshow(img, cmap=c, vmin=r[0], vmax=r[1])

        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        #bname=''
        #if os.path.splitext(os.path.basename(im))[0][-1] == 'A' or os.path.splitext(os.path.basename(im))[0][-1]=='B':
        bname = os.path.splitext(os.path.basename(stringFile))[0][:3]
        #else:
        #    bname = os.path.splitext(os.path.basename(im))[0]
        ax.text(1,1, bname, verticalalignment='top', horizontalalignment='left', color='white', fontsize=5)
        ax.text(1,40, l, verticalalignment='top', horizontalalignment='left', color='white', fontsize=4)

        cbaxes = inset_axes(ax, width="3%", height="30%", loc=1)
        cb = fig.colorbar(img, cax=cbaxes, orientation='vertical', ticks=[math.ceil(r[0]), math.floor(r[1])])
        cb.ax.tick_params(labelsize=6, size=0, pad=0.1, labelcolor='white', labelleft=True, labelright=False)
        #if (ARG.units):
        #cb.ax.set_xlabel(l, fontsize=4.2)
        #cb.ax.xaxis.label.set_color('white')
        #if len(l)>20:
        #    cb.ax.xaxis.label.set_x(-0.2)
    pyplot.savefig("%s%03d.%s"%(name, f, format), format=format, dpi=ARG.dpi, bbox_inches='tight', pad_inches = 0)
    pyplot.close(fig)
