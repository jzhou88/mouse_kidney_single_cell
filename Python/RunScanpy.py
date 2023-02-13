#!/usr/bin/python

import argparse
import datetime
import inspect
import itertools
import math
import matplotlib
import numpy
import os
import pandas
import phate
import re
import scanorama
import scanpy
import scipy
import scprep
import seaborn
import string
import struct
import sys
import time

import Common
import ScanpyAnalysis

## 2020-11-04-ScanpyHarmony14
## resolution 0.9
## annotated on 2021-04-14
cluster2annotation = {'0': 'DCT', \
					  '1': 'PT', \
					  '2': 'PT', \
					  '3': 'PT', \
					  '4': 'ALOH', \
					  '5': 'T lymph', \
					  '6': 'PT', \
					  '7': 'Granul', \
					  '8': 'Macro', \
					  '9': 'Endo', \
					  '10': 'CD PC', \
					  '11': 'PT', \
					  '12': 'PT', \
					  '13': 'PT', \
					  '14': 'Mono', \
					  '15': 'PT', \
					  '16': 'PT', \
					  '17': 'NK', \
					  '18': 'CD IC', \
					  '19': 'PT', \
					  '20': 'B lymph', \
					  '21': 'PT', \
					  '22': 'Novel2', \
					  '23': 'PT', \
					  '24': 'DLOH', \
					  '25': 'DC & pDC', \
					  '26': 'PT', \
					  '27': 'PT', \
					  '28': 'Podo', \
					  '29': 'PT', \
					  '30': 'PT', \
					  '31': 'PT & DCT'}


def runScanpy(args):
	try:
		integ = args.i
		dirOut = os.path.abspath(args.o)
		resol = args.r
		step = args.s
		sufName = args.suf
		markers = args.m
		npcs = args.pc

		scanpy.settings.verbosity = 3
		scanpy.logging.print_header()
		scanpy.settings.set_figure_params(dpi_save=150, vector_friendly=False, fontsize=14, facecolor='white')

		if (1 == step):
			print('\nRunning ScanpyAnalysis.cleanLoom(%s)...' % (dirOut))
			ScanpyAnalysis.cleanLoom(dirOut)

		elif (2 == step):
			print('\nRunning ScanpyAnalysis.preprocess(%s, %s)...' % (dirOut, sufName))
			ScanpyAnalysis.preprocess(dirOut, sufName)

		elif (3 == step):
			print('\nRunning ScanpyAnalysis.pca(%s, %s)...' % (dirOut, npcs))
			ScanpyAnalysis.pca(dirOut, npcs)

		elif (4 == step):
			print('\nRunning ScanpyAnalysis.integration(%s, %s, %s)...' % (dirOut, integ, npcs))
			ScanpyAnalysis.integration(dirOut, integ, npcs)

		elif (5 == step):
			print('\nRunning ScanpyAnalysis.neighbor(%s, %s, %s)...' % (dirOut, integ, npcs))
			ScanpyAnalysis.neighbor(dirOut, integ, npcs)

		elif (6 == step):
			print('\nRunning ScanpyAnalysis.clustering(%s, %s, %s)...' % (dirOut, resol, markers))
			ScanpyAnalysis.clustering(dirOut, resol, markers)

		elif (7 == step):
			print('\nRunning ScanpyAnalysis.annotation(%s, %s, )...' % (dirOut, resol))
			ScanpyAnalysis.annotation(dirOut, resol, cluster2annotation)
			
		else:
			Common.error('could not identify Step \'%s\'' % (step))

		return 0
	except Exception as e:
		Common.exception(e, 'could not run scanpy')


def main():
	parser = argparse.ArgumentParser(description = 'run scanpy', formatter_class = argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', default = 'harmony', help = 'optional: specify which integration approach to run')
	parser.add_argument('-o', required = True, help = 'specify where to save the outputs')
	parser.add_argument('-r', type = float, default = 0.0, help = 'optional: specify resolution to run the clustering')
	parser.add_argument('-s', required = True, type = int, help = 'select which step to run')
	parser.add_argument('--suf', default = '', help = 'optional: specify the suffix of file names')
	parser.add_argument('--pc', type = int, default = 50, help = 'optional: specify the number of PCs to use')
	parser.add_argument('-m', type = int, default = 0, help = 'optional: specify what marker genes to plot after the clustering')
	args = parser.parse_args()
	runScanpy(args)
	return 0


if __name__ == '__main__':
	sys.exit(main())

