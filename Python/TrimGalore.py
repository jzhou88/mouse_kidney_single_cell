#!/usr/bin/python2

import argparse
import datetime
import inspect
import itertools
import math
import numpy
import os
import re
import scipy
import string
import struct
import sys
import time

import Common


def trimGalore(args):
	try:
		dirWork = args.w
		dirIn = args.i
		dirOut = args.o
		pairedEnd = args.p
		submit = args.s
		cores = args.c
		mem = args.m
		denovo = args.d

		dirWork = os.path.abspath(dirWork)
		dirIn = os.path.abspath(dirIn)
		dirOut = os.path.abspath(dirOut)
		Common.mkDir(dirWork)
		Common.dirExist(dirIn)
		Common.mkDir(dirOut)

		dirInEnd = None
		dirOutEnd = None
		cmdComm = cmdComm = '%s --fastqc' % (Common.trimGaloreCurr)
		if pairedEnd:
			dirInEnd = os.path.join(dirIn, 'paired-end')
			dirOutEnd = os.path.join(dirOut, 'paired-end')
			cmdComm = '%s --retain_unpaired --paired' % (cmdComm)
		else:
			dirInEnd = os.path.join(dirIn, 'single-end')
			dirOutEnd = os.path.join(dirOut, 'single-end')
		Common.dirExist(dirInEnd)
		Common.mkDir(dirOutEnd)
		fileModelList = os.path.join(dirInEnd, 'model.list')
		Common.fileExist(fileModelList)
		
		linesModel = Common.readLinesFromFile(fileModelList)
		for lineModel in linesModel:
			model = lineModel.rstrip()
			dirWorkModel = os.path.join(dirWork, model)
			dirInModel = os.path.join(dirInEnd, model)
			dirOutModel = os.path.join(dirOutEnd, model)
			Common.mkDir(dirWorkModel)
			Common.dirExist(dirInModel)
			Common.mkDir(dirOutModel)

			print('\nProcessing %s ...' % (model))
			dirWorkTrim = os.path.join(dirWorkModel, 'trim')
			fileFastqList = os.path.join(dirInModel, 'fastq.list')
			dirOutTrim = os.path.join(dirOutModel, 'trim')
			Common.mkDir(dirWorkTrim)
			Common.fileExist(fileFastqList)
			Common.mkDir(dirOutTrim)

			linesFastq = Common.readLinesFromFile(fileFastqList)
			n = 0
			for lineFastq in linesFastq:
				n += 1
				fastq = lineFastq.rstrip()
				print('\nProcessing %s ...' % (fastq))
				cmdCurr = '%s %s --output_dir %s' % (cmdComm, fastq, dirOutTrim)
				nameJob = 'trim_galore_%s_%s' % (model, n)
				fileJob = os.path.join(dirWorkTrim, 'job_script_%s.sh' % (nameJob))
				Common.submitJob(fileJob, cmdCurr, nameJob, cores, mem, denovo, submit, dirWorkTrim)

		return 0
	except Exception as e:
		Common.exception(e, 'could not run Trim Galore')


def main():
	parser = argparse.ArgumentParser(description = 'run Trim Galore', \
									formatter_class = argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-w', required = True, \
						help = 'directory containing working folders')
	parser.add_argument('-i', required = True, \
						help = 'directory containing input FastQ files')
	parser.add_argument('-o', required = True, \
						help = 'directory containing output folders')
	parser.add_argument('-p', action = 'store_true', \
						help = 'optional: if specified, run Trim Galore on paired-end FastQ files')
	parser.add_argument('-s', action = 'store_true', help = 'optional: if specified, submit jobs')
	parser.add_argument('-c', type = int, default = 1, help = 'optional: allocate this many cores to a job')
	parser.add_argument('-m', type = int, default = 64, help = 'optional: allocate this much memory to a job')
	parser.add_argument('-d', action = 'store_true', help = 'optional: if specified, submit jobs to the \'denovo\' queue')
	args = parser.parse_args()
	trimGalore(args)
	return 0


if __name__ == '__main__':
	sys.exit(main())


