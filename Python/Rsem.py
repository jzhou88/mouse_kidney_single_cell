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


def rsem(args):
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

		cmdComm = '%s --num-threads %s --fragment-length-max 1000 --no-bam-output --estimate-rspd' \
					' --forward-prob 0 --append-names --alignments' % (Common.rsemCurr, cores)
		dirInEnd = None
		dirOutEnd = None
		if pairedEnd:
			dirInEnd = os.path.join(dirIn, 'paired-end')
			dirOutEnd = os.path.join(dirOut, 'paired-end')
			cmdComm = '%s --paired-end' % (cmdComm)
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
			dirWorkRsem = os.path.join(dirWorkModel, 'rsem')
			dirInBam = os.path.join(dirInModel, 'bam')
			dirOutRsem = os.path.join(dirOutModel, 'rsem')
			Common.mkDir(dirWorkRsem)
			Common.dirExist(dirInBam)
			Common.mkDir(dirOutRsem)
			fileBamList = os.path.join(dirInBam, 'bam.list')
			fileBamPref = os.path.join(dirInBam, 'bam.prefix')
			Common.fileExist(fileBamList)
			Common.fileExist(fileBamPref)

			linesBam = Common.readLinesFromFile(fileBamList)
			linesPrefix = Common.readLinesFromFile(fileBamPref)
			n = 0
			for (lineBam, linePrefix) in zip(linesBam, linesPrefix):
				n += 1
				bam = lineBam.rstrip()
				prefix = linePrefix.rstrip()
				print('\nProcessing %s %s ...' % (bam, prefix))
				cmdCurr = '%s %s %s %s' % (cmdComm, bam, Common.prefixRsemIndex, os.path.join(dirOutRsem, prefix))
				nameJob = 'rsem_%s_%s' % (model, n)
				fileJob = os.path.join(dirWorkRsem, 'job_script_%s.sh' % (nameJob))
				Common.submitJob(fileJob, cmdCurr, nameJob, cores, mem, denovo, submit, dirWorkRsem)

		return 0
	except Exception as e:
		Common.exception(e, 'could not run RSEM')


def main():
	parser = argparse.ArgumentParser(description = 'run RSEM', \
									formatter_class = argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-w', required = True, \
						help = 'directory containing working folders')
	parser.add_argument('-i', required = True, \
						help = 'directory containing input BAM files')
	parser.add_argument('-o', required = True, \
						help = 'directory containing output folders')
	parser.add_argument('-p', action = 'store_true', \
						help = 'optional: if specified, run RSEM on paired-end BAM files')
	parser.add_argument('-s', action = 'store_true', help = 'optional: if specified, submit jobs')
	parser.add_argument('-c', type = int, default = 1, help = 'optional: allocate this many cores to a job')
	parser.add_argument('-m', type = int, default = 64, help = 'optional: allocate this much memory to a job')
	parser.add_argument('-d', action = 'store_true', help = 'optional: if specified, submit jobs to the \'denovo\' queue')
	args = parser.parse_args()
	rsem(args)
	return 0


if __name__ == '__main__':
	sys.exit(main())


