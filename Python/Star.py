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


def star(args):
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

		cmdComm = '%s --runMode alignReads --runThreadN %s --genomeDir %s --outFilterMultimapNmax 20' \
					' --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999' \
					' --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 1000000' \
					' --alignMatesGapMax 1000000 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33' \
					' --twopassMode Basic --readFilesCommand zcat --alignSoftClipAtReferenceEnds Yes' \
					' --quantMode TranscriptomeSAM GeneCounts --outSAMtype BAM Unsorted SortedByCoordinate' \
					' --outSAMunmapped Within KeepPairs --chimSegmentMin 15 --chimJunctionOverhangMin 15' \
					' --chimOutType WithinBAM SoftClip --outSAMattributes NH HI AS nM NM MD jM jI XS' \
					' --outSAMattrRGline ID:rg1 SM:sm1' % (Common.starCurr, cores, Common.indexSTAR)
		dirInEnd = None
		dirOutEnd = None
		if pairedEnd:
			dirInEnd = os.path.join(dirIn, 'paired-end')
			dirOutEnd = os.path.join(dirOut, 'paired-end')
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
			dirWorkBam = os.path.join(dirWorkModel, 'bam')
			dirInTrim = os.path.join(dirInModel, 'trim')
			dirOutBam = os.path.join(dirOutModel, 'bam')
			Common.mkDir(dirWorkBam)
			Common.dirExist(dirInTrim)
			Common.mkDir(dirOutBam)
			fileFastqList = os.path.join(dirInTrim, 'fastq.list')
			fileFastqPref = os.path.join(dirInTrim, 'fastq.prefix')
			Common.fileExist(fileFastqList)
			Common.fileExist(fileFastqPref)

			linesFastq = Common.readLinesFromFile(fileFastqList)
			linesPrefix = Common.readLinesFromFile(fileFastqPref)
			n = 0
			for (lineFastq, linePrefix) in zip(linesFastq, linesPrefix):
				n += 1
				fastq = lineFastq.rstrip()
				prefix = linePrefix.rstrip()
				print('\nProcessing %s %s ...' % (fastq, prefix))
				cmdCurr = '%s --readFilesIn %s --outFileNamePrefix %s' % \
							(cmdComm, fastq, os.path.join(dirOutBam, prefix))
				nameJob = 'star_%s_%s' % (model, n)
				fileJob = os.path.join(dirWorkBam, 'job_script_%s.sh' % (nameJob))
				Common.submitJob(fileJob, cmdCurr, nameJob, cores, mem, denovo, submit, dirWorkBam)

		return 0
	except Exception as e:
		Common.exception(e, 'could not run STAR')


def main():
	parser = argparse.ArgumentParser(description = 'run STAR', \
									formatter_class = argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-w', required = True, \
						help = 'directory containing working folders')
	parser.add_argument('-i', required = True, \
						help = 'directory containing input FastQ files')
	parser.add_argument('-o', required = True, \
						help = 'directory containing output folders')
	parser.add_argument('-p', action = 'store_true', \
						help = 'optional: if specified, run STAR on paired-end FastQ files')
	parser.add_argument('-s', action = 'store_true', help = 'optional: if specified, submit jobs')
	parser.add_argument('-c', type = int, default = 1, help = 'optional: allocate this many cores to a job')
	parser.add_argument('-m', type = int, default = 64, help = 'optional: allocate this much memory to a job')
	parser.add_argument('-d', action = 'store_true', help = 'optional: if specified, submit jobs to the \'denovo\' queue')
	args = parser.parse_args()
	star(args)
	return 0


if __name__ == '__main__':
	sys.exit(main())


