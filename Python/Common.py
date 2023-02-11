#!/usr/bin/python

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


annotationExpressedRepeats = '/home/jzhou88/installed/ucsc-genome/mm10_rmsk.gtf'
annotationGenome = '/home/jzhou88/installed/refdata-cellranger/refdata-cellranger-mm10-3.0.0/genes/genes.gtf'
indexSTAR = '/home/jzhou88/projects/datasets/From_Hongbo/bulkRNAseq/index/STAR-2.6.1a/mm10'
prefixRsemIndex = '/home/jzhou88/projects/datasets/From_Hongbo/bulkRNAseq/index/RSEM-1.3.1/mm10/gencode_vM18'
refCellRangerAtacMouse = '/home/jzhou88/installed/refdata-cellranger-atac/refdata-cellranger-arc-mm10-2020-A-2.0.0'
refCellRangerMouseCell = '/home/jzhou88/installed/refdata-cellranger/refdata-cellranger-mm10-3.0.0'
refCellRangerMouseNucleus = '/home/jzhou88/installed/refdata-cellranger/mm10-3.0.0_premrna'
rsemCurr = '/home/jzhou88/installed/RSEM-1.3.3/rsem-calculate-expression'
starCurr = '/home/jzhou88/installed/STAR-2.6.1e/bin/Linux_x86_64_static/STAR'
trimGaloreCurr = '/home/jzhou88/installed/TrimGalore-0.6.6/trim_galore'


def dirExist(dirname):
	if not(os.path.isdir(dirname)):
		error('dir \'%s\' does not exist' % (dirname))


def error(msg):
	print('Error: %s.' % (msg))
	sys.exit(-1)


def exception(e, msg):
	print('Exception:', e)
	error(msg)


def fileExist(filename):
	if not(os.path.isfile(filename)):
		error('file \'%s\' does not exist' % (filename))


def mkDir(pathDir):
	try:
		if not(os.path.isdir(pathDir)):
			os.mkdir(pathDir)
		return 0
	except Exception as e:
		exception(e, 'could not make directory \'%s\'' % (pathDir))


def readLinesFromFile(filename, chkLine = True):
	try:
		fileExist(filename)
		f = open(filename, 'r')
		lines = f.readlines()
		f.close()
		if chkLine and (len(lines) <= 0):
			error('file \'%s\' contains nothing' % (filename))
		return lines
	except Exception as e:
		exception(e, 'could not read file \'%s\'' % (filename))


def submitJob(fileJob, cmdRun, nameJob, cores, mem, denovo = False, submit = False, dirJob = None):
	try:
		f = open(fileJob, 'w')
		f.write('#!/bin/bash\n')
		f.write('#BSUB -J %s # LSF job name\n' % (nameJob))
		f.write('#BSUB -o %s.%%J.out # Name of the job output file\n' % (nameJob))
		f.write('#BSUB -e %s.%%J.error # Name of the job error file\n' % (nameJob))
		if (None != dirJob):
			f.write('#BSUB -cwd %s # Current working directory for the job\n' % (dirJob))
		f.write('\n')
		f.write(cmdRun + '\n')
		f.close()
		cores = int(math.ceil(cores))
		mem = int(math.ceil(mem)) * 1024
		cmd = 'bsub -N -u zhoujia@pennmedicine.upenn.edu'
		if (denovo or ((479 * 1024) < mem)):
			cmd += ' -q denovo'
		cmd += ' -M %s' % (mem)
		if (1 < cores):
			cmd += ' -n %s -R "rusage [mem=%s] span[hosts=1]"' % (cores, mem)
		else:
			cmd += ' -R "rusage [mem=%s]"' % (mem)
		cmd += ' < %s' % (fileJob)
		print('\n')
		os.system('cat %s' % (fileJob))
		print(cmd)
		if submit:
			os.system(cmd)
		return 0
	except Exception as e:
		exception(e, 'could not submit job \'%s\'' % (fileJob))

