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


### 92 datasets in total
## 38 datasets
susztak_lab = ['FA1', 'FA2', 'Notch3', 'Notch4', 'PGC1a-1', 'PGC1a-2', 'UUO1', 'UUO2', \
				'longIRI1d1', 'longIRI1d2', 'longIRI3d1', 'longIRI3d2', 'longIRI14d1', 'longIRI14d2', \
				'shortIRI1d1', 'shortIRI1d2', 'shortIRI3d1', 'shortIRI3d2', 'shortIRI14d1', 'shortIRI14d2', \
				'IRI5-4', 'IRI5-5', 'IRI8-3', \
				'ESRRA', 'G2CA-E', 'G2CA-C', 'MEK', 'MCK', 'Ksp', 'Scl', \
				'norm1', 'norm2', 'norm3', 'norm4', 'norm5', 'Pod', 'MMctrl7', 'MMctrl8']
## 24 datasets
gse139107 = ['GSE139107-IRIsham1-1', 'GSE139107-IRIsham1-2', 'GSE139107-IRIsham2', 'GSE139107-IRIsham3', \
				'GSE139107-IRI4h1', 'GSE139107-IRI4h2', 'GSE139107-IRI4h3', \
				'GSE139107-IRI12h1-1', 'GSE139107-IRI12h1-2', 'GSE139107-IRI12h2', 'GSE139107-IRI12h3', \
				'GSE139107-IRI2d1-1', 'GSE139107-IRI2d1-2', 'GSE139107-IRI2d2-1', 'GSE139107-IRI2d2-2', 'GSE139107-IRI2d3', \
				'GSE139107-IRI14d1-1', 'GSE139107-IRI14d1-2', 'GSE139107-IRI14d2', 'GSE139107-IRI14d3', \
				'GSE139107-IRI6w1-1', 'GSE139107-IRI6w1-2', 'GSE139107-IRI6w2', 'GSE139107-IRI6w3']
## 4 datasets
gse140023 = ['GSE140023-sham', 'GSE140023-uuo2', 'GSE140023-uuo7', 'GSE140023-ruuo']
## 17 datasets
gse146912 = ['GSE146912-control1', 'GSE146912-control2', 'GSE146912-control3', \
				'GSE146912-nephritis1d1', 'GSE146912-nephritis1d2', 'GSE146912-nephritis5d1', 'GSE146912-nephritis5d2', \
				'GSE146912-doxorubicin1', 'GSE146912-doxorubicin2', \
				'GSE146912-CD2AP-WT', 'GSE146912-CD2AP-KO', \
				'GSE146912-BTBRobPlus1', 'GSE146912-BTBRobPlus2', \
				'GSE146912-BTBRobob12wk1', 'GSE146912-BTBRobob12wk2', 'GSE146912-BTBRobob21wk1', 'GSE146912-BTBRobob21wk2']
## 7 datasets
gse151658 = ['GSE151658-LPS0hr', 'GSE151658-LPS1hr', 'GSE151658-LPS4hr', 'GSE151658-LPS16hr', 'GSE151658-LPS27hr', \
				'GSE151658-LPS36hr', 'GSE151658-LPS48hr']
## 6 datasets
gse157292 = ['GSE157292-nk1', 'GSE157292-nk2', 'GSE157292-rk1', 'GSE157292-rk2', 'GSE157292-tk1', 'GSE157292-tk2']

### 71 disease datasets
disease_samples = ['FA1', 'FA2', 'Notch3', 'Notch4', 'PGC1a-1', 'PGC1a-2', 'UUO1', 'UUO2', \
					'longIRI1d1', 'longIRI1d2', 'longIRI3d1', 'longIRI3d2', 'longIRI14d1', 'longIRI14d2', \
					'shortIRI1d1', 'shortIRI1d2', 'shortIRI3d1', 'shortIRI3d2', 'shortIRI14d1', 'shortIRI14d2', \
					'IRI5-4', 'IRI5-5', 'IRI8-3', \
					'ESRRA', 'G2CA-E', 'MEK', 'Ksp', \
					'GSE139107-IRI4h1', 'GSE139107-IRI4h2', 'GSE139107-IRI4h3', \
					'GSE139107-IRI12h1-1', 'GSE139107-IRI12h1-2', 'GSE139107-IRI12h2', 'GSE139107-IRI12h3', \
					'GSE139107-IRI2d1-1', 'GSE139107-IRI2d1-2', 'GSE139107-IRI2d2-1', 'GSE139107-IRI2d2-2', \
					'GSE139107-IRI2d3', \
					'GSE139107-IRI14d1-1', 'GSE139107-IRI14d1-2', 'GSE139107-IRI14d2', 'GSE139107-IRI14d3', \
					'GSE139107-IRI6w1-1', 'GSE139107-IRI6w1-2', 'GSE139107-IRI6w2', 'GSE139107-IRI6w3', \
					'GSE140023-uuo2', 'GSE140023-uuo7', 'GSE140023-ruuo', \
					'GSE146912-nephritis1d1', 'GSE146912-nephritis1d2', 'GSE146912-nephritis5d1', 'GSE146912-nephritis5d2', \
					'GSE146912-doxorubicin1', 'GSE146912-doxorubicin2', \
					'GSE146912-CD2AP-KO', \
					'GSE146912-BTBRobob12wk1', 'GSE146912-BTBRobob12wk2', 'GSE146912-BTBRobob21wk1', \
					'GSE146912-BTBRobob21wk2', \
					'GSE151658-LPS1hr', 'GSE151658-LPS4hr', 'GSE151658-LPS16hr', 'GSE151658-LPS27hr', \
					'GSE151658-LPS36hr', 'GSE151658-LPS48hr', \
					'GSE157292-rk1', 'GSE157292-rk2', 'GSE157292-tk1', 'GSE157292-tk2']
### 25 control datasets
control_samples = ['G2CA-C', 'MCK', 'Scl', 'norm1', 'norm2', 'norm3', 'norm4', 'norm5', 'Pod', 'MMctrl7', 'MMctrl8', \
				   'GSE139107-IRIsham1-1', 'GSE139107-IRIsham1-2', 'GSE139107-IRIsham2', 'GSE139107-IRIsham3', \
				   'GSE140023-sham', \
				   'GSE146912-control1', 'GSE146912-control2', 'GSE146912-control3', \
				   'GSE146912-CD2AP-WT', \
				   'GSE146912-BTBRobPlus1', 'GSE146912-BTBRobPlus2', \
				   'GSE151658-LPS0hr', \
				   'GSE157292-nk1', 'GSE157292-nk2']
## 72 single-cell samples
cell_samples = ['FA1', 'FA2', 'Notch3', 'Notch4', 'PGC1a-1', 'PGC1a-2', 'UUO1', 'UUO2', \
				'longIRI1d1', 'longIRI1d2', 'longIRI3d1', 'longIRI3d2', 'longIRI14d1', 'longIRI14d2', \
				'shortIRI1d1', 'shortIRI1d2', 'shortIRI3d1', 'shortIRI3d2', 'shortIRI14d1', 'shortIRI14d2', \
				'IRI5-4', 'IRI5-5', 'IRI8-3', \
				'ESRRA', 'G2CA-E', 'G2CA-C', 'MEK', 'MCK', 'Ksp', 'Scl', \
				'norm1', 'norm2', 'norm3', 'norm4', 'norm5', 'Pod', 'MMctrl7', 'MMctrl8', \
				'GSE140023-sham', 'GSE140023-uuo2', 'GSE140023-uuo7', 'GSE140023-ruuo', \
				'GSE146912-control1', 'GSE146912-control2', 'GSE146912-control3', \
				'GSE146912-nephritis1d1', 'GSE146912-nephritis1d2', \
				'GSE146912-nephritis5d1', 'GSE146912-nephritis5d2', \
				'GSE146912-doxorubicin1', 'GSE146912-doxorubicin2', \
				'GSE146912-CD2AP-WT', 'GSE146912-CD2AP-KO', \
				'GSE146912-BTBRobPlus1', 'GSE146912-BTBRobPlus2', \
				'GSE146912-BTBRobob12wk1', 'GSE146912-BTBRobob12wk2', \
				'GSE146912-BTBRobob21wk1', 'GSE146912-BTBRobob21wk2', \
				'GSE151658-LPS0hr', \
				'GSE151658-LPS1hr', 'GSE151658-LPS4hr', 'GSE151658-LPS16hr', \
				'GSE151658-LPS27hr', 'GSE151658-LPS36hr', 'GSE151658-LPS48hr', \
				'GSE157292-nk1', 'GSE157292-nk2', 'GSE157292-rk1', 'GSE157292-rk2', 'GSE157292-tk1', 'GSE157292-tk2']
## 24 single-nucleus samples
nucleus_samples = ['GSE139107-IRIsham1-1', 'GSE139107-IRIsham1-2', 'GSE139107-IRIsham2', 'GSE139107-IRIsham3', \
				   'GSE139107-IRI4h1', 'GSE139107-IRI4h2', 'GSE139107-IRI4h3', \
				   'GSE139107-IRI12h1-1', 'GSE139107-IRI12h1-2', 'GSE139107-IRI12h2', 'GSE139107-IRI12h3', \
				   'GSE139107-IRI2d1-1', 'GSE139107-IRI2d1-2', 'GSE139107-IRI2d2-1', 'GSE139107-IRI2d2-2', 'GSE139107-IRI2d3', \
				   'GSE139107-IRI14d1-1', 'GSE139107-IRI14d1-2', 'GSE139107-IRI14d2', 'GSE139107-IRI14d3', \
				   'GSE139107-IRI6w1-1', 'GSE139107-IRI6w1-2', 'GSE139107-IRI6w2', 'GSE139107-IRI6w3']

markersKnownDefault = {'GEC/Endo': ['Nrp1', 'Kdr', 'Ehd3', 'Igfbp3'], \
					   'Podo': ['Nphs1', 'Nphs2'], \
					   'PT(S1/S2/S3)': ['Slc27a2', 'Lrp2', 'Slc5a2', 'Snhg11', 'Slc5a12', 'Slc7a13', 'Atp11a', 'Slc22a30'], \
					   'DLOH': ['Slc4a11', 'Ptgds', 'Cp'], \
					   'ALOH': ['Slc12a1', 'Ppp1r1b'], \
					   'DCT/CNT': ['Slc12a3', 'Pvalb', 'Trpv5'], \
					   'CD PC': ['Aqp2', 'Hsd11b2'], \
					   'CD Trans': ['Insrr', 'Rhbg'], \
					   'CD IC(A-IC/B-IC)': ['Atp6v1g3', 'Atp6v0d2', 'Slc4a1', 'Aqp6', 'Slc26a4', 'Hmx2'], \
					   'Fib': ['Col1a1', 'Col3a1', 'Vim', 'Fn1'], \
					   'Mono(Ly6Chi/Ly6Clo)': ['Plac8', 'S100a4', 'F13a1', 'Chil3', 'Ace', 'Treml4'], \
					   'DC(11b+/11b-)/pDC': ['Clec10a', 'Cd209a', 'Xcr1', 'Siglech', 'Ccr9'], \
					   'Baso': ['Fcer1a', 'Mcpt8', 'Csrp3', 'Cd200r3', 'Cyp11a1', 'Ms4a2'], \
					   'Macro': ['C1qa', 'C1qb', 'C1qc'], \
					   'Granul': ['S100a8', 'S100a9', 'Hp'], \
					   'B lymph(B1/B2)': ['Cd79a', 'Cd79b', 'Igkc', 'Ebf1', 'Ighm', 'Igha', 'Jchain', 'Derl3'], \
					   'T lymph': ['Ltb', 'Cxcr6', 'Lef1', 'Cd28', 'Icos', 'Rora', 'Actn2', 'Ly6g5b', 'Cd3g', 'Cd3d', 'Cd8a'], \
					   'NK': ['Ncr1', 'Ccl5', 'Nkg7', 'Gzma'], \
					   'Novel1/Novel2': ['Mki67', 'Cdca3', 'Stmn1', 'Lockd', 'Top2a'], \
					   'Podocyte': ['Pdpn', 'Wt1', 'Mafb', 'Synpo', 'Cdkn1c', 'Ptpro'], \
					   'PEC': ['Cldn1', 'Pax8'], \
					   'Mesangial': ['Pdgfrb', 'Gata3', 'Des', 'Itga8', 'Plvap', 'Prkca', 'Art3', 'Nt5e'], \
					   'Endothelial': ['Flt1', 'Tie1', 'Pecam1', 'Emcn'], \
					   'SMC': ['Acta2', 'Myh11', 'Tagln', 'Ren1', 'Akr1b7', 'Rgs5', 'Rergl', 'Map3k7cl'], \
					   'Immune': ['Ptprc', 'Lyz1', 'Csf1r', 'Itgam', 'Ms4a1'], \
					   'TEC': ['Fxyd2', 'Slc14a2', 'Aqp1', 'Umod']}

markersKnownPT = {'Pan-PT': ['Slc34a1', 'Lrp2', 'Acsm1', 'Acsm2', 'Cpt1a', 'Acox3', 'Slc26a6', 'Slc9a3', 'Glud1', 'Pck1', 'Aqp8', 'Hnf4a', 'Ppara'], \
				  'PCT': ['Slc5a2', 'Slc5a12', 'Adra1a', 'Slc6a19', 'Slc7a8', 'Slc7a9'], \
				  'PST': ['Atp11a', 'Slc13a3', 'Slc16a9', 'Slc27a2', 'Slc7a13', 'Slc22a6', 'Slc1a1'], \
				  'PT progenitors': ['Notch2', 'Lgr4'], \
				  'Injured PT': ['Havcr1', 'Krt20', 'Hspa1a', 'Vcam1', 'Dcdc2a', 'Sema5a']}


def annotation(dirOut, resol, clu2ann):
	try:
		dirRDS = os.path.join(dirOut, 'RDS')
		dirFig = os.path.join(dirOut, 'figures')
		dirClt = os.path.join(dirFig, 'clusters')
		Common.dirExist(dirOut)
		Common.dirExist(dirRDS)
		Common.dirExist(dirFig)
		Common.dirExist(dirClt)

		print('\nUsing resolution %s...' % (resol))

		fileInput = os.path.join(dirRDS, 'after_clustering_res%s.h5ad' % (resol))
		Common.fileExist(fileInput)
		print('\nReading \'%s\'...' % (fileInput))
		adata = scanpy.read(fileInput, cache=True)
		print(adata)
		print(adata.obs)

		## Annotate the cells.
		adata.obs['celltypes'] = adata.obs['clusters'].map(clu2ann).astype('category')
		print('\nAfter annotation:')
		print(adata)
		print(adata.obs)

		## Save the result
		fileOutput = os.path.join(dirRDS, 'after_annotation_res%s.h5ad' % (resol))
		adata.write(fileOutput)
		Common.fileExist(fileOutput)

		## UMAP
		um = scanpy.pl.umap(adata, color='celltypes', use_raw=True, legend_loc='on data', legend_fontsize=6, \
							title='cell types', return_fig=True)
		fileUMAP = os.path.join(dirClt, 'umap_celltypes_res%s.pdf' % (resol))
		um.savefig(fileUMAP, bbox_inches='tight')
		Common.fileExist(fileUMAP)

		return 0
	except Exception as e:
		Common.exception(e, 'could not do the annotation')


def bbknn(dirOut, npcs=50):
	try:
		dirRDS = os.path.join(dirOut, 'RDS')
		dirFig = os.path.join(dirOut, 'figures')
		Common.dirExist(dirOut)
		Common.dirExist(dirRDS)
		Common.mkDir(dirFig)

		print('\nUsing %s PCs...' % (npcs))

		fileInput = os.path.join(dirRDS, 'after_pca_%spcs.h5ad' % (npcs))
		Common.fileExist(fileInput)
		print('\nReading \'%s\'...' % (fileInput))
		adata = scanpy.read(fileInput, cache=True)
		print(adata)
		print(adata.obs)

		## Integrate data using BBKNN
		scanpy.external.pp.bbknn(adata, batch_key='sample', n_pcs=npcs)
		## Embed the graph in 2 dimensions using UMAP
		scanpy.tl.umap(adata)
		## Save the result
		fileOutput = os.path.join(dirRDS, 'after_bbknn_%spcs.h5ad' % (npcs))
		adata.write(fileOutput)
		Common.fileExist(fileOutput)

		## Visualize distributions across batches
		for sampleCurr in list(adata.obs['sample'].values.categories):
			um = scanpy.pl.umap(adata, color='sample', use_raw=True, groups=[sampleCurr], return_fig=True)
			fileUMAP = os.path.join(dirFig, 'umap_batch_%s.png' % (sampleCurr))
			um.savefig(fileUMAP, bbox_inches='tight')
			Common.fileExist(fileUMAP)

		return 0
	except Exception as e:
		Common.exception(e, 'could not run BBKNN')


def cleanLoom(dirOut):
	try:
		dirRDS = os.path.join(dirOut, 'RDS')
		Common.dirExist(dirOut)
		Common.dirExist(dirRDS)

		fileLoomOrig = os.path.join(dirRDS, 'after_QC.loom')
		Common.fileExist(fileLoomOrig)
		print('\nReading \'%s\'...' % (fileLoomOrig))
		adata = scanpy.read_loom(fileLoomOrig)
		print('\nOriginal \'%s\':' % (fileLoomOrig))
		print(adata)
		print(adata.obs)

		print('\nCleaning \'%s\'...' % (fileLoomOrig))
		adata.obs.drop(columns=['ClusterID', 'ClusterName'], inplace=True)
		adata.var.drop(columns=adata.var_keys(), inplace=True)
		adata.layers.pop('norm_data', None)
		print('\nClean \'%s\':' % (fileLoomOrig))
		print(adata)
		print(adata.obs)

		## The scanpy v1.6.0 does not work well with its dependeny h5py v3.0.0, 
		## where strings and indices are stored as byte string rather than character string.
		## Downgrading h5py to v2.10.0 can circumvent this bug.
		fileH5adClean = os.path.join(dirRDS, 'after_QC_clean.h5ad')
		adata.write(fileH5adClean)
		Common.fileExist(fileH5adClean)

		## Because of the bug mentioned above, also save a loom file as the backup.
		fileLoomClean = os.path.join(dirRDS, 'after_QC_clean.loom')
		adata.write_loom(fileLoomClean, write_obsm_varm=True)
		Common.fileExist(fileLoomClean)

		## Double-check the saved files.
		## h5ad file
		del adata
		print('\nReading \'%s\'...' % (fileH5adClean))
		adata = scanpy.read(fileH5adClean, cache=True)
		print(adata)
		print(adata.obs)
		## loom file
		del adata
		print('\nReading \'%s\'...' % (fileLoomClean))
		adata = scanpy.read_loom(fileLoomClean)
		print(adata)
		print(adata.obs)
		return 0
	except Exception as e:
		Common.exception(e, 'could not clean the loom file')


def clustering(dirOut, resol, markers=0):
	try:
		dirRDS = os.path.join(dirOut, 'RDS')
		dirFig = os.path.join(dirOut, 'figures')
		dirClt = os.path.join(dirFig, 'clusters')
		Common.dirExist(dirOut)
		Common.dirExist(dirRDS)
		Common.dirExist(dirFig)
		Common.mkDir(dirClt)

		print('\nUsing resolution %s...' % (resol))
		print('Using markers of type %s...' % (markers))

		fileInput = os.path.join(dirRDS, 'after_neighbor.h5ad')
		Common.fileExist(fileInput)
		print('\nReading \'%s\'...' % (fileInput))
		adata = scanpy.read(fileInput, cache=True)
		print(adata)
		print(adata.obs)

		## Compute clusters using the leiden method and store the results with the name `clusters`.
		scanpy.tl.leiden(adata, resolution=resol, key_added='clusters')
		print('\nAfter leiden clustering:')
		print(adata)
		print(adata.obs)
		## Save the result.
		fileOutput = os.path.join(dirRDS, 'after_clustering_res%s.h5ad' % (resol))
		adata.write(fileOutput)
		Common.fileExist(fileOutput)

		## UMAP
		titlePlot = ('resolution %s, %s clusters' % (resol, adata.obs.clusters.nunique()))
		um = scanpy.pl.umap(adata, color='clusters', use_raw=True, legend_loc='on data', legend_fontsize=6, \
							title=titlePlot, return_fig=True)
		fileUMAP = os.path.join(dirClt, 'umap_clusters_res%s.png' % (resol))
		um.savefig(fileUMAP, bbox_inches='tight')
		Common.fileExist(fileUMAP)

		## UMAP for each cluster
		for clusterCurr in list(adata.obs.clusters.values.categories):
			um = scanpy.pl.umap(adata, color='clusters', use_raw=True, groups=[clusterCurr], na_in_legend=False, \
								legend_loc='on data', legend_fontsize=6, title=titlePlot, return_fig=True)
			fileUMAP = os.path.join(dirClt, 'umap_clusters_res%s_clt%s.png' % (resol, clusterCurr))
			um.savefig(fileUMAP, bbox_inches='tight')
			Common.fileExist(fileUMAP)

		if (0 == markers):
			markers = markersKnownDefault
		elif (1 == markers):
			markers = markersKnownPT
		else:
			Common.error('could not identify marker type \'%s\'' % (markers))
		## dotplot
		dp = scanpy.pl.dotplot(adata, var_names=markers, groupby='clusters', use_raw=True, dendrogram=False, \
								title=titlePlot, return_fig=True)
		fileDotPlot = os.path.join(dirClt, 'dotplot_clusters_res%s.pdf' % (resol))
		dp.style(cmap='viridis_r', dot_edge_color=None, grid=True).savefig(fileDotPlot, bbox_inches='tight')
		Common.fileExist(fileDotPlot)

		## stacked-violin plot
		sv = scanpy.pl.stacked_violin(adata, var_names=markers, groupby='clusters', use_raw=True, dendrogram=False, \
										title=titlePlot, return_fig=True)
		fileStkViol = os.path.join(dirClt, 'stacked_violin_clusters_res%s.pdf' % (resol))
		sv.style(cmap='viridis_r', linewidth=0).savefig(fileStkViol, bbox_inches='tight')
		Common.fileExist(fileStkViol)

		return 0
	except Exception as e:
		Common.exception(e, 'could not do the clustering')


def harmony(dirOut):
	try:
		dirRDS = os.path.join(dirOut, 'RDS')
		Common.dirExist(dirOut)
		Common.dirExist(dirRDS)

		fileInput = os.path.join(dirRDS, 'after_pca.h5ad')
		Common.fileExist(fileInput)
		print('\nReading \'%s\'...' % (fileInput))
		adata = scanpy.read(fileInput, cache=True)
		print(adata)
		print(adata.obs)

		## Integrate data using Harmony
		scanpy.external.pp.harmony_integrate(adata, key='sample', basis='X_pca', adjusted_basis='X_harmony', max_iter_harmony=100)
		## Save the result
		fileOutput = os.path.join(dirRDS, 'after_harmony.h5ad')
		adata.write(fileOutput)
		Common.fileExist(fileOutput)
		print(adata)
		print(adata.obs)

		return 0
	except Exception as e:
		Common.exception(e, 'could not run Harmony')


def integration(dirOut, integ, npcs=50):
	try:
		if ('bbknn' == integ):
			print('\nRunning bbknn(%s, %s)...' % (dirOut, npcs))
			bbknn(dirOut, npcs)

		elif ('harmony' == integ):
			print('\nRunning harmony(%s)...' % (dirOut))
			harmony(dirOut)

		elif ('scanorama' == integ):
			print('\nRunning scanoramaInteg(%s, %s)...' % (dirOut, npcs))
			scanoramaInteg(dirOut, npcs)

		else:
			Common.error('could not identify integration approach \'%s\'' % (integ))

		return 0
	except Exception as e:
		Common.exception(e, 'could not do the integration')


def neighbor(dirOut, integ, npcs=50):
	try:
		dirRDS = os.path.join(dirOut, 'RDS')
		dirFig = os.path.join(dirOut, 'figures')
		Common.dirExist(dirOut)
		Common.dirExist(dirRDS)
		Common.mkDir(dirFig)

		print('\nUsing %s integration...' % (integ))
		print('Using %s PCs...' % (npcs))

		fileInput = os.path.join(dirRDS, 'after_%s.h5ad' % (integ))
		Common.fileExist(fileInput)
		print('\nReading \'%s\'...' % (fileInput))
		adata = scanpy.read(fileInput, cache=True)
		print(adata)
		print(adata.obs)

		## Compute the neighborhood graph
		repUse = 'X_%s' % (integ)
		print(repUse)
		scanpy.pp.neighbors(adata, n_neighbors=30, n_pcs=npcs, use_rep=repUse)
		## Embed the graph in 2 dimensions using UMAP
		scanpy.tl.umap(adata)
		## Save the result
		fileOutput = os.path.join(dirRDS, 'after_neighbor.h5ad')
		adata.write(fileOutput)
		Common.fileExist(fileOutput)
		print(adata)
		print(adata.obs)

		## Visualize distributions across batches
		for batch in list(adata.obs['sample'].values.categories):
			um = scanpy.pl.umap(adata, color='sample', use_raw=True, groups=[batch], na_in_legend=False, return_fig=True)
			fileUMAP = os.path.join(dirFig, 'umap_batch_%s.png' % (batch))
			um.savefig(fileUMAP, bbox_inches='tight')
			Common.fileExist(fileUMAP)

		return 0
	except Exception as e:
		Common.exception(e, 'could not compute the neighborhood graph')


def pca(dirOut, npcs=50):
	try:
		dirRDS = os.path.join(dirOut, 'RDS')
		Common.dirExist(dirOut)
		Common.dirExist(dirRDS)

		print('\nUsing %s PCs...' % (npcs))

		fileInput = os.path.join(dirRDS, 'after_HVG.h5ad')
		Common.fileExist(fileInput)
		print('\nReading \'%s\'...' % (fileInput))
		adata = scanpy.read(fileInput, cache=True)
		print(adata)
		print(adata.obs)

		## Scale each gene to unit variance. 
		## Clip values exceeding standard deviation 10. Setting this can help reduce 
		## the effects of features that are only expressed in a very small number of cells.
		scanpy.pp.scale(adata, max_value=10)
		## Reduce the dimensionality of the data by running principal component analysis (PCA), 
		## which reveals the main axes of variation and denoises the data.
		scanpy.tl.pca(adata, n_comps=npcs)
		## Save the result.
		fileOutput = os.path.join(dirRDS, 'after_pca.h5ad')
		adata.write(fileOutput)
		Common.fileExist(fileOutput)
		print(adata)
		print(adata.obs)

		return 0
	except Exception as e:
		Common.exception(e, 'could not run PCA')


def preprocess(dirOut, sufName = ''):
	try:
		dirRDS = os.path.join(dirOut, 'RDS')
		Common.dirExist(dirOut)
		Common.dirExist(dirRDS)

		fileInput = os.path.join(dirRDS, 'after_QC_clean%s.h5ad' % (sufName))
		Common.fileExist(fileInput)
		print('\nReading \'%s\'...' % (fileInput))
		adata = scanpy.read(fileInput, cache=True)
		print(adata)
		print(adata.obs)

		## Be used later for identifying highly variable genes.
		# adata.layers['counts'] = adata.X.copy()
		## Keep only samples that have >10 cells.
		## Otherwise, scanpy.pp.highly_variable_genes may not work.
		# print(list(adata.obs['sample'].values.categories))
		# print(adata.obs['sample'].value_counts(sort=False))
		# adataSel = adata[adata.obs['sample'].isin(list(adata.obs['sample'].values.categories[adata.obs['sample'].value_counts(sort=False)>10])),:]
		# print(adataSel)
		# print(adataSel.obs)
		# del adata
		# adata = adataSel
		# del adataSel
		# print(adata)
		# print(adata.obs)

		## Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell, 
		## so that counts become comparable among cells. 
		scanpy.pp.normalize_total(adata, target_sum=1e4)
		## Logarithmize the data.
		scanpy.pp.log1p(adata)
		## Identify highly-variable genes.
		## Highly-variable genes are selected within each batch separately and merged.
		## Avoid the selection of batch-specific genes and acts as a lightweight batch correction method.
		scanpy.pp.highly_variable_genes(adata, flavor='seurat', subset=False, batch_key='sample')
		# scanpy.pp.highly_variable_genes(adata, layer='counts', n_top_genes=2000, flavor='seurat_v3', subset=False, batch_key='sample')

		## Set the .raw attribute of AnnData object to the normalized and logarithmized raw gene expression 
		## for later use in differential testing and visualizations of gene expression. 
		## This simply freezes the state of the AnnData object.
		adata.raw = adata
		## Save the result.
		fileOutput = os.path.join(dirRDS, 'after_HVG.h5ad')
		adata.write(fileOutput)
		Common.fileExist(fileOutput)
		print(adata)
		print(adata.obs)

		return 0
	except Exception as e:
		Common.exception(e, 'could not do the preprocessing')


def scanoramaInteg(dirOut, npcs = 50):
	try:
		dirRDS = os.path.join(dirOut, 'RDS')
		dirFig = os.path.join(dirOut, 'figures')
		Common.dirExist(dirOut)
		Common.dirExist(dirRDS)
		Common.mkDir(dirFig)

		print('\nUsing %s PCs...' % (npcs))

		fileInput = os.path.join(dirRDS, 'after_HVG.h5ad')
		Common.fileExist(fileInput)
		print('\nReading \'%s\'...' % (fileInput))
		adata = scanpy.read(fileInput, cache=True)
		print(adata)
		print(adata.obs)

		## Split AnnData object into list of new objects (one object per sample).
		batches = []
		batches.extend(susztak_lab)
		batches.extend(gse139107)
		batches.extend(gse140023)
		batches.extend(gse146912)
		batches.extend(gse151658)
		batches.extend(gse157292)
		adatas = []
		for batch in batches:
			print(batch)
			ad = adata[adata.obs['sample']==batch,:]
			adatas.append(ad)
		for ad in adatas:
			print(list(ad.obs['sample'].values.categories))
		## Run Scanorama integration.
		## Return a list of numpy.ndarrays.
		print('\nRunning integration...')
		datasets_dimred, genes = scanorama.integrate([ad.X for ad in adatas], [ad.var_names.values for ad in adatas], \
														dimred=npcs, union=True)
		print('\nIntegration completed.')
		## Make into one matrix.
		X_dimred = numpy.concatenate(datasets_dimred)
		## Add to the AnnData object
		adata.obsm['X_scanorama'] = X_dimred
		## Compute the neighborhood graph
		scanpy.pp.neighbors(adata, n_neighbors=30, n_pcs=npcs, use_rep='X_scanorama')
		## Embed the graph in 2 dimensions using UMAP
		scanpy.tl.umap(adata)
		## Save the result
		fileOutput = os.path.join(dirRDS, 'after_scanorama_%spcs.h5ad' % (npcs))
		adata.write(fileOutput)
		Common.fileExist(fileOutput)

		## Visualize distributions across batches
		for batch in list(adata.obs['sample'].values.categories):
			um = scanpy.pl.umap(adata, color='sample', use_raw=True, groups=[batch], return_fig=True)
			fileUMAP = os.path.join(dirFig, 'umap_batch_%s.png' % (batch))
			um.savefig(fileUMAP, bbox_inches='tight')
			Common.fileExist(fileUMAP)

		return 0
	except Exception as e:
		Common.exception(e, 'could not run Scanorama')

