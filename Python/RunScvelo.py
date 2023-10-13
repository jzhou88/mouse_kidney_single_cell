#!/usr/bin/python

import numpy
import os
import pandas
import scipy
import scvelo

scvelo.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scvelo.settings.presenter_view = True  # set max width size for presenter view
scvelo.set_figure_params(style='scvelo')  # for beautified visualization

fileH5ad = '/home/jzhou88/projects/KidneyDisease/2021-11-02/MouseRnaCellPtSeuratHarmonypy.2.2/RDS/after_DDRTree.h5ad'
dirFig = '/home/jzhou88/projects/KidneyDisease/2021-11-02/MouseRnaCellPtSeuratHarmonypy.2.2/figures/trajectory/scvelo'

adata = scvelo.read(fileH5ad, cache=True)
scvelo.utils.show_proportions(adata)
print(adata)
adata.var_names

coresReq = 16
scvelo.pp.filter_and_normalize(adata)
scvelo.pp.moments(adata, n_neighbors=30, n_pcs=None, mode='connectivities', method='umap', use_rep='X_harmony')
scvelo.tl.recover_dynamics(adata, n_jobs=coresReq)
scvelo.tl.velocity(adata, mode='dynamical')
scvelo.tl.velocity_graph(adata, mode_neighbors='connectivities', n_jobs=coresReq)
scvelo.tl.latent_time(adata)

list(adata.obs['cell.type'].values.categories)

celltypes = ['S1', 'S2', 'S2/S3', 'S3', 'Injured1', 'Injured2']
colors7 = ['#619CFF', '#F564E3', '#F8766D', '#B79F00', '#00BA38', '#00BFC4']
for ct in list(adata.obs['cell.type'].values.categories):
	pngFig = os.path.join(dirFig, 'ddrtree_cell.type_%s.png' % (ct))
	if 'S2/S3' == ct:
		pngFig = os.path.join(dirFig, 'ddrtree_cell.type_S2S3.png')
	scvelo.pl.scatter(adata, basis='ddrtree', color='cell.type', palette=colors7, size=None, alpha=None, groups=[ct], sort_order=True, \
					  legend_loc='upper center', legend_fontsize=20, title='', fontsize=None, figsize=None, dpi=400, \
					  frameon=None, show=False, save=pngFig)

for ct in list(adata.obs['cell.type'].values.categories):
	pngFig = os.path.join(dirFig, 'umap_cell.type_%s.png' % (ct))
	if 'S2/S3' == ct:
		pngFig = os.path.join(dirFig, 'umap_cell.type_S2S3.png')
	scvelo.pl.scatter(adata, basis='umap', color='cell.type', palette=colors7, size=None, alpha=None, groups=[ct], sort_order=True, \
					  legend_loc='lower center', legend_fontsize=20, title='', fontsize=None, figsize=None, dpi=400, \
					  frameon=None, show=False, save=pngFig)

list(adata.obs['model.1'].values.categories)

models = ['Control', 'FAN', 'UUO', \
		  'longIRI1d', 'longIRI3d', 'longIRI14d', 'shortIRI1d', 'shortIRI3d', 'shortIRI14d', \
		  'Notch1', 'PGC1a', 'APOL1', 'Esrra', \
		  'LPS1hr', 'LPS4hr', 'LPS16hr', 'LPS27hr', 'LPS36hr', 'LPS48hr']
colors19 = ['#00B5EE', '#F8766D', '#00A7FF', '#E9842C', '#E26EF7', '#7F96FF', '#F863DF', '#FF62BF', '#FF6A9A', '#BC81FF', \
			'#00C0B4', '#00BDD4', '#D69100', '#6FB000', '#BC9D00', '#9CA700', '#00C08E', '#00B813', '#00BD61']
pngFig = os.path.join(dirFig, 'umap_velocity_model.1.png')
scvelo.pl.velocity_embedding_stream(adata, basis='umap', color='model.1', palette=colors19, size=None, alpha=0.3, groups = models, \
									legend_loc='right margin', title='', fontsize=None, figsize=None, dpi=400, frameon=None, show=False, save=pngFig)

markers = ['Havcr1', 'Slc5a2', 'Slc22a30']
for ct in markers:
	pngFig = os.path.join(dirFig, 'umap_%s.png' % (ct))
	dm = scanpy.pl.embedding(adata, basis='umap', color=ct, color_map='gnuplot_r', vmax=4, return_fig=True)
	dm.savefig(pngFig, bbox_inches='tight')

for ct in list(adata.obs['model.1'].values.categories):
	pngFig = os.path.join(dirFig, 'umap_model.1_%s.png' % (ct))
	scvelo.pl.scatter(adata, basis='umap', color='model.1', palette=colors19, size=None, alpha=None, groups=[ct], sort_order=True, \
					  legend_loc='lower center', legend_fontsize=20, title='', fontsize=None, figsize=None, dpi=400, \
					  frameon=None, show=False, save=pngFig)

for ct in list(adata.obs['model.1'].values.categories):
	pngFig = os.path.join(dirFig, 'ddrtree_model.1_%s.png' % (ct))
	scvelo.pl.scatter(adata, basis='ddrtree', color='model.1', palette=colors19, size=None, alpha=None, groups=[ct], sort_order=True, \
					  legend_loc='upper center', legend_fontsize=20, title='', fontsize=None, figsize=None, dpi=400, \
					  frameon=None, show=False, save=pngFig)

