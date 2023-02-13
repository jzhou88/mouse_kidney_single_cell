#!/usr/bin/python

import numpy
import os
import pandas
import scipy
import scvelo

scvelo.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scvelo.settings.presenter_view = True  # set max width size for presenter view
scvelo.set_figure_params(style='scvelo')  # for beautified visualization

fileH5ad = '/home/jzhou88/projects/KidneyDisease/2021-11-02/RnaCellPtSeuratHarmonypy.1/RDS/after_DDRTree.h5ad'
dirFig = '/home/jzhou88/projects/KidneyDisease/2021-11-02/RnaCellPtSeuratHarmonypy.1/figures/trajectory/scvelo'

adata = scvelo.read(fileH5ad, cache=True)
scvelo.utils.show_proportions(adata)
print(adata)

coresReq = 16
scvelo.pp.filter_and_normalize(adata, min_shared_counts=30, n_top_genes=2000)
scvelo.pp.moments(adata, n_neighbors=30, n_pcs=None, mode='connectivities', method='umap', use_rep='X_harmony')
scvelo.tl.recover_dynamics(adata, n_jobs=coresReq)
scvelo.tl.velocity(adata, mode='dynamical')
scvelo.tl.velocity_graph(adata, mode_neighbors='connectivities', n_jobs=coresReq)
scvelo.tl.latent_time(adata)

### '#F8766D','#D89000','#A3A500','#39B600','#00BF7D','#00BFC4','#00B0F6','#9590FF','#E76BF3','#FF62BC'
### P8,       P9,       P4,       P1,       P6,       P10,      P2,       P3,       P5,       P7
### 9         10        5         1         7         2         3         4         6         8
celltypes = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10']
colors10 = ['#39B600', '#00BFC4', '#00B0F6', '#9590FF', '#A3A500', '#E76BF3', '#00BF7D', '#FF62BC', '#F8766D', '#D89000']
pngFig = os.path.join(dirFig, 'ddrtree_velocity_celltype.png')
scvelo.pl.velocity_embedding_stream(adata, basis='ddrtree', color='cell.type', palette=colors10, size=None, alpha=0.3, groups = celltypes, \
                                    legend_loc='right margin', legend_fontsize=20, title='', fontsize=None, figsize=None, dpi=400, frameon=None, \
                                    show=False, save=pngFig)

for ct in list(adata.obs['cell.type'].values.categories):
    pngFig = os.path.join(dirFig, 'ddrtree_celltype_%s.png' % (ct))
    scvelo.pl.scatter(adata, basis='ddrtree', color='cell.type', palette=colors10, size=None, alpha=None, groups=[ct], sort_order=True, \
                      legend_loc='none', title='', fontsize=None, figsize=None, dpi=400, frameon=None, show=False, save=pngFig)

### "#F8766D", "#C49A00", "#53B400", "#00C094",    "#00B6EB",     "#A58AFF", "#FB61D7"
### 'Control', 'FAN',     'UUO',     'longIRI14d', 'shortIRI14d', 'Notch',   'PGC1a'
### 1,         2,         5,         6,            7,             3,         4
models = ['Control', 'FAN', 'UUO', 'longIRI14d', 'shortIRI14d', 'Notch', 'PGC1a']
colors7 = ["#F8766D", "#C49A00", "#A58AFF", "#FB61D7", "#53B400", "#00C094", "#00B6EB"]
pngFig = os.path.join(dirFig, 'ddrtree_velocity_model.png')
scvelo.pl.velocity_embedding_stream(adata, basis='ddrtree', color='model', palette=colors7, size=None, alpha=0.3, groups = models, \
                                    legend_loc='right margin', legend_fontsize=20, title='', fontsize=None, figsize=None, dpi=400, frameon=None, \
                                    show=False, save=pngFig)

for m in list(adata.obs['model'].values.categories):
    pngFig = os.path.join(dirFig, 'ddrtree_model_%s.png' % (m))
    scvelo.pl.scatter(adata, basis='ddrtree', color='model', palette=colors7, size=None, alpha=None, groups=[m], sort_order=True, \
                      legend_loc='none', title='', fontsize=None, figsize=None, dpi=400, frameon=None, show=False, save=pngFig)

colormap = 'gnuplot'
pngFig = os.path.join(dirFig, 'ddrtree_latent_time_unsorted.png')
scvelo.pl.scatter(adata, basis='ddrtree', color='latent_time', color_map=colormap, colorbar=None, size=80, sort_order=False, title='', \
                  fontsize=None, figsize=None, dpi=400, frameon=None, show=False, save=pngFig)
