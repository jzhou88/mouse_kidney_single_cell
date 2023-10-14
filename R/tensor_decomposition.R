# tensor_decomposition.R
base::cat("\nRunning \"tensor_decomposition.R\" ...\n")

sobj <- base::readRDS(file = base::file.path(rds.dir, "after_annotation_res0.8.RDS"))
Seurat::DefaultAssay(object = sobj) <- "RNA"
Seurat::Idents(object = sobj) <- "cell.type"
sobj <- Seurat::RenameIdents(
  object = sobj,
  "PCT" = "PT",
  "PCT/PST" = "PT",
  "PST" = "PT",
  "injured PT" = "PT",
  "CNT/CD PC" = "CD PC",
  "Mono" = "Immune",
  "Macro" = "Immune",
  "Macro/DC" = "Immune",
  "pDC" = "Immune",
  "Baso" = "Immune",
  "Granul" = "Immune",
  "B lymph" = "Immune", 
  "T lymph" = "Immune", 
  "NK" = "Immune"    
)
sobj$ctypes <- Seurat::Idents(object = sobj)
base::table(sobj$cell.type)
base::table(sobj$ctypes)

conditions <- sobj$sample
cell.names <- base::names(x = conditions)
conditions <- base::as.character(x = conditions)
conditions[conditions %in% rna.disease.samples] <- "diseased"
conditions[conditions %in% rna.control.samples] <- "healthy"
base::names(x = conditions) <- cell.names
sobj <- Seurat::AddMetaData(object = sobj, metadata = conditions, col.name = "condition")
base::table(sobj$condition)

curr.rna.control.samples <- base::c("G2CA-C", "GSE151658-LPS0hr", "MMctrl7", "MMctrl8", 
                                    "norm1", "norm2", "norm3", "norm4", "norm5")
models <- sobj$sample
cell.names <- base::names(x = models)
models <- base::as.character(x = models)

models[models %in% curr.rna.control.samples] <- "Control"
models[models %in% base::c("ESRRA")] <- "Esrra"
models[models %in% base::c("FA1", "FA2")] <- "FAN"
models[models %in% base::c("G2CA-E")] <- "APOL1"
models[models %in% base::c("GSE151658-LPS1hr", "GSE151658-LPS4hr", "GSE151658-LPS16hr", 
                           "GSE151658-LPS27hr", "GSE151658-LPS36hr", "GSE151658-LPS48hr")] <- "LPS"
models[models %in% base::c("longIRI1d1", "longIRI1d2", "longIRI3d1", "longIRI3d2", 
                           "longIRI14d1", "longIRI14d2", "shortIRI1d1", "shortIRI1d2", 
                           "shortIRI3d1", "shortIRI3d2", "shortIRI14d1", "shortIRI14d2")] <- "IRI"
models[models %in% base::c("Notch3")] <- "Notch1"
models[models %in% base::c("PGC1a-1", "PGC1a-2")] <- "PGC1a"
models[models %in% base::c("UUO1", "UUO2")] <- "UUO"

base::names(x = models) <- cell.names
sobj <- Seurat::AddMetaData(object = sobj, metadata = models, col.name = "model")
base::table(sobj$model)

strains <- sobj$sample
cell.names <- base::names(x = strains)
strains <- base::as.character(x = strains)
strains[strains %in% select.b6.mouse.scRNA] <- "B6"
strains[strains %in% select.fvb.mouse.scRNA] <- "FVB"
base::names(x = strains) <- cell.names
sobj <- Seurat::AddMetaData(object = sobj, metadata = strains, col.name = "strain")
base::table(sobj$strain)

ages <- sobj$sample
cell.names <- base::names(x = ages)
ages <- base::as.character(x = ages)

ages[ages %in% base::c("ESRRA")] <- "8 weeks"
ages[ages %in% base::c("FA1", "FA2")] <- "5-8 weeks"
ages[ages %in% base::c("G2CA-C", "G2CA-E")] <- "24 weeks"
ages[ages %in% base::c("GSE151658-LPS0hr", "GSE151658-LPS1hr", "GSE151658-LPS4hr", "GSE151658-LPS16hr", 
                       "GSE151658-LPS27hr", "GSE151658-LPS36hr", "GSE151658-LPS48hr")] <- "8-10 weeks"
ages[ages %in% base::c("longIRI1d1", "longIRI1d2", "longIRI3d1", "longIRI3d2", 
                       "longIRI14d1", "longIRI14d2", "shortIRI1d1", "shortIRI1d2", 
                       "shortIRI3d1", "shortIRI3d2", "shortIRI14d1", "shortIRI14d2")] <- "10-12 weeks"
ages[ages %in% base::c("MMctrl7", "MMctrl8")] <- "8-16 weeks"
ages[ages %in% base::c("norm1", "norm2", "norm3", "norm4")] <- "4-8 weeks"
ages[ages %in% base::c("norm5")] <- "10-12 weeks"
ages[ages %in% base::c("Notch3")] <- "6 weeks"
ages[ages %in% base::c("PGC1a-1", "PGC1a-2")] <- "6 weeks"
ages[ages %in% base::c("UUO1", "UUO2")] <- "7 weeks"

base::names(x = ages) <- cell.names
sobj <- Seurat::AddMetaData(object = sobj, metadata = ages, col.name = "age")
base::table(sobj$age)

count.matrix <- Seurat::GetAssayData(object = sobj, slot = "counts", assay = "RNA")
base::dim(x = count.matrix)
metadata <- sobj@meta.data
metadata$donors <- metadata$sample
base::dim(x = metadata)
head(metadata, n=3L)

known.celltypes <- base::c(
  "Endo",
  "Podo",
  "PT",
  "DLOH",
  "ALOH",
  "DCT",
  "CD PC",
  "CD IC",
  "Proliferating",
  "Immune"
)
known.models <- base::c("Control", "FAN", "UUO", "IRI", "Notch1", "PGC1a", "APOL1", "Esrra", "LPS")
known.ages <- base::c("4-8 weeks", "5-8 weeks", "6 weeks", "7 weeks", "8 weeks", 
                      "8-10 weeks", "8-16 weeks", "10-12 weeks", "24 weeks")
curr.rna.control.samples <- base::c("norm1", "norm2", "norm3", "norm4", "norm5", 
                                    "MMctrl7", "MMctrl8", "G2CA-C", "GSE151658-LPS0hr")
curr.rna.disease.samples <- base::c(
  "FA1", "FA2", "UUO1", "UUO2", 
  "longIRI1d1", "longIRI1d2", "longIRI3d1", "longIRI3d2", "longIRI14d1", "longIRI14d2", 
  "shortIRI1d1", "shortIRI1d2", "shortIRI3d1", "shortIRI3d2", "shortIRI14d1", "shortIRI14d2", 
  "Notch3", "PGC1a-1", "PGC1a-2", "G2CA-E", "ESRRA", 
  "GSE151658-LPS1hr", "GSE151658-LPS4hr", "GSE151658-LPS16hr", 
  "GSE151658-LPS27hr", "GSE151658-LPS36hr", "GSE151658-LPS48hr")
known.samples <- base::c(curr.rna.control.samples, curr.rna.disease.samples)

metadata$ctypes <- base::factor(x = metadata$ctypes, levels = known.celltypes)
metadata$condition <- base::factor(x = metadata$condition, levels = base::c("healthy", "diseased"))
metadata$model <- base::factor(x = metadata$model, levels = known.models)
metadata$strain <- base::factor(x = metadata$strain, levels = base::c("B6", "FVB"))
metadata$age <- base::factor(x = metadata$age, levels = known.ages)
metadata$donors <- base::factor(x = metadata$donors, levels = known.samples)

table(metadata$ctypes)
table(metadata$condition)
table(metadata$model)
table(metadata$strain)
table(metadata$age)
table(metadata$donors)

param_list <- scITD::initialize_params(ctypes_use = known.celltypes, ncores = 8, rand_seed = 10)

proj.container <- scITD::make_new_container(params = param_list, count_data = count.matrix, meta_data = metadata, label_donor_sex = F)

proj.container <- scITD::form_tensor(container = proj.container, donor_min_cells = 0, norm_method = 'trim', scale_factor=10000, 
                                     vargenes_method = 'norm_var_pvals', vargenes_thresh = 0.1, batch_var = "sample", 
                                     scale_var = T, var_scale_power = 2, verbose = F)

length(x = proj.container$all_vargenes)

proj.container <- scITD::determine_ranks_tucker(container=proj.container, max_ranks_test=base::c(10,15), shuffle_level='cells', 
                                                num_iter=10, #batch_var = "sample", 
                                                norm_method='trim', scale_factor=10000, scale_var=T, var_scale_power=2)

proj.container$plots$rank_determination_plot

ranks.to.use <- base::c(5, 9)
proj.container <- run_stability_analysis(container=proj.container, ranks=ranks.to.use, subset_type='subset',
                                         sub_prop=0.95, n_iterations=50)

proj.container$plots$stability_plot_dsc

proj.container <- scITD::run_tucker_ica(container=proj.container, ranks=ranks.to.use, tucker_type='regular', rotation_type='hybrid')
proj.container <- scITD::get_meta_associations(container=proj.container, 
                                               vars_test=c('condition', 'model', 'strain', 'age'), 
                                               stat_use='pval')

proj.container <- scITD::plot_donor_matrix(container = proj.container, 
                                           meta_vars=c('condition', 'model', 'strain', 'age'), 
                                           cluster_by_meta = 'model', show_donor_ids = T, 
                                           add_meta_associations='pval')
w <- 8
h <- 12
p <- proj.container$plots$donor_matrix
grDevices::pdf(file = base::file.path(fig.pub.dir, "heatmap_scITD_condition_model_donor_score.pdf"), 
               width = w, height = h, onefile = T, paper = "special")
base::print(x = p)
grDevices::dev.off()

proj.container <- scITD::get_lm_pvals(container = proj.container)
p.thresh <- 0.05
celltypes.to.callout <- base::list(
  "Factor 1" = base::c("PT"),
  "Factor 2" = base::c("PT"),
  "Factor 3" = base::c("PT"),
  "Factor 4" = base::c("PT"),
  "Factor 5" = base::c("PT")
)
proj.container <- scITD::get_all_lds_factor_plots(container = proj.container, use_sig_only=T, nonsig_to_zero=T, 
                                                  annot="sig_genes", sig_thresh=p.thresh, display_genes=F, 
                                                  gene_callouts=T, callout_n_gene_per_ctype=30, callout_ctypes=celltypes.to.callout,
                                                  show_var_explained=T)
p <- render_multi_plots(container = proj.container, data_type='loadings', max_cols = 5)
w <- 42
h <- 7
grDevices::pdf(
  file = base::file.path(fig.pub.dir, 
                         base::sprintf(fmt = "heatmap_scITD_gene_loadings_AllFactors_p%s.pdf", p.thresh)), 
  width = w, height = h, onefile = T, paper = "special"
)
base::print(x = p)
grDevices::dev.off()
