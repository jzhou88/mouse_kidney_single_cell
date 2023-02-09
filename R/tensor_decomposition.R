# tensor_decomposition.R
base::cat("\nRunning \"tensor_decomposition.R\" ...\n")

sobj <- base::readRDS(file = base::file.path(rds.dir, "after_annotation_res1.7.RDS"))
Seurat::DefaultAssay(object = sobj) <- "RNA"
Seurat::Idents(object = sobj) <- "adj.cell.type"
sobj <- Seurat::RenameIdents(object = sobj, "Novel2" = "Proliferating", "PT S1" = "PT", "PT S2" = "PT", "PT S3" = "PT")
sobj$ctypes <- Seurat::Idents(object = sobj)

select.samples <- base::c(base::unlist(x = select.rna.disease.samples, recursive = T, use.names = F), select.rna.control.samples)
select.sobj <- subset(x = sobj, subset = (sample %in% select.samples))

models <- select.sobj$sample
cell.names <- base::names(x = models)
models <- base::as.character(x = models)
disease.models <- base::names(x = select.rna.disease.samples)
for (dm in disease.models) {
  models[models %in% select.rna.disease.samples[[dm]]] <- dm
}
models[models %in% select.rna.control.samples] <- "Control"
base::names(x = models) <- cell.names
select.sobj <- Seurat::AddMetaData(object = select.sobj, metadata = models, col.name = "model")

conditions <- select.sobj$sample
cell.names <- base::names(x = conditions)
conditions <- base::as.character(x = conditions)
conditions[conditions %in% rna.disease.samples] <- "diseased"
conditions[conditions %in% rna.control.samples] <- "healthy"
base::names(x = conditions) <- cell.names
select.sobj <- Seurat::AddMetaData(object = select.sobj, metadata = conditions, col.name = "condition")

orig.sel.sobj <- select.sobj

count.matrix <- Seurat::GetAssayData(object = select.sobj, slot = "counts", assay = "RNA")
metadata <- select.sobj@meta.data
metadata$donors <- metadata$sample

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
  "Mono", 
  "DC", 
  "Macro", 
  "Granul", 
  "B lymph", 
  "T lymph", 
  "NK"
)
known.models <- base::c("Control", base::names(x = select.rna.disease.samples))
known.samples <- base::c(select.rna.control.samples, base::unlist(x = select.rna.disease.samples, recursive = T, use.names = F))

metadata$ctypes <- base::factor(x = metadata$ctypes, levels = known.celltypes)
metadata$model <- base::factor(x = metadata$model, levels = known.models)
metadata$condition <- base::factor(x = metadata$condition, levels = base::c("healthy", "diseased"))
metadata$donors <- base::factor(x = metadata$donors, levels = known.samples)

param_list <- scITD::initialize_params(ctypes_use = known.celltypes, ncores = 16, rand_seed = 10)

proj.container <- scITD::make_new_container(params = param_list, count_data = count.matrix, meta_data = metadata, label_donor_sex = F)

proj.container <- scITD::form_tensor(container = proj.container, donor_min_cells = 0, norm_method = 'trim', scale_factor=10000, 
                                     vargenes_method = 'norm_var_pvals', vargenes_thresh = 0.1, batch_var = "sample", 
                                     scale_var = T, var_scale_power = 2, verbose = F)

proj.container <- scITD::determine_ranks_tucker(container=proj.container, max_ranks_test=base::c(10,15), shuffle_level='cells', 
                                                num_iter=10, norm_method='trim', scale_factor=10000, scale_var=T, var_scale_power=2)
proj.container$plots$rank_determination_plot

ranks.to.use <- base::c(4, 6)
proj.container <- run_stability_analysis(container=proj.container, ranks=ranks.to.use, subset_type='subset', sub_prop=0.95, n_iterations=50)
proj.container$plots$stability_plot_dsc

proj.container <- scITD::run_tucker_ica(container=proj.container, ranks=ranks.to.use, tucker_type='regular', rotation_type='hybrid')
proj.container <- scITD::get_meta_associations(container=proj.container, vars_test=c('condition', 'model'), stat_use='pval')
orig.proj.container <- proj.container

proj.container <- orig.proj.container
proj.container <- scITD::plot_donor_matrix(container = proj.container, meta_vars=c('condition', 'model'), cluster_by_meta = 'model', 
                                           show_donor_ids = F, add_meta_associations='pval')
proj.container$plots$donor_matrix

proj.container <- scITD::get_lm_pvals(container = proj.container)
orig.proj.container.1 <- proj.container

proj.container <- orig.proj.container.1
p.thresh <- 0.001
proj.container <- scITD::get_all_lds_factor_plots(container = proj.container, use_sig_only=T, nonsig_to_zero=T, 
                                                  annot="sig_genes", sig_thresh=p.thresh, display_genes=F, 
                                                  gene_callouts=T, callout_n_gene_per_ctype=5, 
                                                  show_var_explained=T)
for (i in base::names(x = proj.container$plots$all_lds_plots)) {
  p <- proj.container$plots$all_lds_plots[[i]]
  base::print(x = p)
}
render_multi_plots(container = proj.container, data_type='loadings', max_cols = 4)

orig.proj.container.2 <- proj.container

proj.container <- orig.proj.container.2
curr.factor <- 1
curr.celltype <- "PT"
curr.gene <- "Il1b"

meta <- proj.container$scMinimal_full$metadata[,c('donors', 'model')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL

d_exp <- proj.container[["scMinimal_ctype"]][[curr.celltype]][["pseudobulk"]][,curr.gene]
dsc <- proj.container$tucker_results[[1]][,curr.factor]
tmp <- cbind.data.frame(dsc[names(d_exp)],d_exp,meta[names(d_exp),1])
colnames(tmp) <- c('dscore','expres','model')

# add regression line
lmres <- lm(expres~dscore,data=tmp)
line_range <- seq(min(tmp$dscore),max(tmp$dscore),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')

mycol <- scales::hue_pal()(length(known.models))
names(mycol) <- known.models
f <- 18
ps <- 3.5
ggplot(tmp,aes(x=dscore,y=expres)) +
geom_point(mapping=aes(color=model),alpha = 0.75,pch=19,size=ps) +
ggpubr::stat_cor(method="pearson",label.x.npc="left",label.y.npc="top",output.type='expression',p.accuracy=0.001) +
geom_line(data=line_df,aes(x=myx,y=myy,color="#000000")) +
scale_color_manual(values=mycol) + 
ylab(base::sprintf(fmt = '%s expression (%s)', curr.gene, curr.celltype)) +
xlab(base::sprintf(fmt = 'Factor %s sample scores', curr.factor)) +
theme_bw() + 
ggplot2::theme(
  text = ggplot2::element_text(size = f),
  axis.title = ggplot2::element_text(face = "bold"),
  legend.title = ggplot2::element_text(face = "bold"))

proj.container <- orig.proj.container.2
curr.factor <- 2
curr.celltype <- "PT"
curr.gene <- "Isg15"

meta <- proj.container$scMinimal_full$metadata[,c('donors', 'model')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL

d_exp <- proj.container[["scMinimal_ctype"]][[curr.celltype]][["pseudobulk"]][,curr.gene]
dsc <- proj.container$tucker_results[[1]][,curr.factor]
tmp <- cbind.data.frame(dsc[names(d_exp)],d_exp,meta[names(d_exp),1])
colnames(tmp) <- c('dscore','expres','model')

# add regression line
lmres <- lm(expres~dscore,data=tmp)
line_range <- seq(min(tmp$dscore),max(tmp$dscore),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')

mycol <- scales::hue_pal()(length(known.models))
names(mycol) <- known.models
f <- 18
ps <- 3.5
ggplot(tmp,aes(x=dscore,y=expres)) +
geom_point(mapping=aes(color=model),alpha = 0.75,pch=19,size=ps) +
ggpubr::stat_cor(method="pearson",label.x.npc="left",label.y.npc="top",output.type='expression',p.accuracy=0.001) +
geom_line(data=line_df,aes(x=myx,y=myy,color="#000000")) +
scale_color_manual(values=mycol) + 
ylab(base::sprintf(fmt = '%s expression (%s)', curr.gene, curr.celltype)) +
xlab(base::sprintf(fmt = 'Factor %s sample scores', curr.factor)) +
theme_bw() + 
ggplot2::theme(
  text = ggplot2::element_text(size = f),
  axis.title = ggplot2::element_text(face = "bold"),
  legend.title = ggplot2::element_text(face = "bold"))
