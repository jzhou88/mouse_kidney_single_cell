# monocle2_analysis.R
base::cat("\nRunning \"monocle2_analysis.R\" ...\n")

### Define variables ------
sobj <- NULL
cds <- NULL
downsample.class <- NULL
orig.sobj <- NULL
clustering <- NULL
annotating <- NULL
root.state <- NULL
branch.points <- NULL
lineages <- NULL

act.assay <- "RNA"
fn.suf <- ""
max.cell.n <- 10000 # maximum number of cells used in trajectory analysis
set.w <- 0.003
set.h <- 0.1

lower.mito <- F
model.only <- F
calc.deg <- F
mito.ribo.only <- F
excl.mito.ribo <- F
stop.after.ordering.genes <- F

do.traj <- T
do.reorder <- F
do.branch <- F
do.lineage <- F

plot.cluster <- F
plot.celltype <- T
plot.sample <- F
plot.species <- F
plot.state <- T
plot.pseudotime <- T
plot.condition <- T
plot.model <- T
plot.disease <- F
plot.cycle <- F
plot.branch <- F
plot.lineage <- F

export.ddrtree <- F
dump.h5ad <- F

if (do.traj) {
  sobj <- base::readRDS(file = base::file.path(rds.dir, "after_monocle2_downsampling.RDS"))
  
} else {
  cds <- base::readRDS(file = base::file.path(rds.dir, "monocle2_DDRTree.RDS"))
}
if (plot.celltype) {
  annotating <- "cell.type"
}

if (do.traj) {
  ### Construct trajectories ------
  
  orig.cell.n <- base::ncol(x = sobj)
  if (lower.mito) {
    ### Keep cells with lower mitochondrial gene ratio
    base::cat("\n")
    base::print(x = sobj)
    base::cat("\n")
    base::print(x = utils::head(x = sobj[[]], n = 3L))
    base::print(x = utils::tail(x = sobj[[]], n = 3L))
    max.mito <- 0.2
    sobj <- subset(x = sobj, subset = (mitoRatio < max.mito))
    base::cat("\n")
    base::print(x = sobj)
    base::cat("\n")
    base::print(x = utils::head(x = sobj[[]], n = 3L))
    base::print(x = utils::tail(x = sobj[[]], n = 3L))
  }
  if (model.only) {
    ### Keep cells of select disease models and control samples
    select.samples <- base::c(base::unlist(x = select.rna.disease.samples, recursive = T, use.names = F), select.rna.control.samples)
    base::cat("\n")
    base::print(x = select.samples)
    Seurat::Idents(object = sobj) <- "sample"
    base::cat("\n")
    base::print(x = base::table(Seurat::Idents(object = sobj)))
    sobj <- subset(x = sobj, subset = (sample %in% select.samples))
    base::cat("\n")
    base::print(x = base::table(Seurat::Idents(object = sobj)))
  }
  if ((!base::is.null(x = downsample.class)) && (base::ncol(x = sobj) > max.cell.n)) {
    ### Downsample cells
    Seurat::Idents(object = sobj) <- downsample.class
    base::cat("\n")
    base::print(x = sobj)
    base::cat("\n")
    base::print(x = utils::head(x = sobj[[]], n = 3L))
    base::print(x = utils::tail(x = sobj[[]], n = 3L))
    base::cat("\n")
    base::print(x = base::table(Seurat::Idents(object = sobj)))
    class.n <- base::length(x = base::levels(x = sobj))
    avg.cell.n <- base::floor(x = max.cell.n / class.n)
    base::cat(base::sprintf(fmt = "\nmax.cell.n = %s\n", max.cell.n))
    base::cat(base::sprintf(fmt = "class.n = %s\n", class.n))
    base::cat(base::sprintf(fmt = "avg.cell.n = %s\n", avg.cell.n))
    sobj <- subset(x = sobj, downsample = avg.cell.n, seed = 1)
    base::cat("\n")
    base::print(x = base::table(Seurat::Idents(object = sobj)))
    #sobj[[downsample.class]] <- NULL
    base::cat("\n")
    base::print(x = sobj)
    base::cat("\n")
    base::print(x = utils::head(x = sobj[[]], n = 3L))
    base::print(x = utils::tail(x = sobj[[]], n = 3L))
  }
  if (base::ncol(x = sobj) < orig.cell.n) {
    ### Save downsampled cells
    base::saveRDS(object = sobj, file = base::file.path(rds.dir, "after_monocle2_downsampling.RDS"))
  }
  ### Store Data in a CellDataSet Object
  base::cat("\n")
  base::print(x = utils::head(x = sobj@meta.data))
  pd <- methods::new(Class = "AnnotatedDataFrame", data = sobj@meta.data)
  fd <- base::data.frame(gene_short_name = base::rownames(x = sobj), row.names = base::rownames(x = sobj))
  base::cat("\n")
  base::print(x = utils::head(x = fd))
  fd <- methods::new(Class = "AnnotatedDataFrame", data = fd)
  count.matrix <- Seurat::GetAssayData(object = sobj, slot = "counts", assay = act.assay)
  sobj <- NULL
  base::gc()
  base::cat("\n")
  base::print(x = count.matrix[1:6,1:6])
  ## Got lowerDetectionLimit = 0.5 from the tutorial
  cds <- monocle::newCellDataSet(cellData = count.matrix, phenoData = pd, featureData = fd, 
                                 lowerDetectionLimit = 0.5, expressionFamily = VGAM::negbinomial.size())
  base::cat("\n")
  base::print(x = cds)
  base::cat("\n")
  base::print(x = utils::head(x = Biobase::pData(object = cds), n = 3L))
  ### Estimate size factors and dispersions
  cds <- BiocGenerics::estimateSizeFactors(object = cds)
  cds <- BiocGenerics::estimateDispersions(object = cds)
  base::cat("\n")
  base::print(x = utils::head(x = Biobase::pData(object = cds), n = 3L))
  ### Order cells in pseudotime along a trajectory
  ordering_genes <- NULL
  if (calc.deg) {
    tmp <- monocle::differentialGeneTest(cds = cds, fullModelFormulaStr = "~cell.type + mitoRatio", 
                                         reducedModelFormulaStr = "~mitoRatio", relative_expr = T, cores = req.cores, verbose = F)
    celltype_DEG_genes <- tmp[base::order(tmp$qval, decreasing = F),]
    base::cat("\n")
    base::print(x = base::nrow(x = celltype_DEG_genes))
    base::print(x = utils::head(x = celltype_DEG_genes))
    ordering_genes <- base::row.names(x = celltype_DEG_genes)[1:1000]
    base::cat("\nOrdering cells based on differentially expressed genes ...\n")
    base::cat("\n")
    base::print(x = base::length(x = ordering_genes))
    base::print(x = utils::head(x = ordering_genes, n = 3L))
    base::print(x = utils::tail(x = ordering_genes, n = 3L))
    
  } else {
    disp_table <- monocle::dispersionTable(cds = cds)
    base::cat("\n")
    base::print(x = base::nrow(x = disp_table))
    base::print(x = utils::head(x = disp_table))
    ## Got these thresholds from Jihwan's Cell Metabolism 2021 paper
    # ordering_genes <- subset(x = disp_table, subset = ((mean_expression >= 0.05) & (dispersion_empirical >= 2)))$gene_id
    ## Got mean_expression >= 0.5 and dispersion_empirical >= (1 * dispersion_fit) from the tutorial
    # ordering_genes <- base::as.character(x = subset(x = disp_table, subset = ((mean_expression >= mean.expr.cutoff) & 
    #                                                                       (dispersion_empirical >= (1 * dispersion_fit))))$gene_id)
    ordering_genes <- base::as.character(x = subset(x = disp_table, subset = ((mean_expression >= 0.01) & (dispersion_empirical >= 2)))$gene_id)
    base::cat("\nOrdering cells based on select genes with high dispersion ...\n")
    base::cat("\n")
    base::print(x = base::length(x = ordering_genes))
    base::print(x = utils::head(x = ordering_genes, n = 3L))
    base::print(x = utils::tail(x = ordering_genes, n = 3L))
  }
  if (mito.ribo.only) {
    all.mito.ribo.genes <- base::readRDS(file = "/home/jzhou88/projects/KidneyDisease/2020-11-04/ScanpyHarmony14/RDS/mitoANDribo_genes.RDS")
    ordering_genes <- subset(x = ordering_genes, subset = (ordering_genes %in% all.mito.ribo.genes))
    base::cat("\nOrdering cells based on mitochondrial and ribosomal protein ordering genes ...\n")
    
  } else if (excl.mito.ribo) {
    all.mito.ribo.genes <- base::readRDS(file = "/home/jzhou88/projects/KidneyDisease/2020-11-04/ScanpyHarmony14/RDS/mitoANDribo_genes.RDS")
    ordering_genes <- subset(x = ordering_genes, subset = !(ordering_genes %in% all.mito.ribo.genes))
    base::cat("\nOrdering cells based on ordering genes that are neither mitochondrial nor ribosomal protein genes ...\n")
    
  } else {
    base::cat("\nOrdering cells based on all identified ordering genes ...\n")
  }
  base::print(x = base::length(x = ordering_genes))
  base::print(x = utils::head(x = ordering_genes, n = 3L))
  base::print(x = utils::tail(x = ordering_genes, n = 3L))
  cds <- monocle::setOrderingFilter(cds = cds, ordering_genes = ordering_genes)
  base::cat("\n")
  base::print(x = base::table(Biobase::fData(object = cds)$use_for_ordering))
  ogp <- monocle::plot_ordering_genes(cds = cds)
  grDevices::pdf(file = base::file.path(fig.mo2.dir, "ordering_genes_mean_vs_dispersion.pdf"), 
                 width = 7, height = 7, onefile = T, paper = "special")
  base::print(x = ogp)
  grDevices::dev.off()
  if (stop.after.ordering.genes) {
    base::stop("stop after ordering genes.")
  }
  ## residualModelFormulaStr = NULL, because batch correction does not give reasonable results
  cds <- monocle::reduceDimension(cds = cds, max_components = 2, reduction_method = "DDRTree", norm_method = "log", 
                                  residualModelFormulaStr = NULL, pseudo_expr = 1, verbose = F)
  cds <- monocle::orderCells(cds = cds, reverse = F)
  base::saveRDS(object = cds, file = base::file.path(rds.dir, "monocle2_DDRTree.RDS"))
  
} else if (do.reorder) {
  ### Reorder cells ------
  
  base::cat("\n")
  base::print(x = root.state)
  cds <- monocle::orderCells(cds = cds, root_state = root.state)
  base::saveRDS(object = cds, file = base::file.path(rds.dir, "monocle2_DDRTree.RDS"))
}

if (do.branch) {
  ### Analyze branches ------
  
  genes.to.test <- base::row.names(x = Biobase::fData(object = cds))
  base::cat("\n")
  base::print(x = base::length(x = genes.to.test))
  base::print(x = utils::head(x = genes.to.test, n = 3L))
  base::print(x = utils::tail(x = genes.to.test, n = 3L))
  if (mito.ribo.only) {
    all.mito.ribo.genes <- base::readRDS(file = "/home/jzhou88/projects/KidneyDisease/2020-11-04/ScanpyHarmony14/RDS/mitoANDribo_genes.RDS")
    genes.to.test <- subset(x = genes.to.test, subset = (genes.to.test %in% all.mito.ribo.genes))
    base::cat("\nRunning BEAM on mitochondrial and ribosomal protein genes ...\n")
    
  } else if (excl.mito.ribo) {
    all.mito.ribo.genes <- base::readRDS(file = "/home/jzhou88/projects/KidneyDisease/2020-11-04/ScanpyHarmony14/RDS/mitoANDribo_genes.RDS")
    genes.to.test <- subset(x = genes.to.test, subset = !(genes.to.test %in% all.mito.ribo.genes))
    base::cat("\nRunning BEAM on genes that are neither mitochondrial nor ribosomal protein genes ...\n")
    
  } else {
    base::cat("\nRunning BEAM on all the genes ...\n")
  }
  base::print(x = base::length(x = genes.to.test))
  base::print(x = utils::head(x = genes.to.test, n = 3L))
  base::print(x = utils::tail(x = genes.to.test, n = 3L))
  sub.cds <- cds[genes.to.test,]
  base::cat("\n")
  base::print(x = sub.cds)
  
  BEAM_res <- base::list()
  for (bp in branch.points) {
    bp.name <- base::sprintf(fmt = "BP %s", bp)
    base::cat(base::sprintf(fmt = "\nProcessing Branch point %s ...\n", bp))
    tmp <- monocle::BEAM(cds = sub.cds, branch_point = bp, relative_expr = T, verbose = F, cores = req.cores)
    BEAM_res[[bp.name]] <- tmp[base::order(tmp$qval, decreasing = F),]
  }
  base::saveRDS(object = BEAM_res, file = base::file.path(rds.dir, base::sprintf(fmt = "branch_dependent_genes%s.RDS", fn.suf)))
  openxlsx::write.xlsx(x = BEAM_res, file = base::file.path(rst.mo2.dir, base::sprintf(fmt = "branch_dependent_genes%s.xlsx", fn.suf)), 
                       startCol = 1, startRow = 1, colNames = T, rowNames = T, overwrite = T)
}

if (do.lineage) {
  ### Analyze lineages ------
  
  genes.to.test <- base::row.names(x = Biobase::fData(object = cds))
  base::cat("\n")
  base::print(x = base::length(x = genes.to.test))
  base::print(x = utils::head(x = genes.to.test, n = 3L))
  base::print(x = utils::tail(x = genes.to.test, n = 3L))
  if (mito.ribo.only) {
    all.mito.ribo.genes <- base::readRDS(file = "/home/jzhou88/projects/KidneyDisease/2020-11-04/ScanpyHarmony14/RDS/mitoANDribo_genes.RDS")
    genes.to.test <- subset(x = genes.to.test, subset = (genes.to.test %in% all.mito.ribo.genes))
    base::cat("\nRunning differential gene test on mitochondrial and ribosomal protein genes ...\n")
    
  } else if (excl.mito.ribo) {
    all.mito.ribo.genes <- base::readRDS(file = "/home/jzhou88/projects/KidneyDisease/2020-11-04/ScanpyHarmony14/RDS/mitoANDribo_genes.RDS")
    genes.to.test <- subset(x = genes.to.test, subset = !(genes.to.test %in% all.mito.ribo.genes))
    base::cat("\nRunning differential gene test on genes that are neither mitochondrial nor ribosomal protein genes ...\n")
    
  } else {
    base::cat("\nRunning differential gene test on all the genes ...\n")
  }
  base::print(x = base::length(x = genes.to.test))
  base::print(x = utils::head(x = genes.to.test, n = 3L))
  base::print(x = utils::tail(x = genes.to.test, n = 3L))
  
  diff_test_res <- base::list()
  for (l.name in base::names(x = lineages)) {
    base::cat(base::sprintf(fmt = "\nProcessing Lineage %s ...\n", l.name))
    lineage.cells <- base::row.names(x = subset(x = Biobase::pData(object = cds), subset = (State %in% lineages[[l.name]])))
    base::cat("\n")
    base::print(x = base::length(x = lineage.cells))
    sub.cds <- cds[genes.to.test, lineage.cells]
    base::cat("\n")
    base::print(x = sub.cds)
    base::cat("\n")
    base::print(x = base::unique(x = Biobase::pData(object = sub.cds)$State))
    
    tmp <- monocle::differentialGeneTest(cds = sub.cds, fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                         relative_expr = T, cores = req.cores, verbose = F)
    diff_test_res[[l.name]] <- tmp[base::order(tmp$qval, decreasing = F),]
  }
  base::saveRDS(object = diff_test_res, file = base::file.path(rds.dir, base::sprintf(fmt = "pseudotime_dependent_genes%s.RDS", fn.suf)))
  openxlsx::write.xlsx(x = diff_test_res, file = base::file.path(rst.mo2.dir, base::sprintf(fmt = "pseudotime_dependent_genes%s.xlsx", fn.suf)), 
                       startCol = 1, startRow = 1, colNames = T, rowNames = T, overwrite = T)
}

if (plot.cluster) {
  ### Plot cluster ------
  ### Visualize the trajectory in the reduced dimensional space, coloring the cells by clusters
  PlotMonocle2(cds = cds, feature.to.plot = clustering, save.fig.dir = fig.mo2.dir, feature.name = "cluster", fn.suf = fn.suf, 
               plot.together = T, plot.split = T)
}
if (plot.celltype) {
  ### Plot cell type ------
  ### Visualize the trajectory in the reduced dimensional space, coloring the cells by cell types
  PlotMonocle2(cds = cds, feature.to.plot = annotating, save.fig.dir = fig.mo2.dir, feature.name = "celltype", fn.suf = fn.suf, 
               plot.together = T, plot.split = T)
}
if (plot.sample) {
  ### Plot sample ------
  ### Visualize the trajectory in the reduced dimensional space, coloring the cells by samples
  PlotMonocle2(cds = cds, feature.to.plot = "sample", save.fig.dir = fig.mo2.dir, feature.name = "sample", fn.suf = fn.suf, 
               plot.together = F, plot.split = T)
}
if (plot.species) {
  ### Plot species ------
  ### Visualize the trajectory in the reduced dimensional space, coloring the cells by species
  PlotMonocle2(cds = cds, feature.to.plot = "species", save.fig.dir = fig.mo2.dir, feature.name = "species", fn.suf = fn.suf, 
               plot.together = T, plot.split = T)
}
if (plot.state) {
  ### Plot state ------
  ### Visualize the trajectory in the reduced dimensional space, coloring the cells by states
  PlotMonocle2(cds = cds, feature.to.plot = "State", save.fig.dir = fig.mo2.dir, feature.name = "state", fn.suf = fn.suf, 
               plot.together = T, plot.split = T)
}
if (plot.pseudotime) {
  ### Plot pseudotime ------
  ### Visualize the trajectory in the reduced dimensional space, coloring the cells by pseudotime
  PlotMonocle2(cds = cds, feature.to.plot = "Pseudotime", save.fig.dir = fig.mo2.dir, feature.name = "pseudotime", fn.suf = fn.suf, 
               plot.together = T, plot.split = F)
}
if (plot.condition) {
  ### Plot condition ------
  ### Visualize the trajectory in the reduced dimensional space, coloring the cells by condition
  conditions <- base::as.character(x = Biobase::pData(object = cds)$sample)
  base::cat("\n")
  base::print(x = base::table(conditions))
  conditions[conditions %in% rna.disease.samples] <- "disease"
  conditions[conditions %in% rna.control.samples] <- "control"
  base::cat("\n")
  base::print(x = base::table(conditions))
  Biobase::pData(object = cds)$condition <- base::as.character(x = conditions)
  base::cat("\n")
  base::print(x = utils::head(x = Biobase::pData(object = cds), n = 3L))
  base::print(x = utils::tail(x = Biobase::pData(object = cds), n = 3L))
  
  PlotMonocle2(cds = cds, feature.to.plot = "condition", save.fig.dir = fig.mo2.dir, feature.name = "condition", fn.suf = fn.suf, 
               plot.together = T, plot.split = T)
}

if (plot.model) {
  ### Plot model ------
  ### Visualize the trajectory in the reduced dimensional space, coloring the cells by disease models
  models <- base::as.character(x = Biobase::pData(object = cds)$sample)
  base::cat("\n")
  base::print(x = base::table(models))
  for (m in base::names(x = select.rna.disease.samples)) {
    models[models %in% select.rna.disease.samples[[m]]] <- m
  }
  models[models %in% select.rna.control.samples] <- "Control"
  base::cat("\n")
  base::print(x = base::table(models))
  Biobase::pData(object = cds)$model <- base::as.character(x = models)
  base::cat("\n")
  base::print(x = utils::head(x = Biobase::pData(object = cds), n = 3L))
  base::print(x = utils::tail(x = Biobase::pData(object = cds), n = 3L))
  PlotMonocle2(cds = cds, feature.to.plot = "model", save.fig.dir = fig.mo2.dir, feature.name = "model", fn.suf = fn.suf, 
               plot.together = T, plot.split = T)
}

if (plot.disease) {
  ### Plot disease ------
  
  diseases <- base::as.character(x = Biobase::pData(object = cds)$sample)
  base::cat("\n")
  base::print(x = base::table(diseases))
  diseases[diseases %in% base::c("FA1", "FA2", "Notch3", "Notch4", "PGC1a-1", "PGC1a-2", "UUO1", "UUO2")] <- "CKD"
  diseases[diseases %in% base::c("longIRI14d1", "longIRI14d2", "shortIRI14d1", "shortIRI14d2")] <- "AKI"
  diseases[diseases %in% select.rna.control.samples] <- "Control"
  base::cat("\n")
  base::print(x = base::table(diseases))
  Biobase::pData(object = cds)$disease <- base::as.character(x = diseases)
  base::cat("\n")
  base::print(x = utils::head(x = Biobase::pData(object = cds), n = 3L))
  base::print(x = utils::tail(x = Biobase::pData(object = cds), n = 3L))
  PlotMonocle2(cds = cds, feature.to.plot = "disease", save.fig.dir = fig.mo2.dir, feature.name = "disease", fn.suf = fn.suf, 
               plot.together = T, plot.split = T)
}

if (plot.cycle) {
  ### Plot cell cycle ------
  ### Visualize the trajectory in the reduced dimensional space, coloring the cells by cell cycle
  PlotMonocle2(cds = cds, feature.to.plot = "Phase", save.fig.dir = fig.mo2.dir, feature.name = "cellcycle", fn.suf = fn.suf, 
               plot.together = T, plot.split = T)
}

if (plot.branch) {
  ### Plot branch ------
  BEAM_res <- base::readRDS(file = base::file.path(rds.dir, base::sprintf(fmt = "branch_dependent_genes%s.RDS", fn.suf)))
  hms <- base::list()
  for (bp in branch.points) {
    bp.name <- base::sprintf(fmt = "BP %s", bp)
    base::cat(base::sprintf(fmt = "\nPlotting Branch point %s ...\n", bp))
    tmp <- BEAM_res[[bp.name]]
    tmp <- tmp[1:100, base::c("gene_short_name", "pval", "qval")]
    sub.cds <- cds[base::row.names(x = tmp),]
    base::cat("\n")
    base::print(x = sub.cds)
    w <- set.w * base::ncol(x = sub.cds) + 6
    h <- set.h * base::nrow(x = sub.cds) + 3
    grDevices::pdf(file = base::file.path(fig.mo2.dir, base::sprintf(fmt = "branch%s_dependent_heatmap%s.pdf", bp, fn.suf)), 
                   width = w, height = h, onefile = T, paper = "special")
    hms[[bp.name]] <- monocle::plot_genes_branched_heatmap(cds_subset = sub.cds, 
                                                           branch_point = bp, 
                                                           cluster_rows = T, 
                                                           num_clusters = 2, 
                                                           show_rownames = T, 
                                                           use_gene_short_name = T, 
                                                           return_heatmap = T, 
                                                           cores = req.cores)
    grDevices::dev.off()
  }
  base::saveRDS(object = hms, file = base::file.path(rds.dir, base::sprintf(fmt = "branch_dependent_heatmaps%s.RDS", fn.suf)))
}

if (plot.lineage) {
  ### Plot lineage ------
  diff_test_res <- base::readRDS(file = base::file.path(rds.dir, base::sprintf(fmt = "pseudotime_dependent_genes%s.RDS", fn.suf)))
  hms <- base::list()
  for (l.name in base::names(x = lineages)) {
    base::cat(base::sprintf(fmt = "\nPlotting Lineage %s ...\n", l.name))
    lineage.cells <- base::row.names(x = subset(x = Biobase::pData(object = cds), subset = (State %in% lineages[[l.name]])))
    base::cat("\n")
    base::print(x = base::length(x = lineage.cells))
    tmp <- diff_test_res[[l.name]]
    tmp <- tmp[1:100, base::c("gene_short_name", "pval", "qval")]
    sub.cds <- cds[base::row.names(x = tmp), lineage.cells]
    base::cat("\n")
    base::print(x = sub.cds)
    base::cat("\n")
    base::print(x = base::unique(x = Biobase::pData(object = sub.cds)$State))
    w <- set.w * base::ncol(x = sub.cds) + 6
    h <- set.h * base::nrow(x = sub.cds) + 3
    grDevices::pdf(file = base::file.path(fig.mo2.dir, base::sprintf(fmt = "pseudotime_%s_dependent_heatmap%s.pdf", l.name, fn.suf)), 
                   width = w, height = h, onefile = T, paper = "special")
    hms[[l.name]] <- monocle::plot_pseudotime_heatmap(cds_subset = sub.cds, 
                                                      cluster_rows = T, 
                                                      num_clusters = 2, 
                                                      show_rownames = T, 
                                                      use_gene_short_name = T, 
                                                      return_heatmap = T, 
                                                      cores = req.cores)
    grDevices::dev.off()
  }
  base::saveRDS(object = hms, file = base::file.path(rds.dir, base::sprintf(fmt = "pseudotime_dependent_heatmaps%s.RDS", fn.suf)))
}

if (export.ddrtree) {
  ### Export DDRTree ------
  
  cell.ids <- base::as.character(x = Biobase::pData(object = cds)$cellID)
  states <- base::as.character(x = Biobase::pData(object = cds)$State)
  base::cat("\n")
  base::print(x = base::length(x = states))
  base::print(x = utils::head(x = states))
  base::names(x = states) <- cell.ids
  base::cat("\n")
  base::print(x = utils::head(x = states))
  base::saveRDS(object = states, file = base::file.path(rds.dir, "monocle2_State.RDS"))
  ddrtree.embed <- monocle::reducedDimS(cds = cds)
  base::cat("\n")
  base::print(x = base::dim(x = ddrtree.embed))
  ddrtree.embed <- base::t(x = ddrtree.embed)
  base::cat("\n")
  base::print(x = base::dim(x = ddrtree.embed))
  base::cat("\n")
  base::print(x = utils::head(x = ddrtree.embed))
  base::rownames(x = ddrtree.embed) <- cell.ids
  base::colnames(x = ddrtree.embed) <- base::paste0("Component_", 1:base::ncol(x = ddrtree.embed))
  base::cat("\n")
  base::print(x = utils::head(x = ddrtree.embed))
  
  base::cat("\n")
  base::print(x = orig.sobj)
  select.sobj <- subset(x = orig.sobj, subset = (cellID %in% cell.ids))
  base::cat("\n")
  base::print(x = select.sobj)
  cell.ids <- base::as.character(x = select.sobj$cellID)
  states <- states[cell.ids]
  base::cat("\n")
  base::print(x = base::length(x = states))
  base::print(x = utils::head(x = states))
  if (base::ncol(x = select.sobj) != base::sum(base::names(x = states) == cell.ids)) {
    base::stop("numbers of cells are not the same when exporting states.")
  }
  base::names(x = states) <- base::rownames(x = select.sobj@meta.data)
  select.sobj <- Seurat::AddMetaData(object = select.sobj, metadata = states, col.name = "State")
  base::cat("\n")
  base::print(x = base::table(select.sobj$State))
  ddrtree.embed <- ddrtree.embed[cell.ids,]
  base::cat("\n")
  base::print(x = base::dim(x = ddrtree.embed))
  base::cat("\n")
  base::print(x = utils::head(x = ddrtree.embed))
  if (base::ncol(x = select.sobj) != base::sum(base::rownames(x = ddrtree.embed) == cell.ids)) {
    base::stop("numbers of cells are not the same when exporting DDRTree embeddings.")
  }
  base::rownames(x = ddrtree.embed) <- base::rownames(x = select.sobj@meta.data)
  base::cat("\n")
  base::print(x = utils::head(x = ddrtree.embed))
  select.sobj[["ddrtree"]] <- Seurat::CreateDimReducObject(embeddings = ddrtree.embed, assay = act.assay, key = "Component_")
  base::cat("\n")
  base::print(x = select.sobj@reductions)
  base::saveRDS(object = select.sobj, file = base::file.path(rds.dir, "after_DDRTree.RDS"))
  
  if (dump.h5ad) {
    ## Remove the "scale.data" slot.
    Seurat::DefaultAssay(object = select.sobj) <- "RNA"
    select.sobj <- Seurat::SetAssayData(object = select.sobj, slot = "scale.data", new.data = methods::new(Class = "matrix"), assay = "RNA")
    base::cat("\n")
    base::print(x = select.sobj$RNA@scale.data)
    
    ## Add model info.
    models <- select.sobj$sample
    base::cat("\n")
    base::print(x = base::table(models))
    base::print(x = utils::head(x = models, n = 3L))
    base::print(x = utils::tail(x = models, n = 3L))
    cell.names <- base::names(x = models)
    base::cat("\n")
    base::print(x = utils::head(x = cell.names, n = 3L))
    base::print(x = utils::tail(x = cell.names, n = 3L))
    models <- base::as.character(x = models)
    base::cat("\n")
    base::print(x = base::table(models))
    base::print(x = utils::head(x = models, n = 3L))
    base::print(x = utils::tail(x = models, n = 3L))
    disease.models <- base::names(x = select.rna.disease.samples)
    for (dm in disease.models) {
      models[models %in% select.rna.disease.samples[[dm]]] <- dm
    }
    models[models %in% select.rna.control.samples] <- "Control"
    base::cat("\n")
    base::print(x = base::table(models))
    base::names(x = models) <- cell.names
    base::cat("\n")
    base::print(x = utils::head(x = models, n = 3L))
    base::print(x = utils::tail(x = models, n = 3L))
    select.sobj <- Seurat::AddMetaData(object = select.sobj, metadata = models, col.name = "model")
    base::cat("\n")
    base::print(x = utils::head(x = select.sobj[[]], n = 3L))
    base::print(x = utils::tail(x = select.sobj[[]], n = 3L))
    Seurat::Idents(object = select.sobj) <- "model"
    base::cat("\n")
    base::print(x = base::table(Seurat::Idents(object = select.sobj)))
    
    ## Add disease info.
    #diseases <- select.sobj$sample
    #base::cat("\n")
    #base::print(x = base::table(diseases))
    #base::print(x = utils::head(x = diseases, n = 3L))
    #base::print(x = utils::tail(x = diseases, n = 3L))
    #cell.names <- base::names(x = diseases)
    #base::cat("\n")
    #base::print(x = utils::head(x = cell.names, n = 3L))
    #base::print(x = utils::tail(x = cell.names, n = 3L))
    #diseases <- base::as.character(x = diseases)
    #base::cat("\n")
    #base::print(x = base::table(diseases))
    #base::print(x = utils::head(x = diseases, n = 3L))
    #base::print(x = utils::tail(x = diseases, n = 3L))
    #diseases[diseases %in% base::c("FA1", "FA2", "Notch3", "Notch4", "PGC1a-1", "PGC1a-2", "UUO1", "UUO2")] <- "CKD"
    #diseases[diseases %in% base::c("longIRI14d1", "longIRI14d2", "shortIRI14d1", "shortIRI14d2")] <- "AKI"
    #diseases[diseases %in% select.rna.control.samples] <- "Control"
    #base::cat("\n")
    #base::print(x = base::table(diseases))
    #base::names(x = diseases) <- cell.names
    #base::cat("\n")
    #base::print(x = utils::head(x = diseases, n = 3L))
    #base::print(x = utils::tail(x = diseases, n = 3L))
    #select.sobj <- Seurat::AddMetaData(object = select.sobj, metadata = diseases, col.name = "disease")
    #base::cat("\n")
    #base::print(x = utils::head(x = select.sobj[[]], n = 3L))
    #base::print(x = utils::tail(x = select.sobj[[]], n = 3L))
    #Seurat::Idents(object = select.sobj) <- "disease"
    #base::cat("\n")
    #base::print(x = base::table(Seurat::Idents(object = select.sobj)))
    
    ## Add condition info.
    conditions <- select.sobj$sample
    base::cat("\n")
    base::print(x = base::table(conditions))
    base::print(x = utils::head(x = conditions, n = 3L))
    base::print(x = utils::tail(x = conditions, n = 3L))
    cell.names <- base::names(x = conditions)
    base::cat("\n")
    base::print(x = utils::head(x = cell.names, n = 3L))
    base::print(x = utils::tail(x = cell.names, n = 3L))
    conditions <- base::as.character(x = conditions)
    base::cat("\n")
    base::print(x = base::table(conditions))
    base::print(x = utils::head(x = conditions, n = 3L))
    base::print(x = utils::tail(x = conditions, n = 3L))
    conditions[conditions %in% rna.disease.samples] <- "diseased"
    conditions[conditions %in% rna.control.samples] <- "healthy"
    base::cat("\n")
    base::print(x = base::table(conditions))
    base::names(x = conditions) <- cell.names
    base::cat("\n")
    base::print(x = utils::head(x = conditions, n = 3L))
    base::print(x = utils::tail(x = conditions, n = 3L))
    select.sobj <- Seurat::AddMetaData(object = select.sobj, metadata = conditions, col.name = "condition")
    base::cat("\n")
    base::print(x = utils::head(x = select.sobj[[]], n = 3L))
    base::print(x = utils::tail(x = select.sobj[[]], n = 3L))
    Seurat::Idents(object = select.sobj) <- "condition"
    base::cat("\n")
    base::print(x = base::table(Seurat::Idents(object = select.sobj)))
    ## Convert Seurat into AnnData(H5AD).
    trans.file = base::file.path(rds.dir, "after_DDRTree.h5Seurat")
    SeuratDisk::SaveH5Seurat(object = select.sobj, filename = trans.file, overwrite = T, verbose = T)
    SeuratDisk::Convert(source = trans.file, dest = "h5ad", assay = "RNA", overwrite = T, verbose = T)
    if (base::file.exists(base::file.path(rds.dir, "after_DDRTree.h5ad"))) {
      base::file.remove(trans.file)
    }
  }
}
