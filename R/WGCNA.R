# WGCNA.R
base::cat("\nRunning \"WGCNA.R\" ...\n")

### Mouse scRNA PT ------

sobj <- base::readRDS(file = base::file.path(rds.dir, "after_DDRTree.RDS"))
curr.sobj <- sobj

curr.sobj[["spliced"]] <- NULL
curr.sobj[["unspliced"]] <- NULL

all.models <- base::unique(x = base::as.character(x = curr.sobj$model.1))
all.celltypes <- base::unique(x = base::as.character(x = curr.sobj$cell.type))
all.models
all.celltypes

curr.sobj$model.cell.type <- base::paste(
  base::as.character(x = curr.sobj$model.1), 
  base::as.character(x = curr.sobj$cell.type), 
  sep = "_"
)

nfeatures <- 4000
curr.sobj <- Seurat::FindVariableFeatures(object = curr.sobj, assay = "RNA", selection.method = "vst",
                                          nfeatures = nfeatures, verbose = F)
curr.sobj
head(curr.sobj[[]], n=3L)

HVGs <- Seurat::VariableFeatures(object = curr.sobj, assay = "RNA")
length(HVGs)
head(HVGs)

k <- 10
min.cells <- 20
split.sobj <- Seurat::SplitObject(object = curr.sobj, split.by = "model.cell.type")
select.sobj <- base::list()
for (curr.model in all.models) {
  for (curr.celltype in all.celltypes) {
    curr.group <- base::paste(curr.model, curr.celltype, sep = "_")
    if (!utils::hasName(x = split.sobj, name = curr.group)) {
      base::cat(base::sprintf(fmt = "\n%s does not exist.\n", curr.group))
      next
    }
    curr.sobj <- split.sobj[[curr.group]]
    if (base::ncol(x = curr.sobj) < min.cells) {
      base::cat(base::sprintf(fmt = "\n%s has %s cells.\nSkip %s.\n", curr.group, base::ncol(x = curr.sobj), curr.group))
      next
    }
    base::cat(base::sprintf(fmt = "\nProcessing %s ...\n", curr.group))
    curr.sobj <- curr.sobj[HVGs,]
    curr.metacell.sobj <- hdWGCNA::construct_metacells(seurat_obj = curr.sobj, name = curr.group, k = k, 
                                                       reduction = "umap", assay = "RNA", slot = "data")
    curr.metacell.sobj$model <- curr.model
    curr.metacell.sobj$cell.type <- curr.celltype
    select.sobj[[curr.group]] <- curr.metacell.sobj
  }
}
metacell.sobj <- merge(x = select.sobj[[1]], y = select.sobj[2:base::length(x = select.sobj)], 
                       add.cell.ids = NULL, merge.data = T)
metacell.sobj
head(x = metacell.sobj[[]], n = 3L)
tail(x = metacell.sobj[[]], n = 3L)
table(metacell.sobj$model)
table(metacell.sobj$cell.type)
base::saveRDS(
  object = metacell.sobj, 
  file = base::file.path(rds.dir, base::sprintf(fmt = "metacell_model_celltype_%sHVGs_k%s.RDS", length(HVGs), k))
)
split.sobj <- NULL
select.sobj <- NULL
curr.sobj <- NULL
base::gc()

orig.metacell.sobj <- metacell.sobj

metacell.sobj <- orig.metacell.sobj
metacell.sobj <- Seurat::ScaleData(object = metacell.sobj, assay = "RNA", 
                                   features = base::row.names(x = metacell.sobj), verbose = F)
npcs <- 51
metacell.sobj <- Seurat::RunPCA(object = metacell.sobj, assay = "RNA", 
                                features = base::row.names(x = metacell.sobj), npcs = npcs, verbose = F)
metadata <- metacell.sobj@meta.data
pca.embed <- Seurat::Embeddings(object = metacell.sobj, reduction = "pca")
reticulate::use_python(python = "/home/jzhou88/installed/miniconda3/envs/r4/bin/python")
harmonypy <- reticulate::import(module = "harmonypy")
ho <- harmonypy$run_harmony(pca.embed, metadata, c("model"), max_iter_harmony = 100L, verbose = T)
the <- ho$Z_corr
ho <- NULL
base::gc()
harmony.embed <- base::t(x = the)
base::rownames(x = harmony.embed) <- base::rownames(x = pca.embed)
base::colnames(x = harmony.embed) <- base::paste0("harmony_", 1:npcs)
metacell.sobj[["harmony"]] <- Seurat::CreateDimReducObject(embeddings = harmony.embed, assay = "RNA", 
                                                           key = "harmony_")
metacell.sobj <- Seurat::RunUMAP(object = metacell.sobj, dims = 1:npcs, reduction = "harmony", 
                                 umap.method = "uwot", return.model = T, metric = "cosine", 
                                 verbose = F, reduction.name = "umap", reduction.key = "UMAP_")

known.models <- base::c("Control", "FAN", "UUO",
                        "longIRI1d", "longIRI3d", "longIRI14d",
                        "shortIRI1d", "shortIRI3d", "shortIRI14d",
                        "Notch1", "PGC1a", "APOL1", "Esrra",
                        "LPS1hr", "LPS4hr", "LPS16hr", "LPS27hr", "LPS36hr", "LPS48hr")
metacell.sobj$model <- base::factor(x = metacell.sobj$model, levels = known.models)
Seurat::Idents(object = metacell.sobj) <- "model"
w <- 6
h <- 5
dp <- Seurat::DimPlot(object = metacell.sobj, reduction = "umap", label = F) +
  Seurat::NoAxes()
grDevices::pdf(file = base::file.path(fig.wgc.dir,
                                      base::sprintf(fmt = "umap_metacell_model_celltype_%sHVGs_k%s.pdf", length(HVGs), k)), 
               width = w, height = h, onefile = T, paper = "special")
base::print(x = dp)
grDevices::dev.off()

known.celltypes <- base::c("S1", "S2", "S2/S3", "S3", "Injured1", "Injured2")
metacell.sobj$cell.type <- base::factor(x = metacell.sobj$cell.type, levels = known.celltypes)
Seurat::Idents(object = metacell.sobj) <- "cell.type"
w <- 6
h <- 5
dp <- Seurat::DimPlot(object = metacell.sobj, reduction = "umap", label = F) + 
  Seurat::NoAxes()
grDevices::pdf(
  file = base::file.path(fig.wgc.dir,
                         base::sprintf(fmt = "umap_metacell_model_celltype_%sHVGs_k%s_GroupedByCellTypes.pdf", length(HVGs), k)),
  width = w, height = h, onefile = T, paper = "special"
)
base::print(x = dp)
grDevices::dev.off()

metacell.sobj <- orig.metacell.sobj
base::options(stringsAsFactors = FALSE)
WGCNA::enableWGCNAThreads()
datExpr <- base::as.data.frame(x = Seurat::GetAssayData(object = metacell.sobj, slot = "data", assay = "RNA"))
datExpr <- base::as.data.frame(x = base::t(x = datExpr))
datExpr <- datExpr[, WGCNA::goodGenes(datExpr = datExpr)]
base::dim(x = datExpr)

powers <- base::c(base::seq(from = 1, to = 10, by = 1), base::seq(from = 12, to = 30, by = 2))
powerTable <- base::list(
  data = WGCNA::pickSoftThreshold(
    data = datExpr, 
    powerVector = powers, 
    blockSize = 20000,
    corFnc = "bicor", 
    networkType = "signed", 
    verbose = 2
  )[[2]]
)

colors <- base::c("red")
plotCols <- base::c(2, 5, 6, 7)
colNames <- base::c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity")
ylim <- base::matrix(data = NA, nrow = 2, ncol = 4)
for (col in 1:base::length(x = plotCols)) {
  ylim[1, col] <- base::min(ylim[1, col], powerTable$data[, plotCols[col]], na.rm = TRUE)
  ylim[2, col] <- base::max(ylim[2, col], powerTable$data[, plotCols[col]], na.rm = TRUE)
}
par(mfcol = base::c(2,2))
par(mar = base::c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7
for (col in 1:base::length(x = plotCols)){
  plot(powerTable$data[,1], -sign(powerTable$data[,3])*powerTable$data[,2], xlab="Soft Threshold (power)",
       ylab=colNames[col],type="n", ylim = ylim[, col], main = colNames[col])
  addGrid()
  
  if (col==1) {
    text(powerTable$data[,1], -sign(powerTable$data[,3])*powerTable$data[,2], labels=powers,cex=cex1,col=colors[1])
  } else {
    text(powerTable$data[,1], powerTable$data[,plotCols[col]], labels=powers,cex=cex1,col=colors[1])
  }
}

multiExpr <- list()
multiExpr[["PT"]] <- list(data=datExpr)
checkSets(multiExpr)
head(multiExpr$PT$data)

cor.type <- "pearson"
min.module.size <- 50

softPower <- 5
net <- WGCNA::blockwiseConsensusModules(multiExpr, 
                                        blocks = NULL, 
                                        maxBlockSize = 20000, ## This should be set to a smaller size if the user has limited RAM
                                        corType = cor.type, ## no use for bicor
                                        power = softPower,
                                        consensusQuantile = 0.3,
                                        networkType = "signed",
                                        TOMType = "unsigned",
                                        TOMDenom = "min",
                                        scaleTOMs = TRUE, scaleQuantile = 0.8,
                                        sampleForScaling = TRUE, sampleForScalingFactor = 1000,
                                        useDiskCache = T, chunkSize = NULL,
                                        deepSplit = 4,
                                        pamStage = F,
                                        detectCutHeight = 0.995, minModuleSize = min.module.size,
                                        mergeCutHeight = 0.2,
                                        saveConsensusTOMs = F,
                                        consensusTOMFilePattern = "ConsensusTOM-block.%b.rda", 
                                        numericLabels = F, 
                                        nThreads = 8)
base::saveRDS(
  object = net, 
  file = base::file.path(
    rds.dir, 
    base::sprintf(fmt = "WGCNA_consensus_modules_%s_MinModSize%s.RDS", cor.type, min.module.size)
  )
)
table(net$colors)

consTree <- net$dendrograms[[1]]
moduleLabels <- net$colors;
moduleColors <- base::as.character(x = moduleLabels)
w <- 8
h <- 5
grDevices::pdf(
  file = base::file.path(
    fig.wgc.dir, 
    base::sprintf(fmt = "WGCNA_consensus_gene_dendrogram_%s_MinModSize%s.pdf", cor.type, min.module.size)
  ), width = w, height = h, onefile = T, paper = "special"
)
plotDendroAndColors(dendro = consTree, 
                    colors = moduleColors, 
                    groupLabels = "Module", 
                    autoColorHeight = T, 
                    dendroLabels = F, 
                    hang = 0.03, 
                    addGuide = T, 
                    guideHang = 0.05, 
                    main = "")
grDevices::dev.off()

MEs <- WGCNA::moduleEigengenes(expr = multiExpr[[1]]$data, colors = moduleColors, nPC = 1)$eigengenes
MEs <- WGCNA::orderMEs(MEs = MEs)
KMEs <- WGCNA::signedKME(datExpr = datExpr, datME = MEs, outputColumnName = "kME", corFnc = "bicor")
geneInfo <- base::as.data.frame(x = base::cbind(base::colnames(datExpr), moduleLabels, moduleColors, KMEs))
base::colnames(geneInfo)[1] <- "GeneSymbol"
base::colnames(geneInfo)[2] <- "ModuleLabel"
base::colnames(geneInfo)[3] <- "Initially.Assigned.Module.Color"
utils::write.csv(
  x = geneInfo,
  file = base::file.path(
    rst.wgc.dir, 
    base::sprintf(fmt = "WGCNA_geneInfoSigned_%s_MinModSize%s.csv", cor.type, min.module.size)
  )
)

n_genes <- 25
module_labels <- unique(geneInfo$ModuleLabel)
modules <- unique(geneInfo$Initially.Assigned.Module.Color)
module_labels <- module_labels[order(modules)]
modules <- modules[order(modules)]
mod_colors <- modules
names(mod_colors) <- module_labels
module_labels
modules
mod_colors
module_list <- lapply(modules, function(mod){
  cur <- subset(geneInfo, Initially.Assigned.Module.Color == mod)
  cur[,c('GeneSymbol', paste0('kME', mod))] %>%
    top_n(n_genes) %>% .$GeneSymbol
})
names(module_list) <- module_labels
module_list

test.sobj <- curr.sobj.bef.wgcna
test.sobj <- Seurat::AddModuleScore(object = test.sobj, features = module_list, pool = base::row.names(x = test.sobj), 
                                    nbin = 24, k = F, assay = "RNA", name = "WGCNA_module")
head(test.sobj[[]])

# get module scores from scRNA data
cur_celltype <- "WGCNA"
module_df <- test.sobj@meta.data
features <- names(module_df)[grepl(cur_celltype, names(module_df))]
feat2mod <- module_labels
names(feat2mod) <- features
feat2mod
module_df <- module_df[,c('model.1', features)]
head(module_df)

# compute average module score for each cluster:
tmp <- lapply(unique(module_df$model.1), function(i){
  cur_df <- module_df %>% subset(model.1 == i)
  data.frame(
    value=as.numeric(colSums(cur_df[,features]) / nrow(cur_df)),
    cluster = i,
    feature = features
  )
})
plot_df <- Reduce(rbind, tmp)
plot_df$feature <- feat2mod[plot_df$feature]
tmp
plot_df

# remove gray module:
plot_df <- subset(plot_df, feature != 'grey')
plot_df
mod_colors <- mod_colors[module_labels != 'grey']
module_labels <- module_labels[module_labels != 'grey']
module_labels
names(mod_colors) <- module_labels
mod_colors
# reshape to wide format
plot_df <- reshape2::dcast(plot_df, feature ~ cluster)
rownames(plot_df) <- plot_df$feature
plot_df
plot_df <- plot_df %>% dplyr::select(-c(feature))

module.order <- base::c("black", "blue", "brown", "green", "magenta", "pink", "red", "turquoise", "yellow")
known.models <- base::c("Control", "FAN", "UUO",
                        "longIRI1d", "longIRI3d", "longIRI14d", "shortIRI1d", "shortIRI3d", "shortIRI14d",
                        "Notch1", "PGC1a", "APOL1", "Esrra",
                        "LPS1hr", "LPS4hr", "LPS16hr", "LPS27hr", "LPS36hr", "LPS48hr")
model.order <- base::c("Control",
                       "LPS1hr", "LPS4hr", "LPS16hr", "LPS27hr", "LPS36hr", "LPS48hr",
                       "FAN", "UUO", "PGC1a", "Notch1", "Esrra", "APOL1",
                       "longIRI1d", "longIRI3d", "longIRI14d", "shortIRI1d", "shortIRI3d", "shortIRI14d")

zScore <- function(x){(x - mean(x)) /sd(x)}
matrix_z <- apply(plot_df, 1, zScore)
matrix_z <- matrix_z[model.order,]
#matrix_z
# keep full values to plot onto the heatmap
matrix_full <- matrix_z
#matrix_anno <- ifelse(matrix_full >= quantile(matrix_full, 0.80), signif(matrix_full,2), '')
matrix_anno <- ifelse(matrix_full >= quantile(matrix_full, 0.75), signif(matrix_full,2), '')
#matrix_anno
# change the range
range(matrix_z)
matrix_z <- ifelse(matrix_z >= 2, 2, matrix_z)
matrix_z <- ifelse(matrix_z <= -2, -2, matrix_z)
#matrix_z
# rename modules
colnames(matrix_z) <- module_labels
colnames(matrix_anno) <- module_labels
matrix_anno <- matrix_anno[,module.order]
matrix_z <- matrix_z[,module.order]
# color the rows by scWGCNA module color, cols by single-cell clusters
cluster_colors <- scales::hue_pal()(base::length(x = known.models))
base::names(x = cluster_colors) <- known.models
cluster_colors <- cluster_colors[rownames(matrix_z)]

column_ha <- HeatmapAnnotation(
  Module = colnames(matrix_z),
  col = list(
    Module = mod_colors
  )
)
row_ha <- rowAnnotation(
  Model = rownames(matrix_z),
  col = list(
    Model = cluster_colors
  )
)

w <- 16
h <- 8
w.1 <- unit(11, "in") 
h.1 <- unit(6, "in")
cluster.rows <- F
cluster.cols <- T
grDevices::pdf(
  file = base::file.path(
    fig.wgc.dir, 
    base::sprintf(fmt = "WGCNA_module_score_heatmap_%s_MinModSize%s_transposed.pdf", cor.type, min.module.size)
  ), width = w, height = h, onefile = T, paper = "special"
)
ComplexHeatmap::Heatmap(
  matrix_z,
  cluster_rows = cluster.rows,
  cluster_columns = cluster.cols,
  top_annotation = column_ha,
  left_annotation = row_ha,
  width = w.1, 
  height = h.1,
  column_names_rot = 45,
  show_heatmap_legend = T,
  use_raster = T,
  cell_fun = function(j,i,x,y,w,h,col){
    grid.text(matrix_anno[i,j], x, y)
  }
)
grDevices::dev.off()

module.labels <- base::unique(x = geneInfo$ModuleLabel)
modules <- base::list()
for (curr.label in module.labels) {
  base::cat(base::sprintf(fmt = "\nProcessing Module %s ...\n", curr.label))
  modules[[curr.label]] <- geneInfo[geneInfo$ModuleLabel == curr.label,]
}
openxlsx::write.xlsx(
  x = modules, 
  file = base::file.path(
    rst.wgc.dir, base::sprintf(fmt = "WGCNA_SignedGeneInfo_PerModule_%s_MinModSize%s.xlsx", cor.type, min.module.size)
  ), startCol = 1, startRow = 1, colNames = T, rowNames = F, overwrite = T
)

max.p.val <- 0.05
fn <- "WGCNA_Module_KEGG_GOTERM_BP"
gene.sets <- Xlsx2List(
  xlsx.file = base::file.path(rst.wgc.dir, base::sprintf(fmt = "WebGestalt/%s.xlsx", fn)), 
  rowNames = F, cols = base::c(1, 2, 7, 9, 11)
)

all.modules <- base::c("black", "blue", "brown", "green", "magenta", "pink", "red", "turquoise", "yellow")
bar.colors <- base::c("white", "blue", "brown", "green", "magenta", "pink", "red", "turquoise", "yellow")
text.colors <- base::c("black", "white", "white", "black", "black", "black", "black", "black", "black")
base::names(x = bar.colors) <- all.modules
base::names(x = text.colors) <- all.modules

top.n <- 3
f <- 16
w <- 16
h <- 2
for (curr.module in all.modules) {
  base::cat(base::sprintf(fmt = "\nProcessing %s ...\n", curr.module))
  curr.gs <- gene.sets[[curr.module]]
  curr.gs <- curr.gs[base::which(curr.gs$database == "pathway_KEGG"),]
  curr.gs <- curr.gs[base::which(curr.gs$FDR < max.p.val),]
  curr.gs <- curr.gs[base::order(curr.gs$enrichmentRatio, decreasing = T),]
  top.n.gs <- curr.gs
  if (base::nrow(x = curr.gs) > top.n) {
    top.n.gs <- curr.gs[1:top.n,]
  }
  
  # Create the bar plot using ggplot2
  p <- ggplot2::ggplot(
    data = top.n.gs,
    mapping = ggplot2::aes(x = stats::reorder(description, enrichmentRatio), y = enrichmentRatio)
  ) +
    ggplot2::geom_bar(stat = "identity", fill = bar.colors[curr.module]) +
    ggplot2::coord_flip() + # Flip the x and y axes
    ggplot2::labs(x = NULL, y = NULL) + # Set the x and y axis labels
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    ) + 
    ggplot2::geom_text(
      mapping = ggplot2::aes(label = description, y = 0),
      color = text.colors[curr.module],
      fontface = "bold",
      hjust = 0,
      size = 9)
  # Display the plot
  grDevices::pdf(
    file = base::file.path(
      fig.pub.dir, 
      base::sprintf(fmt = "barplot_WGCNA_Module_Top%sKEGG_%s.pdf", top.n, curr.module)
    ), width = w, height = h, onefile = T, paper = "special")
  base::print(x = p)
  grDevices::dev.off()
}

### Mouse DKD snRNA PT ------

sobj <- base::readRDS(file = base::file.path(rds.dir, "after_neighbor.RDS"))
curr.sobj <- sobj

all.conditions <- base::unique(x = base::as.character(x = curr.sobj$condition))
all.celltypes <- base::unique(x = base::as.character(x = curr.sobj$cell.type))
all.conditions
all.celltypes

curr.sobj$condition.cell.type <- base::paste(
  base::as.character(x = curr.sobj$condition), 
  base::as.character(x = curr.sobj$cell.type), 
  sep = "_"
)
table(curr.sobj$condition.cell.type)
sum(table(curr.sobj$condition.cell.type))

nfeatures <- 4000
curr.sobj <- Seurat::FindVariableFeatures(object = curr.sobj, assay = "RNA", selection.method = "vst",
                                          nfeatures = nfeatures, verbose = F)
curr.sobj
head(curr.sobj[[]], n=3L)

HVGs <- Seurat::VariableFeatures(object = curr.sobj, assay = "RNA")
length(HVGs)
head(HVGs)

k <- 10
min.cells <- 20
split.sobj <- Seurat::SplitObject(object = curr.sobj, split.by = "condition.cell.type")
select.sobj <- base::list()
for (curr.condition in all.conditions) {
  for (curr.celltype in all.celltypes) {
    curr.group <- base::paste(curr.condition, curr.celltype, sep = "_")
    if (!utils::hasName(x = split.sobj, name = curr.group)) {
      base::cat(base::sprintf(fmt = "\n%s does not exist.\n", curr.group))
      next
    }
    curr.sobj <- split.sobj[[curr.group]]
    if (base::ncol(x = curr.sobj) < min.cells) {
      base::cat(base::sprintf(fmt = "\n%s has %s cells.\nSkip %s.\n", curr.group, base::ncol(x = curr.sobj), curr.group))
      next
    }
    base::cat(base::sprintf(fmt = "\nProcessing %s ...\n", curr.group))
    curr.sobj <- curr.sobj[HVGs,]
    curr.metacell.sobj <- hdWGCNA::construct_metacells(seurat_obj = curr.sobj, name = curr.group, k = k, 
                                                       reduction = "umap", assay = "RNA", slot = "data")
    curr.metacell.sobj$condition <- curr.condition
    curr.metacell.sobj$cell.type <- curr.celltype
    curr.metacell.sobj$condition.cell.type <- curr.group
    select.sobj[[curr.group]] <- curr.metacell.sobj
  }
}
metacell.sobj <- merge(x = select.sobj[[1]], y = select.sobj[2:base::length(x = select.sobj)], 
                       add.cell.ids = NULL, merge.data = T)
metacell.sobj
head(x = metacell.sobj[[]], n = 3L)
tail(x = metacell.sobj[[]], n = 3L)
table(metacell.sobj$condition)
table(metacell.sobj$cell.type)
table(metacell.sobj$condition.cell.type)
base::saveRDS(
  object = metacell.sobj, 
  file = base::file.path(rds.dir, base::sprintf(fmt = "metacell_condition_celltype_%sHVGs_k%s.RDS", length(HVGs), k))
)
split.sobj <- NULL
select.sobj <- NULL
curr.sobj <- NULL
base::gc()

orig.metacell.sobj <- metacell.sobj

metacell.sobj <- orig.metacell.sobj
metacell.sobj <- Seurat::ScaleData(object = metacell.sobj, assay = "RNA", 
                                   features = base::row.names(x = metacell.sobj), verbose = F)
npcs <- 51
metacell.sobj <- Seurat::RunPCA(object = metacell.sobj, assay = "RNA", 
                                features = base::row.names(x = metacell.sobj), npcs = npcs, verbose = F)
metadata <- metacell.sobj@meta.data
pca.embed <- Seurat::Embeddings(object = metacell.sobj, reduction = "pca")
reticulate::use_python(python = "/home/jzhou88/installed/miniconda3/envs/r4/bin/python")
harmonypy <- reticulate::import(module = "harmonypy")
ho <- harmonypy$run_harmony(pca.embed, metadata, c("condition"), max_iter_harmony = 100L, verbose = T)
the <- ho$Z_corr
ho <- NULL
base::gc()
harmony.embed <- base::t(x = the)
base::rownames(x = harmony.embed) <- base::rownames(x = pca.embed)
base::colnames(x = harmony.embed) <- base::paste0("harmony_", 1:npcs)
metacell.sobj[["harmony"]] <- Seurat::CreateDimReducObject(embeddings = harmony.embed, assay = "RNA", 
                                                           key = "harmony_")
metacell.sobj <- Seurat::RunUMAP(object = metacell.sobj, dims = 1:npcs, reduction = "harmony", 
                                 umap.method = "uwot", return.model = T, metric = "cosine", 
                                 verbose = F, reduction.name = "umap", reduction.key = "UMAP_")

metacell.sobj$condition <- base::factor(x = metacell.sobj$condition, levels = base::c("Control", "DKD"))
Seurat::Idents(object = metacell.sobj) <- "condition"
w <- 6
h <- 5
dp <- Seurat::DimPlot(object = metacell.sobj, reduction = "umap", label = F) + 
  Seurat::NoAxes()
grDevices::pdf(file = base::file.path(fig.wgc.dir,
                                      base::sprintf(fmt = "umap_metacell_condition_celltype_%sHVGs_k%s.pdf", length(HVGs), k)), 
               width = w, height = h, onefile = T, paper = "special")
base::print(x = dp)
grDevices::dev.off()

metacell.sobj$cell.type <- base::factor(x = metacell.sobj$cell.type, levels = base::c("PT", "injured PT"))
Seurat::Idents(object = metacell.sobj) <- "cell.type"
w <- 6
h <- 5
dp <- Seurat::DimPlot(object = metacell.sobj, reduction = "umap", label = F) + 
  Seurat::NoAxes()
grDevices::pdf(
  file = base::file.path(fig.wgc.dir,
                         base::sprintf(fmt = "umap_metacell_condition_celltype_%sHVGs_k%s_GroupedByCellTypes.pdf", length(HVGs), k)),
  width = w, height = h, onefile = T, paper = "special"
)
base::print(x = dp)
grDevices::dev.off()

metacell.sobj <- orig.metacell.sobj
base::options(stringsAsFactors = FALSE)
WGCNA::enableWGCNAThreads()
datExpr <- base::as.data.frame(x = Seurat::GetAssayData(object = metacell.sobj, slot = "data", assay = "RNA"))
datExpr <- base::as.data.frame(x = base::t(x = datExpr))
datExpr <- datExpr[, WGCNA::goodGenes(datExpr = datExpr)]
base::dim(x = datExpr)

powers <- base::c(base::seq(from = 1, to = 10, by = 1), base::seq(from = 12, to = 30, by = 2))
powerTable <- base::list(
  data = WGCNA::pickSoftThreshold(
    data = datExpr, 
    powerVector = powers, 
    blockSize = 20000,
    corFnc = "bicor", 
    networkType = "signed", 
    verbose = 2
  )[[2]]
)

colors <- base::c("red")
plotCols <- base::c(2, 5, 6, 7)
colNames <- base::c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity")
ylim <- base::matrix(data = NA, nrow = 2, ncol = 4)
for (col in 1:base::length(x = plotCols)) {
  ylim[1, col] <- base::min(ylim[1, col], powerTable$data[, plotCols[col]], na.rm = TRUE)
  ylim[2, col] <- base::max(ylim[2, col], powerTable$data[, plotCols[col]], na.rm = TRUE)
}
par(mfcol = base::c(2,2))
par(mar = base::c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7
for (col in 1:base::length(x = plotCols)){
  plot(powerTable$data[,1], -sign(powerTable$data[,3])*powerTable$data[,2], xlab="Soft Threshold (power)", ylab=colNames[col],type="n", 
       ylim = ylim[, col], main = colNames[col])
  addGrid()
  
  if (col==1) {
    text(powerTable$data[,1], -sign(powerTable$data[,3])*powerTable$data[,2], labels=powers,cex=cex1,col=colors[1])
  } else {
    text(powerTable$data[,1], powerTable$data[,plotCols[col]], labels=powers,cex=cex1,col=colors[1])
  }
}

multiExpr <- list()
multiExpr[["PT"]] <- list(data=datExpr)
checkSets(multiExpr)
head(multiExpr$PT$data)

cor.type <- "pearson"
min.module.size <- 50

softPower <- 7
net <- WGCNA::blockwiseConsensusModules(multiExpr, 
                                        blocks = NULL, 
                                        maxBlockSize = 20000, ## This should be set to a smaller size if the user has limited RAM
                                        corType = cor.type, ## no use for bicor
                                        power = softPower,
                                        consensusQuantile = 0.3,
                                        networkType = "signed",
                                        TOMType = "unsigned",
                                        TOMDenom = "min",
                                        scaleTOMs = TRUE, scaleQuantile = 0.8,
                                        sampleForScaling = TRUE, sampleForScalingFactor = 1000,
                                        useDiskCache = T, chunkSize = NULL,
                                        deepSplit = 4,
                                        pamStage = F,
                                        detectCutHeight = 0.995, minModuleSize = min.module.size,
                                        mergeCutHeight = 0.2,
                                        saveConsensusTOMs = F,
                                        consensusTOMFilePattern = "ConsensusTOM-block.%b.rda", 
                                        numericLabels = F, 
                                        nThreads = 8)
base::saveRDS(
  object = net, 
  file = base::file.path(
    rds.dir, 
    base::sprintf(fmt = "WGCNA_consensus_modules_%s_MinModSize%s.RDS", cor.type, min.module.size)
  )
)
table(net$colors)

consTree <- net$dendrograms[[1]]
moduleLabels <- net$colors;
moduleColors <- base::as.character(x = moduleLabels)
w <- 8
h <- 5
grDevices::pdf(
  file = base::file.path(
    fig.wgc.dir, 
    base::sprintf(fmt = "WGCNA_consensus_gene_dendrogram_%s_MinModSize%s.pdf", cor.type, min.module.size)
  ), width = w, height = h, onefile = T, paper = "special"
)
plotDendroAndColors(dendro = consTree, 
                    colors = moduleColors, 
                    groupLabels = "Module", 
                    autoColorHeight = T, 
                    dendroLabels = F, 
                    hang = 0.03, 
                    addGuide = T, 
                    guideHang = 0.05, 
                    main = "")
grDevices::dev.off()

MEs <- WGCNA::moduleEigengenes(expr = multiExpr[[1]]$data, colors = moduleColors, nPC = 1)$eigengenes
MEs <- WGCNA::orderMEs(MEs = MEs)
KMEs <- WGCNA::signedKME(datExpr = datExpr, datME = MEs, outputColumnName = "kME", corFnc = "bicor")
geneInfo <- base::as.data.frame(x = base::cbind(base::colnames(datExpr), moduleLabels, moduleColors, KMEs))
base::colnames(geneInfo)[1] <- "GeneSymbol"
base::colnames(geneInfo)[2] <- "ModuleLabel"
base::colnames(geneInfo)[3] <- "Initially.Assigned.Module.Color"
utils::write.csv(
  x = geneInfo,
  file = base::file.path(
    rst.wgc.dir, 
    base::sprintf(fmt = "WGCNA_geneInfoSigned_%s_MinModSize%s.csv", cor.type, min.module.size)
  )
)

module.labels <- base::unique(x = geneInfo$ModuleLabel)
modules <- base::list()
for (curr.label in module.labels) {
  base::cat(base::sprintf(fmt = "\nProcessing Module %s ...\n", curr.label))
  modules[[curr.label]] <- geneInfo[geneInfo$ModuleLabel == curr.label,]
}
openxlsx::write.xlsx(
  x = modules, 
  file = base::file.path(
    rst.wgc.dir, base::sprintf(fmt = "WGCNA_SignedGeneInfo_PerModule_%s_MinModSize%s.xlsx", cor.type, min.module.size)
  ), startCol = 1, startRow = 1, colNames = T, rowNames = F, overwrite = T
)

mouse.datExpr <- datExpr
mouse.colors <- net$colors

human.sobj <- base::readRDS(
  file = "/home/jzhou88/projects/KidneyDisease/2022-04-21/HumanDkdRnaNucleusPtSeuratHarmonypy/RDS/after_neighbor.RDS"
)
human.datExpr <- base::as.data.frame(x = Seurat::GetAssayData(object = human.sobj, slot = "data", assay = "RNA"))
human.datExpr <- base::as.data.frame(x = base::t(x = human.datExpr))
human.datExpr <- human.datExpr[, WGCNA::goodGenes(datExpr = human.datExpr)]
base::dim(x = human.datExpr)

## Converting mouse genes into human genes does not serve the purpose, 
## because different mouse genes can be mapped into the same human gene.
h2m.genes <- nichenetr::convert_human_to_mouse_symbols(
  symbols = base::colnames(x = human.datExpr), 
  version = 1
)

base::print(x = base::length(x = h2m.genes))
base::print(x = base::sum(base::is.na(x = h2m.genes)))
base::print(x = utils::head(x = h2m.genes))

h2m.genes[base::is.na(x = h2m.genes)] <- base::colnames(x = human.datExpr)[base::is.na(x = h2m.genes)]

base::print(x = base::length(x = h2m.genes))
base::print(x = base::sum(base::is.na(x = h2m.genes)))
base::print(x = utils::head(x = h2m.genes))

base::colnames(x = human.datExpr) <- h2m.genes
head(colnames(human.datExpr))

setLabels <- c("Mouse", "Human")
multiExpr <- base::list(Mouse = base::list(data = mouse.datExpr), Human = base::list(data = human.datExpr))
multiColor <- list(Mouse = mouse.colors)

mp <- WGCNA::modulePreservation(
  multiData = multiExpr,
  multiColor = multiColor,
  referenceNetworks = 1,
  nPermutations = 200,
  randomSeed = 1,
  quickCor = 0,
  verbose = 3
)
base::saveRDS(object = mp, file = base::file.path(rds.dir, "WGCNA_module_preservation_WholeHumanGenes.RDS"))

ref <- 1
test <- 2
statsObs <- base::cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ <-base::cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])

base::print(x = base::cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
                            signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

# Module labels and module sizes are also contained in the results
modColors <- base::rownames(x = mp$preservation$observed[[ref]][[test]])
moduleSizes <- mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods <- !(modColors %in% c("grey", "gold"));
# Text labels for points
text <- modColors[plotMods];
# Auxiliary convenience variable
plotData <- base::cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains <- base::c("Preservation Median rank", "Preservation Zsummary");

max(moduleSizes[plotMods])
min(moduleSizes[plotMods])

w <- 10
h <- 5
# Start the plot
sizeGrWindow(w, h);
grDevices::pdf(file = base::file.path(fig.wgc.dir, "WGCNA_mouse_module_preservation_in_WholeHumanGenes_Zsummary_medianRank.pdf"),
               width = w, height = h, onefile = T, paper = "special")
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2) {
  curr.min <- base::min(plotData[, p], na.rm = TRUE)
  curr.max <- base::max(plotData[, p], na.rm = TRUE)
  # Adjust ploting ranges appropriately
  if (2 == p) {
    if (curr.min > (-curr.max/10)) {
      curr.min <- (-curr.max/10)
    }
    ylim <- base::c(curr.min - 0.1 * (curr.max - curr.min), curr.max + 0.1 * (curr.max - curr.min))
  } else {
    ylim <- base::c(curr.max + 0.1 * (curr.max - curr.min), curr.min - 0.1 * (curr.max - curr.min))
  }
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = base::c(20, 300), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08)
  # For Zsummary, add threshold lines
  if (2 == p) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
# If plotting into a file, close it
grDevices::dev.off()


### Human DKD snRNA PT ------

sobj <- base::readRDS(file = base::file.path(rds.dir, "after_neighbor.RDS"))
curr.sobj <- sobj

all.conditions <- base::unique(x = base::as.character(x = curr.sobj$condition))
all.celltypes <- base::unique(x = base::as.character(x = curr.sobj$cell.type))
all.conditions
all.celltypes

curr.sobj$condition.cell.type <- base::paste(
  base::as.character(x = curr.sobj$condition), 
  base::as.character(x = curr.sobj$cell.type), 
  sep = "_"
)
table(curr.sobj$condition.cell.type)
sum(table(curr.sobj$condition.cell.type))

nfeatures <- 4000
curr.sobj <- Seurat::FindVariableFeatures(object = curr.sobj, assay = "RNA", selection.method = "vst",
                                          nfeatures = nfeatures, verbose = F)
curr.sobj
head(curr.sobj[[]], n=3L)

HVGs <- Seurat::VariableFeatures(object = curr.sobj, assay = "RNA")
length(HVGs)
head(HVGs)

k <- 10
min.cells <- 20
split.sobj <- Seurat::SplitObject(object = curr.sobj, split.by = "condition.cell.type")
select.sobj <- base::list()
for (curr.condition in all.conditions) {
  for (curr.celltype in all.celltypes) {
    curr.group <- base::paste(curr.condition, curr.celltype, sep = "_")
    if (!utils::hasName(x = split.sobj, name = curr.group)) {
      base::cat(base::sprintf(fmt = "\n%s does not exist.\n", curr.group))
      next
    }
    curr.sobj <- split.sobj[[curr.group]]
    if (base::ncol(x = curr.sobj) < min.cells) {
      base::cat(base::sprintf(fmt = "\n%s has %s cells.\nSkip %s.\n", curr.group, base::ncol(x = curr.sobj), curr.group))
      next
    }
    base::cat(base::sprintf(fmt = "\nProcessing %s ...\n", curr.group))
    curr.sobj <- curr.sobj[HVGs,]
    curr.metacell.sobj <- hdWGCNA::construct_metacells(seurat_obj = curr.sobj, name = curr.group, k = k, 
                                                       reduction = "umap", assay = "RNA", slot = "data")
    curr.metacell.sobj$condition <- curr.condition
    curr.metacell.sobj$cell.type <- curr.celltype
    curr.metacell.sobj$condition.cell.type <- curr.group
    select.sobj[[curr.group]] <- curr.metacell.sobj
  }
}
metacell.sobj <- merge(x = select.sobj[[1]], y = select.sobj[2:base::length(x = select.sobj)], 
                       add.cell.ids = NULL, merge.data = T)
metacell.sobj
head(x = metacell.sobj[[]], n = 3L)
tail(x = metacell.sobj[[]], n = 3L)
table(metacell.sobj$condition)
table(metacell.sobj$cell.type)
table(metacell.sobj$condition.cell.type)
base::saveRDS(
  object = metacell.sobj, 
  file = base::file.path(rds.dir, base::sprintf(fmt = "metacell_condition_celltype_%sHVGs_k%s.RDS", length(HVGs), k))
)
split.sobj <- NULL
select.sobj <- NULL
curr.sobj <- NULL
base::gc()

orig.metacell.sobj <- metacell.sobj

metacell.sobj <- orig.metacell.sobj
metacell.sobj <- Seurat::ScaleData(object = metacell.sobj, assay = "RNA", 
                                   features = base::row.names(x = metacell.sobj), verbose = F)
npcs <- 51
metacell.sobj <- Seurat::RunPCA(object = metacell.sobj, assay = "RNA", 
                                features = base::row.names(x = metacell.sobj), npcs = npcs, verbose = F)
metadata <- metacell.sobj@meta.data
pca.embed <- Seurat::Embeddings(object = metacell.sobj, reduction = "pca")
reticulate::use_python(python = "/home/jzhou88/installed/miniconda3/envs/r4/bin/python")
harmonypy <- reticulate::import(module = "harmonypy")
ho <- harmonypy$run_harmony(pca.embed, metadata, c("condition"), max_iter_harmony = 100L, verbose = T)
the <- ho$Z_corr
ho <- NULL
base::gc()
harmony.embed <- base::t(x = the)
base::rownames(x = harmony.embed) <- base::rownames(x = pca.embed)
base::colnames(x = harmony.embed) <- base::paste0("harmony_", 1:npcs)
metacell.sobj[["harmony"]] <- Seurat::CreateDimReducObject(embeddings = harmony.embed, assay = "RNA", 
                                                           key = "harmony_")
metacell.sobj <- Seurat::RunUMAP(object = metacell.sobj, dims = 1:npcs, reduction = "harmony", 
                                 umap.method = "uwot", return.model = T, metric = "cosine", 
                                 verbose = F, reduction.name = "umap", reduction.key = "UMAP_")

metacell.sobj$condition <- base::factor(x = metacell.sobj$condition, levels = base::c("Healthy", "DKD"))
Seurat::Idents(object = metacell.sobj) <- "condition"
w <- 6
h <- 5
dp <- Seurat::DimPlot(object = metacell.sobj, reduction = "umap", label = F) + 
  Seurat::NoAxes()
grDevices::pdf(file = base::file.path(fig.wgc.dir,
                                      base::sprintf(fmt = "umap_metacell_condition_celltype_%sHVGs_k%s.pdf", length(HVGs), k)), 
               width = w, height = h, onefile = T, paper = "special")
base::print(x = dp)
grDevices::dev.off()

metacell.sobj$cell.type <- base::factor(x = metacell.sobj$cell.type, levels = base::c("PT", "injured PT"))
Seurat::Idents(object = metacell.sobj) <- "cell.type"
w <- 6
h <- 5
dp <- Seurat::DimPlot(object = metacell.sobj, reduction = "umap", label = F) + 
  Seurat::NoAxes()
grDevices::pdf(
  file = base::file.path(fig.wgc.dir,
                         base::sprintf(fmt = "umap_metacell_condition_celltype_%sHVGs_k%s_GroupedByCellTypes.pdf", length(HVGs), k)),
  width = w, height = h, onefile = T, paper = "special"
)
base::print(x = dp)
grDevices::dev.off()

metacell.sobj <- orig.metacell.sobj
base::options(stringsAsFactors = FALSE)
WGCNA::enableWGCNAThreads()
datExpr <- base::as.data.frame(x = Seurat::GetAssayData(object = metacell.sobj, slot = "data", assay = "RNA"))
datExpr <- base::as.data.frame(x = base::t(x = datExpr))
datExpr <- datExpr[, WGCNA::goodGenes(datExpr = datExpr)]
base::dim(x = datExpr)

powers <- base::c(base::seq(from = 1, to = 10, by = 1), base::seq(from = 12, to = 30, by = 2))
powerTable <- base::list(
  data = WGCNA::pickSoftThreshold(
    data = datExpr, 
    powerVector = powers, 
    blockSize = 20000,
    corFnc = "bicor", 
    networkType = "signed", 
    verbose = 2
  )[[2]]
)

colors <- base::c("red")
plotCols <- base::c(2, 5, 6, 7)
colNames <- base::c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity", "Max connectivity")
ylim <- base::matrix(data = NA, nrow = 2, ncol = 4)
for (col in 1:base::length(x = plotCols)) {
  ylim[1, col] <- base::min(ylim[1, col], powerTable$data[, plotCols[col]], na.rm = TRUE)
  ylim[2, col] <- base::max(ylim[2, col], powerTable$data[, plotCols[col]], na.rm = TRUE)
}
par(mfcol = base::c(2,2))
par(mar = base::c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7
for (col in 1:base::length(x = plotCols)){
  plot(powerTable$data[,1], -sign(powerTable$data[,3])*powerTable$data[,2], xlab="Soft Threshold (power)",
       ylab=colNames[col],type="n", ylim = ylim[, col], main = colNames[col])
  addGrid()
  
  if (col==1) {
    text(powerTable$data[,1], -sign(powerTable$data[,3])*powerTable$data[,2], labels=powers,cex=cex1,col=colors[1])
  } else {
    text(powerTable$data[,1], powerTable$data[,plotCols[col]], labels=powers,cex=cex1,col=colors[1])
  }
}

multiExpr <- list()
multiExpr[["PT"]] <- list(data=datExpr)
checkSets(multiExpr)
head(multiExpr$PT$data)

cor.type <- "pearson"
min.module.size <- 50

softPower <- 6
net <- WGCNA::blockwiseConsensusModules(multiExpr, 
                                        blocks = NULL, 
                                        maxBlockSize = 20000, ## This should be set to a smaller size if the user has limited RAM
                                        corType = cor.type, ## no use for bicor
                                        power = softPower,
                                        consensusQuantile = 0.3,
                                        networkType = "signed",
                                        TOMType = "unsigned",
                                        TOMDenom = "min",
                                        scaleTOMs = TRUE, scaleQuantile = 0.8,
                                        sampleForScaling = TRUE, sampleForScalingFactor = 1000,
                                        useDiskCache = T, chunkSize = NULL,
                                        deepSplit = 4,
                                        pamStage = F,
                                        detectCutHeight = 0.995, minModuleSize = min.module.size,
                                        mergeCutHeight = 0.2,
                                        saveConsensusTOMs = F,
                                        consensusTOMFilePattern = "ConsensusTOM-block.%b.rda", 
                                        numericLabels = F, 
                                        nThreads = 8)
base::saveRDS(
  object = net, 
  file = base::file.path(
    rds.dir, 
    base::sprintf(fmt = "WGCNA_consensus_modules_%s_MinModSize%s.RDS", cor.type, min.module.size)
  )
)
table(net$colors)

consTree <- net$dendrograms[[1]]
moduleLabels <- net$colors;
moduleColors <- base::as.character(x = moduleLabels)
w <- 8
h <- 5
grDevices::pdf(
  file = base::file.path(
    fig.wgc.dir, 
    base::sprintf(fmt = "WGCNA_consensus_gene_dendrogram_%s_MinModSize%s.pdf", cor.type, min.module.size)
  ), width = w, height = h, onefile = T, paper = "special"
)
plotDendroAndColors(dendro = consTree, 
                    colors = moduleColors, 
                    groupLabels = "Module", 
                    autoColorHeight = T, 
                    dendroLabels = F, 
                    hang = 0.03, 
                    addGuide = T, 
                    guideHang = 0.05, 
                    main = "")
grDevices::dev.off()

MEs <- WGCNA::moduleEigengenes(expr = multiExpr[[1]]$data, colors = moduleColors, nPC = 1)$eigengenes
MEs <- WGCNA::orderMEs(MEs = MEs)
KMEs <- WGCNA::signedKME(datExpr = datExpr, datME = MEs, outputColumnName = "kME", corFnc = "bicor")
geneInfo <- base::as.data.frame(x = base::cbind(base::colnames(datExpr), moduleLabels, moduleColors, KMEs))
base::colnames(geneInfo)[1] <- "GeneSymbol"
base::colnames(geneInfo)[2] <- "ModuleLabel"
base::colnames(geneInfo)[3] <- "Initially.Assigned.Module.Color"
utils::write.csv(
  x = geneInfo,
  file = base::file.path(
    rst.wgc.dir, 
    base::sprintf(fmt = "WGCNA_geneInfoSigned_%s_MinModSize%s.csv", cor.type, min.module.size)
  )
)

module.labels <- base::unique(x = geneInfo$ModuleLabel)
modules <- base::list()
for (curr.label in module.labels) {
  base::cat(base::sprintf(fmt = "\nProcessing Module %s ...\n", curr.label))
  modules[[curr.label]] <- geneInfo[geneInfo$ModuleLabel == curr.label,]
}
openxlsx::write.xlsx(
  x = modules, 
  file = base::file.path(
    rst.wgc.dir, base::sprintf(fmt = "WGCNA_SignedGeneInfo_PerModule_%s_MinModSize%s.xlsx", cor.type, min.module.size)
  ), startCol = 1, startRow = 1, colNames = T, rowNames = F, overwrite = T
)

human.datExpr <- datExpr
human.colors <- net$colors

## Converting mouse genes into human genes does not serve the purpose, 
## because different mouse genes can be mapped into the same human gene.
h2m.genes <- nichenetr::convert_human_to_mouse_symbols(
  symbols = base::colnames(x = human.datExpr), 
  version = 1
)

base::print(x = base::length(x = h2m.genes))
base::print(x = base::sum(base::is.na(x = h2m.genes)))
base::print(x = utils::head(x = h2m.genes))

h2m.genes[base::is.na(x = h2m.genes)] <- base::colnames(x = human.datExpr)[base::is.na(x = h2m.genes)]

base::print(x = base::length(x = h2m.genes))
base::print(x = base::sum(base::is.na(x = h2m.genes)))
base::print(x = utils::head(x = h2m.genes))

base::colnames(x = human.datExpr) <- h2m.genes
head(colnames(human.datExpr))

## Converting mouse genes into human genes does not serve the purpose, 
## because different mouse genes can be mapped into the same human gene.
h2m.genes <- nichenetr::convert_human_to_mouse_symbols(
  symbols = base::names(x = human.colors), 
  version = 1
)

base::print(x = base::length(x = h2m.genes))
base::print(x = base::sum(base::is.na(x = h2m.genes)))
base::print(x = utils::head(x = h2m.genes))

h2m.genes[base::is.na(x = h2m.genes)] <- base::names(x = human.colors)[base::is.na(x = h2m.genes)]

base::print(x = base::length(x = h2m.genes))
base::print(x = base::sum(base::is.na(x = h2m.genes)))
base::print(x = utils::head(x = h2m.genes))

base::names(x = human.colors) <- h2m.genes
head(names(human.colors))

mouse.sobj <- base::readRDS(
  file = "/home/jzhou88/projects/KidneyDisease/2022-04-21/MouseDkdRnaNucleusPtSeuratHarmonypy/RDS/after_neighbor.RDS"
)

mouse.datExpr <- base::as.data.frame(x = Seurat::GetAssayData(object = mouse.sobj, slot = "data", assay = "RNA"))
mouse.datExpr <- base::as.data.frame(x = base::t(x = mouse.datExpr))
mouse.datExpr <- mouse.datExpr[, WGCNA::goodGenes(datExpr = mouse.datExpr)]
base::dim(x = mouse.datExpr)

setLabels <- c("Human", "Mouse")
multiExpr <- base::list(Human = base::list(data = human.datExpr), Mouse = base::list(data = mouse.datExpr))
multiColor <- list(Human = human.colors)

mp <- WGCNA::modulePreservation(
  multiData = multiExpr,
  multiColor = multiColor,
  referenceNetworks = 1,
  nPermutations = 200,
  randomSeed = 1,
  quickCor = 0,
  verbose = 3
)
base::saveRDS(object = mp, file = base::file.path(rds.dir, "WGCNA_module_preservation_WholeMouseGenes.RDS"))

ref <- 1
test <- 2
statsObs <- base::cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ <-base::cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])

base::print(x = base::cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
                            signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

# Module labels and module sizes are also contained in the results
modColors <- base::rownames(x = mp$preservation$observed[[ref]][[test]])
moduleSizes <- mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods <- !(modColors %in% c("grey", "gold"));
# Text labels for points
text <- modColors[plotMods];
# Auxiliary convenience variable
plotData <- base::cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains <- base::c("Preservation Median rank", "Preservation Zsummary")

max(moduleSizes[plotMods])
min(moduleSizes[plotMods])

w <- 10
h <- 5
# Start the plot
sizeGrWindow(w, h);
grDevices::pdf(file = base::file.path(fig.wgc.dir, "WGCNA_human_module_preservation_in_WholeMouseGenes_Zsummary_medianRank.pdf"),
               width = w, height = h, onefile = T, paper = "special")
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2) {
  curr.min <- base::min(plotData[, p], na.rm = TRUE)
  curr.max <- base::max(plotData[, p], na.rm = TRUE)
  # Adjust ploting ranges appropriately
  if (2 == p) {
    if (curr.min > (-curr.max/10)) {
      curr.min <- (-curr.max/10)
    }
    ylim <- base::c(curr.min - 0.1 * (curr.max - curr.min), curr.max + 0.1 * (curr.max - curr.min))
  } else {
    ylim <- base::c(curr.max + 0.1 * (curr.max - curr.min), curr.min - 0.1 * (curr.max - curr.min))
  }
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = base::c(40, 520), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08)
  # For Zsummary, add threshold lines
  if (2 == p) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
# If plotting into a file, close it
grDevices::dev.off()
