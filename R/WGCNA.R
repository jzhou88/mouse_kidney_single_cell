# WGCNA.R
base::cat("\nRunning \"WGCNA.R\" ...\n")

sobj <- base::readRDS(file = base::file.path(rds.dir, "after_neighbor.RDS"))

models <- sobj$sample
cell.names <- base::names(x = models)
models <- base::as.character(x = models)
disease.models <- base::names(x = select.rna.disease.samples)
for (dm in disease.models) {
  models[models %in% select.rna.disease.samples[[dm]]] <- dm
}
models[models %in% select.rna.control.samples] <- "Control"
base::names(x = models) <- cell.names
sobj <- Seurat::AddMetaData(object = sobj, metadata = models, col.name = "model")
Seurat::Idents(object = sobj) <- "cell.type"

HVGs <- Seurat::VariableFeatures(object = sobj, assay = "RNA")
k <- 10
split.sobj <- Seurat::SplitObject(object = sobj, split.by = "model")
for (curr.model in base::names(x = split.sobj)) {
  base::cat(base::sprintf(fmt = "\nProcessing %s ...\n", curr.model))
  curr.sobj <- split.sobj[[curr.model]]
  
  curr.sobj <- curr.sobj[HVGs,]
  
  curr.metacell.sobj <- hdWGCNA::construct_metacells(seurat_obj = curr.sobj, name = curr.model, k = k, reduction = "umap", assay = "RNA", slot = "data")
  curr.metacell.sobj$model <- curr.model
  split.sobj[[curr.model]] <- curr.metacell.sobj
}
metacell.sobj <- merge(x = split.sobj[[1]], y = split.sobj[2:base::length(x = split.sobj)], 
                       add.cell.ids = NULL, merge.data = T)

orig.metacell.sobj <- metacell.sobj
metacell.sobj <- Seurat::ScaleData(object = metacell.sobj, assay = "RNA", features = base::row.names(x = metacell.sobj), verbose = F)
npcs <- 51
metacell.sobj <- Seurat::RunPCA(object = metacell.sobj, assay = "RNA", features = base::row.names(x = metacell.sobj), npcs = npcs, verbose = F)
metadata <- metacell.sobj@meta.data
pca.embed <- Seurat::Embeddings(object = metacell.sobj, reduction = "pca")
harmonypy <- reticulate::import(module = "harmonypy")
ho <- harmonypy$run_harmony(pca.embed, metadata, c("model"), max_iter_harmony = 100L, verbose = T)
the <- ho$Z_corr
ho <- NULL
base::gc()
harmony.embed <- base::t(x = the)
base::rownames(x = harmony.embed) <- base::rownames(x = pca.embed)
base::colnames(x = harmony.embed) <- base::paste0("harmony_", 1:npcs)
metacell.sobj[["harmony"]] <- Seurat::CreateDimReducObject(embeddings = harmony.embed, assay = "RNA", key = "harmony_")
metacell.sobj <- Seurat::RunUMAP(object = metacell.sobj, dims = 1:npcs, reduction = "harmony", umap.method = "uwot", return.model = T, metric = "cosine", 
                                 verbose = F, reduction.name = "umap", reduction.key = "UMAP_")

Seurat::Idents(object = metacell.sobj) <- "model"
known.models <- base::c("Control", base::names(x = select.rna.disease.samples))
Seurat::Idents(object = metacell.sobj) <- base::factor(x = Seurat::Idents(object = metacell.sobj), levels = known.models)

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

softPower <- 9
net <- WGCNA::blockwiseConsensusModules(multiExpr, 
                                        blocks = NULL, 
                                        maxBlockSize = 20000, ## This should be set to a smaller size if the user has limited RAM
                                        corType = "bicor",
                                        power = softPower,
                                        networkType = "signed",
                                        TOMType = "unsigned",
                                        TOMDenom = "min",
                                        useDiskCache = T, chunkSize = NULL,
                                        deepSplit = 4,
                                        pamStage = F,
                                        detectCutHeight = 0.995, minModuleSize = 50,
                                        mergeCutHeight = 0.2,
                                        saveConsensusTOMs = F,
                                        consensusTOMFilePattern = "ConsensusTOM-block.%b.rda", 
                                        numericLabels = F, 
                                        nThreads = 16)

consTree <- net$dendrograms[[1]]
moduleLabels <- net$colors;
moduleColors <- base::as.character(x = moduleLabels)
plotDendroAndColors(dendro = consTree, 
                    colors = moduleColors, 
                    groupLabels = "Module", 
                    autoColorHeight = T, 
                    dendroLabels = F, 
                    hang = 0.03, 
                    addGuide = T, 
                    guideHang = 0.05, 
                    main = "")

MEs <- WGCNA::moduleEigengenes(expr = multiExpr[[1]]$data, colors = moduleColors, nPC = 1)$eigengenes
MEs <- WGCNA::orderMEs(MEs = MEs)
KMEs <- WGCNA::signedKME(datExpr = datExpr, datME = MEs, outputColumnName = "kME", corFnc = "bicor")
geneInfo <- base::as.data.frame(x = base::cbind(base::colnames(datExpr), moduleLabels, moduleColors, KMEs))
base::colnames(geneInfo)[1] <- "GeneSymbol"
base::colnames(geneInfo)[2] <- "ModuleLabel"
base::colnames(geneInfo)[3] <- "Initially.Assigned.Module.Color"
utils::write.csv(x = geneInfo, file = base::file.path(rst.wgc.dir, "WGCNA_geneInfoSigned.csv"))

n_genes <- 25
module_labels <- unique(geneInfo$ModuleLabel)
modules <- unique(geneInfo$Initially.Assigned.Module.Color)
module_labels <- module_labels[order(modules)]
modules <- modules[order(modules)]
mod_colors <- modules
names(mod_colors) <- module_labels
module_list <- lapply(modules, function(mod){
  cur <- subset(geneInfo, Initially.Assigned.Module.Color == mod)
  cur[,c('GeneSymbol', paste0('kME', mod))] %>%
    top_n(n_genes) %>% .$GeneSymbol
})
names(module_list) <- module_labels

test.sobj <- sobj
test.sobj <- Seurat::AddModuleScore(object = test.sobj, features = module_list, pool = base::row.names(x = test.sobj), nbin = 24, k = F, assay = "RNA", name = "PT_module")

# get module scores from scRNA data
cur_celltype <- "PT"
module_df <- test.sobj@meta.data
features <- names(module_df)[grepl(cur_celltype, names(module_df))]
feat2mod <- module_labels
names(feat2mod) <- features
module_df <- module_df[,c('model', features)]

# compute average module score for each cluster:
tmp <- lapply(unique(module_df$model), function(i){
  cur_df <- module_df %>% subset(model == i)
  data.frame(
    value=as.numeric(colSums(cur_df[,features]) / nrow(cur_df)),
    cluster = i,
    feature = features
  )
})
plot_df <- Reduce(rbind, tmp)
plot_df$feature <- feat2mod[plot_df$feature]

# remove gray module:
plot_df <- subset(plot_df, feature != 'grey')
mod_colors <- mod_colors[module_labels != 'grey']
module_labels <- module_labels[module_labels != 'grey']
names(mod_colors) <- module_labels
# reshape to wide format
plot_df <- reshape2::dcast(plot_df, feature ~ cluster)
rownames(plot_df) <- plot_df$feature
plot_df <- plot_df %>% dplyr::select(-c(feature))

module.order <- base::c("blue", "turquoise", "yellow", "brown")
model.order <- base::c("FAN", "UUO", "Notch", "PGC1a", "longIRI14d", "shortIRI14d", "Control")

zScore <- function(x){(x - mean(x)) /sd(x)}
matrix_z <- apply(plot_df, 1, zScore) %>% t()
matrix_z <- matrix_z[,model.order]
# keep full values to plot onto the heatmap
matrix_full <- matrix_z
matrix_anno <- ifelse(matrix_full >= quantile(matrix_full, 0.80), signif(matrix_full,2), '')
# change the range
range(matrix_z)
matrix_z <- ifelse(matrix_z >= 2, 2, matrix_z)
matrix_z <- ifelse(matrix_z <= -2, -2, matrix_z)
# rename modules
rownames(matrix_z) <- module_labels
rownames(matrix_anno) <- module_labels
matrix_anno <- matrix_anno[module.order,]
matrix_z <- matrix_z[module.order,]
# color the rows by scWGCNA module color, cols by single-cell clusters
cluster_colors <- scales::hue_pal()(base::length(x = known.models))
base::names(x = cluster_colors) <- known.models
cluster_colors <- cluster_colors[colnames(matrix_z)]

column_ha <- HeatmapAnnotation(
  Cluster = colnames(matrix_z),
  col = list(
    Cluster = cluster_colors
  )
)

row_ha <- rowAnnotation(
  module = rownames(matrix_z),
  col = list(
    module = mod_colors
  )
)

ComplexHeatmap::Heatmap(
  matrix_z,
  cluster_rows = F,
  cluster_columns = F,
  top_annotation = column_ha,
  left_annotation = row_ha,
  width = unit(3.5, "in"), 
  height = unit(2, "in"),
  show_heatmap_legend = T,
  use_raster = T,
  cell_fun = function(j,i,x,y,w,h,col){
    grid.text(matrix_anno[i,j], x, y)
  }
)

mm.orgdb <- "org.Mm.eg.db"
go.keytype <- "ENTREZID"
gene.ont <- "BP"
mm.go.file <- base::file.path(rds.dir, base::sprintf(fmt = "%s_%s_%s_GO.RDS", mm.orgdb, go.keytype, gene.ont))
mm.go <- NULL
if (base::file.exists(mm.go.file)) {
  mm.go <- base::readRDS(file = mm.go.file)
  
} else {
  mm.go <- GOSemSim::godata(OrgDb = mm.orgdb, keytype = go.keytype, ont = gene.ont)
  base::saveRDS(object = mm.go, file = mm.go.file)
}
dist.method <- "Rel"
sim.thresh <- 0.7
fn <- "WGCNA_Module_GOTERM_BP_FAT"
go.terms <- Xlsx2List(xlsx.file = base::file.path(rst.wgc.dir, base::sprintf(fmt = "david2021/%s.xlsx", fn)), rowNames = F, cols = base::c(2, 13))

modules <- base::c("blue", "turquoise", "yellow", "brown")
max.p.val <- 0.001
for (curr.module in modules) {
  curr.go <- go.terms[[curr.module]]
  curr.go$ID <- stringr::str_split_fixed(string = curr.go$Term, pattern = "~", n = 2)[,1]
  base::row.names(x = curr.go) <- curr.go$ID
  curr.go <- curr.go[base::order(curr.go$FDR, decreasing = F),]
  curr.go <- curr.go[base::which(curr.go$FDR < max.p.val),]
  curr.go$Term <- NULL
  simMatrix <- rrvgo::calculateSimMatrix(
    x = curr.go$ID, 
    orgdb = mm.orgdb, 
    keytype = go.keytype, 
    semdata = mm.go, 
    ont = gene.ont, 
    method = dist.method
  )
  scores <- stats::setNames(object = (-1) * base::log10(x = curr.go$FDR), nm = curr.go$ID)
  reducedTerms <- rrvgo::reduceSimMatrix(
    simMatrix = simMatrix,
    scores = scores,
    threshold = sim.thresh,
    orgdb = mm.orgdb,
    keytype = go.keytype
  )
  treemap::treemap(dtf=reducedTerms, index=base::c("parentTerm", "term"), vSize="score", type="index", title="", 
                   palette=scales::hue_pal(direction = 1)(base::length(x=base::unique(x=reducedTerms$parent))), 
                   fontsize.labels=11, 
                   fontcolor.labels=base::c("white", "black"), 
                   fontfamily.labels="sans", 
                   border.col=NA,
                   border.lwds=0,
                   inflate.labels=T,
                   bg.labels=0,
                   aspRatio = w/h
  )
}
