# seurat_analysis.R
base::cat("\nRunning \"seurat_analysis.R\" ...\n")

sobj <- base::readRDS(file = "/home/jzhou88/projects/KidneyDisease/2020-11-04/ScanpyHarmony14/RDS/after_QC.RDS")
split.seurat.obj <- T
stop.after.hvf <- F
npcs <- 121 # 438,686 cells in total

### Standard Seurat workflow
if (split.seurat.obj) {
  split.sobj <- Seurat::SplitObject(object = sobj, split.by = "sample")
  sobj <- NULL
  base::gc()
  hvf <- NULL
  for (curr.mouse in curr.mice) {
    if (!utils::hasName(x = split.sobj, name = curr.mouse)) {
      base::cat(base::sprintf(fmt = "\nDo not have %s.\n", curr.mouse))
      next
    }
    ## Normalize the data
    split.sobj[[curr.mouse]] <- Seurat::NormalizeData(object = split.sobj[[curr.mouse]], assay = "RNA", 
                                                      normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
    ## Identify highly variable genes (feature selection)
    split.sobj[[curr.mouse]] <- Seurat::FindVariableFeatures(object = split.sobj[[curr.mouse]], assay = "RNA", 
                                                             selection.method = "vst", nfeatures = 2000, verbose = F)
    hvf <- base::c(hvf, Seurat::VariableFeatures(object = split.sobj[[curr.mouse]], assay = "RNA"))
  }
  hvf <- base::unique(x = hvf)
  base::saveRDS(object = hvf, file = base::file.path(rds.dir, "hvf.RDS"))
  if (stop.after.hvf) {
    base::stop("stop after hvf.")
  }
  sobj <- merge(x = split.sobj[[1]], y = split.sobj[2:base::length(x = split.sobj)], 
                add.cell.ids = NULL, merge.data = T, project = proj.name)
  base::cat("\n")
  base::print(x = utils::head(x = sobj[[]], n = 3L))
  base::print(x = utils::tail(x = sobj[[]], n = 3L))
  split.sobj <- NULL
  base::gc()
  Seurat::VariableFeatures(object = sobj, assay = "RNA") <- hvf
  features <- base::rownames(x = sobj)
  vst.variable <- base::rep(x = FALSE, times = base::length(x = features))
  vst.variable[features %in% hvf] <- TRUE
  base::names(x = vst.variable) <- features
  sobj[["RNA"]] <- Seurat::AddMetaData(object = sobj[["RNA"]], metadata = vst.variable, col.name = "vst.variable")
  base::cat("\n")
  base::print(x = utils::head(x = sobj[["RNA"]][[]]))
  
} else {
  ## Normalize the data
  sobj <- Seurat::NormalizeData(object = sobj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000, verbose = F)
  ## Identify highly variable genes (feature selection)
  sobj <- Seurat::FindVariableFeatures(object = sobj, assay = "RNA", selection.method = "vst", nfeatures = 2000, verbose = F)
}

## Scale the data
sobj <- Seurat::ScaleData(object = sobj, assay = "RNA", features = NULL, verbose = F)

## Perform linear dimensional reduction
## ~3k cells use 50 PCs and ~250k cells use 100 PCs.
## 20k~40k cells use 61, 50k~70k cells use 71, 80k~100k cells use 81, 110k~190k cells use 91
## 200k~290k cells use 101, 300k~390k cells use 111, 400k~490k cells use 121
## 500k~590k cells use 131, 600k~690k cells use 141, 700k~790k cells use 151
base::cat(base::sprintf(fmt = "\nCalculating %s PCs ...\n", npcs))
sobj <- Seurat::RunPCA(object = sobj, assay = "RNA", features = NULL, npcs = npcs, verbose = F)
base::saveRDS(object = sobj, file = base::file.path(rds.dir, "after_pca.RDS"))
PlotIntegration(seurat.obj = sobj, act.assay = "RNA", features.to.plot = "PC_1", save.dir = fig.pc.dir, fn.suf = "_PC1")

## Integrate the data (or rather, correct batch effect) using harmonypy
base::cat("\n")
base::print(x = sobj@reductions)
metadata <- sobj@meta.data
base::cat("\n")
base::print(x = base::dim(x = metadata))
pca.embed <- Seurat::Embeddings(object = sobj, reduction = "pca")
base::print(x = base::dim(x = pca.embed))
harmonypy <- reticulate::import(module = "harmonypy")
ho <- harmonypy$run_harmony(pca.embed, metadata, c("sample"), max_iter_harmony = 100L, verbose = T)
the <- ho$Z_corr
ho <- NULL
base::gc()
base::print(x = base::dim(x = the))
base::cat("\n")
base::print(x = utils::head(x = the[,1:6]))
harmony.embed <- base::t(x = the)
base::cat("\n")
base::print(x = base::dim(x = harmony.embed))
base::rownames(x = harmony.embed) <- base::rownames(x = pca.embed)
base::colnames(x = harmony.embed) <- base::paste0("harmony_", 1:npcs)
base::cat("\n")
base::print(x = utils::head(x = harmony.embed[,1:6]))
sobj[["harmony"]] <- Seurat::CreateDimReducObject(embeddings = harmony.embed, assay = "RNA", key = "harmony_")
base::cat("\n")
base::print(x = sobj@reductions)
base::saveRDS(object = sobj, file = base::file.path(rds.dir, "after_harmony.RDS"))
PlotIntegration(seurat.obj = sobj, act.assay = "RNA", features.to.plot = "harmony_1", save.dir = fig.pc.dir, fn.suf = "_harmony1")
PlotPrincipalComponents(seurat.obj = sobj, save.dir = fig.pc.dir, dim.to.plot = npcs, red.to.use = "harmony", do.elbow = F)

pcs <- 121
red <- "harmony"
## Run non-linear dimensional reduction
sobj <- Seurat::RunUMAP(object = sobj, dims = 1:pcs, reduction = red, umap.method = "uwot", return.model = T, metric = "cosine", 
                        verbose = F, reduction.name = "umap", reduction.key = "UMAP_")
## Cluster the cells
sobj <- Seurat::FindNeighbors(object = sobj, reduction = red, dims = 1:pcs, nn.eps = 0, verbose = F)
base::saveRDS(object = sobj, file = base::file.path(rds.dir, "after_neighbor.RDS"))

act.assay <- "RNA"
alg.to.use <- 4
red.to.plot <- "umap"
markers.to.plot <- default.known.markers

Seurat::DefaultAssay(object = sobj) <- act.assay
base::cat("\n")
base::print(x = sobj)
### Cluster the cells
## ~3k cells use resolution 0.4~1.4 and ~250k cells use resolution 3.
## 1k~9k cells use resolutions up to 1.5, 10k~90k cells use resolutions up to 2, 
## 100k~190k cells use resolutions up to 2.5, 200k~290k cells use resolutions up to 3
sobj <- Seurat::FindClusters(object = sobj, resolution = input.res, method = "igraph", algorithm = alg.to.use, n.start = 100, n.iter = 10, 
                             random.seed = 1, verbose = F)
base::saveRDS(object = sobj, file = base::file.path(rds.dir, base::sprintf(fmt = "after_clustering_res%s.RDS", input.res)))

PlotClusters(seurat.obj = sobj, assay.clustered = act.assay, resols = base::c(input.res), save.dir = fig.clu.dir, 
             red.to.plot = red.to.plot, do.red = T, split.clusters = T, do.split = T, var.to.split = "sample", 
             do.bar = T, var.to.bar = "sample", do.dot = T, assay.to.plot = "RNA", markers.to.plot = markers.to.plot, 
             fig.suf = "")
