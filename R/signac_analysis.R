# signac_analysis.R
base::cat("\nRunning \"signac_analysis.R\" ...\n")


sobj <- base::readRDS(file = base::file.path(curr.out.dir, "../HqAdult1/RDS/after_SVD.RDS"))
ssv <- 2
nsv <- 60

base::cat(base::sprintf(fmt = "\nProcessing %s ...\n", picked.name))
base::cat("\n")
base::print(x = sobj)
orig.sobj <- sobj

## Perform non-linear dimension reduction.
sobj <- Seurat::RunUMAP(object = sobj, reduction = "lsi", dims = ssv:nsv, verbose = F)

red.to.plot <- "umap"
var.to.group <- "sample"
var.to.split <- "sample"
f <- 18
h <- 7
dp <- Seurat::DimPlot(object = sobj, reduction = red.to.plot, group.by = var.to.group) + 
  Seurat::FontSize(x.text = f, y.text = f, x.title = f, y.title = f)
w1 <- 9
grDevices::pdf(file = base::file.path(fig.int.dir, base::sprintf("%s_BI_G%s.pdf", red.to.plot, var.to.group)), 
               width = w1, height = h, onefile = T, paper = "special")
base::print(x = dp)
grDevices::dev.off()

sdp <- Seurat::DimPlot(object = sobj, reduction = red.to.plot, group.by = var.to.group, split.by = var.to.split) + 
  Seurat::FontSize(x.text = f, y.text = f, x.title = f, y.title = f) + 
  Seurat::NoLegend()
w2 <- 7 * base::length(x = base::unique(x = sobj@meta.data[[var.to.split]]))
grDevices::pdf(file = base::file.path(fig.int.dir, base::sprintf("%s_BI_G%s_S%s.pdf", red.to.plot, var.to.group, var.to.split)), 
               width = w2, height = h, onefile = T, paper = "special")
base::print(x = sdp)
grDevices::dev.off()

## Perform integration using Harmonypy.
sobj <- orig.sobj
base::cat("\n")
base::print(x = sobj@reductions)
PlotIntegration(seurat.obj = sobj, act.assay = "ATAC", features.to.plot = "LSI_1", save.dir = fig.int.dir, fn.suf = "_LSI1")

metadata <- sobj@meta.data
base::cat("\n")
base::print(x = base::dim(x = metadata))
lsi.embed <- Seurat::Embeddings(object = sobj, reduction = "lsi")
base::print(x = base::dim(x = lsi.embed))
harmonypy <- reticulate::import(module = "harmonypy")
ho <- harmonypy$run_harmony(lsi.embed, metadata, c("sample"), max_iter_harmony = 100L, verbose = T)
the <- ho$Z_corr
ho <- NULL
base::gc()
base::print(x = base::dim(x = the))
base::cat("\n")
base::print(x = utils::head(x = the[,1:6]))
harmony.embed <- base::t(x = the)
base::cat("\n")
base::print(x = base::dim(x = harmony.embed))
base::rownames(x = harmony.embed) <- base::rownames(x = lsi.embed)
base::colnames(x = harmony.embed) <- base::paste0("harmony_", 1:nsv)
base::cat("\n")
base::print(x = utils::head(x = harmony.embed[,1:6]))
sobj[["harmony"]] <- Seurat::CreateDimReducObject(embeddings = harmony.embed, assay = "ATAC", key = "harmony_")
PlotIntegration(seurat.obj = sobj, act.assay = "ATAC", features.to.plot = "harmony_1", save.dir = fig.int.dir, fn.suf = "_harmony1")

base::cat("\n")
base::print(x = sobj)
base::cat("\n")
base::print(x = sobj@reductions)
base::saveRDS(object = sobj, file = base::file.path(rds.dir, "after_integration.RDS"))

## Create a new UMAP using the integrated embeddings.
sobj <- Seurat::RunUMAP(object = sobj, reduction = "harmony", dims = ssv:nsv, verbose = F)

dp <- Seurat::DimPlot(object = sobj, reduction = red.to.plot, group.by = var.to.group) + 
  Seurat::FontSize(x.text = f, y.text = f, x.title = f, y.title = f)
grDevices::pdf(file = base::file.path(fig.int.dir, base::sprintf("%s_AI_G%s.pdf", red.to.plot, var.to.group)), 
               width = w1, height = h, onefile = T, paper = "special")
base::print(x = dp)
grDevices::dev.off()

sdp <- Seurat::DimPlot(object = sobj, reduction = red.to.plot, group.by = var.to.group, split.by = var.to.split) + 
  Seurat::FontSize(x.text = f, y.text = f, x.title = f, y.title = f) + 
  Seurat::NoLegend()
grDevices::pdf(file = base::file.path(fig.int.dir, base::sprintf("%s_AI_G%s_S%s.pdf", red.to.plot, var.to.group, var.to.split)), 
               width = w2, height = h, onefile = T, paper = "special")
base::print(x = sdp)
grDevices::dev.off()

sobj <- Seurat::FindNeighbors(object = sobj, reduction = "harmony", dims = ssv:nsv, nn.eps = 0, verbose = F)

## Quantify the activity of each gene in the genome by assessing the chromatin accessibility associated with each gene, 
## and create a new gene activity assay derived from the scATAC-seq data.
gene.activities <- Signac::GeneActivity(object = sobj, verbose = T)
## Add the gene activity matrix to the Seurat object as a new assay and normalize it.
sobj[["RNA"]] <- Seurat::CreateAssayObject(counts = gene.activities, min.cells = 0, min.features = 0)
Seurat::DefaultAssay(object = sobj) <- "RNA"
base::cat("\n")
base::print(x = sobj)
sobj <- Seurat::NormalizeData(
  object = sobj, 
  assay = "RNA", 
  normalization.method = "LogNormalize", 
  scale.factor = stats::median(x = sobj$nCount_RNA), 
  verbose = F
)
Seurat::DefaultAssay(object = sobj) <- "ATAC"
base::cat("\n")
base::print(x = sobj)
base::saveRDS(object = sobj, file = base::file.path(rds.dir, "after_neighbor.RDS"))

act.assay <- "ATAC"
alg.to.use <- 3
save.clusters <- T
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
