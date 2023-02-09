# ATAC_QC.R
base::cat("\nRunning \"ATAC_QC.R\" ...\n")
base::options(future.globals.maxSize = ((1024 ^ 3) * 128)) # for 128 GB RAM
future::plan(strategy = "multicore", workers = req.cores)

min.prf <- 3000
max.prf <- 20000
min.prp <- 30
max.br <- 0.05
max.ns <- 4
min.tsse <- 2
min.pw <- 20
max.pw <- 10000
species.to.use <- "Mus_musculus"
blacklist.to.use <- blacklist_mm10
min.noc <- 20
npcs <- 60
sobj <- base::readRDS(file = base::file.path(curr.out.dir, "../HqAdult/RDS/after_QCcomp.RDS"))
peaks.called <- base::readRDS(file = base::file.path(curr.out.dir, "../HqAdult/RDS/called_peaks.RDS"))

base::cat(base::sprintf(fmt = "\nProcessing %s ...\n", picked.name))
base::cat("\n")
base::print(x = sobj)
## Remove cells that are outliers for QC metrics.
sobj <- subset(
  x = sobj, 
  subset = (nCount_ATAC > min.prf) & 
    (nCount_ATAC < max.prf) & 
    (pct_reads_in_peaks > min.prp) & 
    (blacklist_ratio < max.br) & 
    (nucleosome_signal < max.ns) & 
    (TSS.enrichment > min.tsse)
)
base::cat("\n")
base::print(x = sobj)

base::cat("\n")
base::print(x = peaks.called)
## Filter out bad peaks based on length.
peakwidths <- BiocGenerics::width(x = peaks.called)
peaks.called <- peaks.called[(peakwidths > min.pw) & (peakwidths < max.pw)]
## Remove peaks on nonstandard chromosomes and in genomic blacklist regions.
peaks.called <- GenomeInfoDb::keepStandardChromosomes(x = peaks.called, species = species.to.use, pruning.mode = "coarse")
peaks.called <- IRanges::subsetByOverlaps(x = peaks.called, ranges = blacklist.to.use, invert = T)
base::cat("\n")
base::print(x = peaks.called)
peaks.to.use <- Signac::GRangesToString(grange = peaks.called)
sobj <- subset(x = sobj, features = peaks.to.use)
base::cat("\n")
base::print(x = sobj)
base::saveRDS(object = sobj, file = base::file.path(rds.dir, "after_QC.RDS"))
PlotQualityMetrics(seurat.obj = sobj, save.fig.dir = fig.qc.dir, var.to.split = "sample", unsplit.box = T, fn.suf = "_after", all.metrics = T)

### Perform normalization and linear dimensional reduction.
## Normalization
sobj <- Signac::RunTFIDF(object = sobj, method = 1, scale.factor = 10000, verbose = T)
## Feature selection
sobj <- Signac::FindTopFeatures(object = sobj, min.cutoff = min.noc, verbose = T)
## Dimension reduction
sobj <- Signac::RunSVD(object = sobj, n = npcs, verbose = T)
base::saveRDS(object = sobj, file = base::file.path(rds.dir, "after_SVD.RDS"))
## Assess correlation between each LSI component and sequencing depth.
dc <- Signac::DepthCor(object = sobj, n = 20)
grDevices::pdf(file = base::file.path(fig.pc.dir, "DC.pdf"), width = 14, height = 7, onefile = T, paper = "special")
base::print(x = dc)
grDevices::dev.off()
