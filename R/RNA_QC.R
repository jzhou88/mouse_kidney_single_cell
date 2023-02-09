# RNA_QC.R
base::cat("\nRunning \"RNA_QC.R\" ...\n")

curr.mouse <- "FA1"
## Parse cell IDs
names.field <- 1
names.delim <- "_"
mt.pattern <- "^mt-"
cell.id.type <- 1
## QC criteria
min.umi <- 500
min.gene <- 250
max.gene <- 2500
max.mito <- 0.5
min.complex <- 0.25
## Run SoupX and/or DoubletFinder
run.soupx <- T
run.doubletfinder <- T
## Do normal QC
min.cells <- 10
do.qc <- T
## Cluster the cells
do.clustering.1 <- T
do.clustering.2 <- T
features.to.add <- NULL

## Set the paths to raw droplets and filtered counts, respectively.
curr.raw.dir <- base::file.path(cellranger.dir, curr.mouse, "raw_feature_bc_matrix")
curr.filtered.dir <- base::file.path(cellranger.dir, curr.mouse, "filtered_feature_bc_matrix")
## Load filtered counts data.
curr.counts <- Seurat::Read10X(data.dir = curr.filtered.dir)

## Profile the filtered counts data.
rvals <- StandardSeurat(counts = curr.counts, project = curr.mouse, assay = assay.to.use, names.field = names.field, 
                        names.delim = names.delim, mt.pattern = mt.pattern, 
                        cell.id.type = 0, # make sure the clustering info has the identical cell IDs as the original counts.
                        do.qc = F, # do not do QC here.
                        rds.before.qc = F, # do not save cells here.
                        fig.before.qc = F, fig.after.qc = F, 
                        do.clustering = do.clustering.1) # cluster the cells if required (e.g., by SoupX).
orig.ncell <- base::ncol(x = rvals$seurat.obj)
doublet.rate <- rvals$doublet.rate
orig.total.umi <- rvals$total.umi

if (run.soupx) {
  ## Remove ambient RNA using SoupX.
  cluster.labels <- rvals$cluster.label
  rvals <- NULL
  base::gc()
  soupx.corrected.counts <- base::file.path(rds.dir, base::sprintf(fmt = "soupx_corrected_counts_%s.RDS", curr.mouse))
  if (base::file.exists(soupx.corrected.counts)) {
    base::cat(base::sprintf(fmt = "\nReading %s...\n", soupx.corrected.counts))
    curr.counts <- base::readRDS(file = soupx.corrected.counts)
    
  } else {
    curr.counts <- StandardSoupX(raw.droplets.dir = curr.raw.dir, filtered.counts.dir = curr.filtered.dir, 
                                 cluster.labels = cluster.labels, save.fig.dir = qc.dir, 
                                 fn.suf = base::sprintf(fmt = "_%s", curr.mouse), forceAccept = T, roundToInt = T)
    base::saveRDS(object = curr.counts, file = soupx.corrected.counts)
  }
  ## Count the total number of UMIs after ambient RNA removal.
  rvals <- StandardSeurat(counts = curr.counts, project = curr.mouse, assay = assay.to.use, names.field = names.field, 
                          names.delim = names.delim, mt.pattern = mt.pattern, cell.id.type = 0, 
                          do.qc = F, # do not do QC here.
                          rds.before.qc = F, # do not save cells here.
                          fig.before.qc = F, fig.after.qc = F, 
                          do.clustering = F) # do not cluster cells here.
  corrected.total.umi <- rvals$total.umi
  ambient.rna <- orig.total.umi - corrected.total.umi
  base::cat(base::sprintf(fmt = "\n%s ambient.RNA %s %s %s\n", curr.mouse, ambient.rna, orig.total.umi, ambient.rna / orig.total.umi))
}
rvals <- NULL
base::gc()

## Do the gene- and cell-level QC.
rvals <- StandardSeurat(counts = curr.counts, project = curr.mouse, assay = assay.to.use, names.field = names.field, 
                        names.delim = names.delim, mt.pattern = mt.pattern, 
                        min.cells = min.cells, # QC: gene-level filtering
                        cell.id.type = cell.id.type, # each cell gets a new unique cell ID (i.e., sample name + original cell ID).
                        do.qc = do.qc, # QC: cell-level filtering
                        min.umi = min.umi, min.gene = min.gene, max.gene = max.gene, max.mito = max.mito, 
                        min.complex = min.complex, 
                        rds.before.qc = T, # save the cells before cell-level QC.
                        save.rds.dir = rds.dir, fig.before.qc = F, fig.after.qc = F, fn.suf = base::sprintf(fmt = "_%s", curr.mouse), 
                        do.clustering = do.clustering.2) # cluster the cells if required (e.g., by DoubletFinder).
sobj <- rvals$seurat.obj
rvals <- NULL
base::gc()
filtered.cell.ids <- base::as.character(x = sobj$cellID)
filtered.ncell <- base::ncol(x = sobj)

if (run.doubletfinder) {
  ## Remove doublets using DoubletFinder.
  ncell <- base::ncol(x = sobj)
  filtered.cell.ids <- StandardDoubletFinder(seurat.obj = sobj, doublet.rate = doublet.rate, save.fig.dir = qc.dir, 
                                             fn.suf = base::sprintf(fmt = "_%s", curr.mouse), save.rst.dir = rst.dir)
  sobj <- subset(x = sobj, subset = (cellID %in% filtered.cell.ids))
  filtered.ncell <- base::ncol(x = sobj)
  if (base::length(x = filtered.cell.ids) != filtered.ncell) {
    base::stop(base::sprintf("%s did not have the right number of filtered cells.", curr.mouse))
  }
  ndoublet <- ncell - filtered.ncell
  base::cat(base::sprintf(fmt = "\n%s found.doublets %s %s %s\n", curr.mouse, ndoublet, ncell, ndoublet / ncell))
}

## Save cell IDs of filtered cells.
base::saveRDS(object = filtered.cell.ids, file = base::file.path(rds.dir, base::sprintf(fmt = "filtered_cells_%s.RDS", curr.mouse)))
base::cat(base::sprintf(fmt = "\n%s filtered.cells %s %s %s\n", curr.mouse, filtered.ncell, orig.ncell, filtered.ncell / orig.ncell))
