# auxFunc.R
base::cat("\nRunning \"auxFunc.R\" ...\n")

CalcClusterDaps <- function(seurat.obj, act.assay, var.as.cluster, save.rds.dir, save.rst.dir, known.clusters, sort.clusters, fn.pre, 
                            against.cluster = NULL, rename.clusters = NULL, logfc.thresh = 0.25, test.to.use = "LR", min.pct = 0.05, 
                            only.pos = F, latent.vars = NULL, p.thresh = 0.05, fn.suf = "", min.cells.per.cluster = 10) {
  if (!base::is.null(x = var.as.cluster)) {
    Seurat::Idents(object = seurat.obj) <- var.as.cluster
  }
  base::cat("\n")
  base::print(x = base::table(Seurat::Idents(object = seurat.obj)))
  clusters <- base::levels(x = seurat.obj)
  if (sort.clusters) {
    clusters <- known.clusters[base::sort(x = base::match(x = clusters, table = known.clusters), decreasing = F)]
  } else {
    clusters <- known.clusters
  }
  Seurat::Idents(object = seurat.obj) <- base::factor(x = Seurat::Idents(object = seurat.obj), levels = clusters)
  num.cells.per.cluster <- base::table(Seurat::Idents(object = seurat.obj))
  base::cat("\n")
  base::print(x = num.cells.per.cluster)
  valid.clusters <- base::names(x = num.cells.per.cluster[num.cells.per.cluster >= min.cells.per.cluster])
  base::cat("\n")
  base::print(x = valid.clusters)
  if ((!base::is.null(x = against.cluster)) && (!(against.cluster %in% valid.clusters))) {
    base::stop(base::sprintf(fmt = "less than %s cells from %s group, not enough for DAP calculation.", min.cells.per.cluster, against.cluster))
  }
  
  if (!base::is.null(x = against.cluster)) {
    clusters <- clusters[!(clusters %in% base::c(against.cluster))]
    base::cat("\n")
    base::print(x = clusters)
  }
  
  sl <- "data"
  mdp <- -Inf
  
  cluster.daps <- base::list()
  cluster.p.daps <- base::list()
  for (cluster in clusters) {
    cluster.name <- cluster
    if (!base::is.null(x = rename.clusters)) {
      cluster.name <- base::sprintf(fmt = "%s %s", rename.clusters, cluster)
    }
    if (!(cluster %in% valid.clusters)) {
      base::cat(base::sprintf(fmt = "\nWarning: less than %s cells from %s group, not enough for DAP calculation.\n", min.cells.per.cluster, cluster.name))
      base::cat(base::sprintf(fmt = "%s skipped.\n", cluster.name))
      next
    }
    base::cat(base::sprintf(fmt = "\nProcessing %s...\n", cluster.name))
    tmp <- Seurat::FindMarkers(object = seurat.obj, ident.1 = cluster, ident.2 = against.cluster, assay = act.assay, slot = sl, 
                               logfc.threshold = logfc.thresh, test.use = test.to.use, min.pct = min.pct, min.diff.pct = mdp, 
                               verbose = F, only.pos = only.pos, latent.vars = latent.vars)
    tmp$diff.pct <- tmp$pct.1 - tmp$pct.2
    cluster.daps[[cluster.name]] <- tmp[base::order(tmp$avg_log2FC, decreasing = T),]
    
    p.tmp <- tmp[base::which(tmp$p_val_adj < p.thresh),]
    cluster.p.daps[[cluster.name]] <- p.tmp[base::order(p.tmp$avg_log2FC, decreasing = T),]
  }
  
  base::saveRDS(object = cluster.daps, file = base::file.path(save.rds.dir, base::sprintf(fmt = "%s_%s_DAPs%s.RDS", fn.pre, test.to.use, fn.suf)))
  openxlsx::write.xlsx(x = cluster.daps, file = base::file.path(save.rst.dir, base::sprintf(fmt = "%s_%s_DAPs%s.xlsx", fn.pre, test.to.use, fn.suf)), 
                       startCol = 1, startRow = 1, colNames = T, rowNames = T, overwrite = T)
  openxlsx::write.xlsx(x = cluster.p.daps, 
                       file = base::file.path(save.rst.dir, base::sprintf(fmt = "%s_%s_DAPs%s_p%s.xlsx", fn.pre, test.to.use, fn.suf, p.thresh)), 
                       startCol = 1, startRow = 1, colNames = T, rowNames = T, overwrite = T)
  
  base::return(cluster.daps)
}

CalcClusterDegs <- function(seurat.obj, act.assay, var.as.cluster, save.rds.dir, save.rst.dir, save.fig.dir, known.clusters, sort.clusters, 
                            fn.pre, against.cluster = NULL, rename.clusters = NULL, logfc.thresh = 0.25, test.to.use = "MAST", min.pct = 0.1, 
                            only.pos = F, p.thresh = 0.05, do.dot = T, do.heatmap = T, top.genes = 20, cells.per.cluster = 100, fn.suf = "", 
                            organism = "Mus musculus", min.cells.per.cluster = 10) {
  if (!base::is.null(x = var.as.cluster)) {
    Seurat::Idents(object = seurat.obj) <- var.as.cluster
  }
  base::cat("\n")
  base::print(x = base::table(Seurat::Idents(object = seurat.obj)))
  clusters <- base::levels(x = seurat.obj)
  if (sort.clusters) {
    clusters <- known.clusters[base::sort(x = base::match(x = clusters, table = known.clusters), decreasing = F)]
  } else {
    clusters <- known.clusters
  }
  Seurat::Idents(object = seurat.obj) <- base::factor(x = Seurat::Idents(object = seurat.obj), levels = clusters)
  num.cells.per.cluster <- base::table(Seurat::Idents(object = seurat.obj))
  base::cat("\n")
  base::print(x = num.cells.per.cluster)
  valid.clusters <- base::names(x = num.cells.per.cluster[num.cells.per.cluster >= min.cells.per.cluster])
  base::cat("\n")
  base::print(x = valid.clusters)
  if ((!base::is.null(x = against.cluster)) && (!(against.cluster %in% valid.clusters))) {
    base::stop(base::sprintf(fmt = "less than %s cells from %s group, not enough for DEG calculation.", min.cells.per.cluster, against.cluster))
  }
  
  if (!base::is.null(x = against.cluster)) {
    clusters <- clusters[!(clusters %in% base::c(against.cluster))]
    base::cat("\n")
    base::print(x = clusters)
  }
  
  sl <- "data"
  mdp <- -Inf
  
  cluster.degs <- base::list()
  cluster.p.degs <- base::list()
  genes.to.plot <- NULL
  for (cluster in clusters) {
    cluster.name <- cluster
    if (!base::is.null(x = rename.clusters)) {
      cluster.name <- base::sprintf(fmt = "%s %s", rename.clusters, cluster)
    }
    if (!(cluster %in% valid.clusters)) {
      base::cat(base::sprintf(fmt = "\nWarning: less than %s cells from %s group, not enough for DEG calculation.\n", min.cells.per.cluster, cluster.name))
      base::cat(base::sprintf(fmt = "%s skipped.\n", cluster.name))
      next
    }
    base::cat(base::sprintf(fmt = "\nProcessing %s ...\n", cluster.name))
    tmp <- Seurat::FindMarkers(object = seurat.obj, ident.1 = cluster, ident.2 = against.cluster, assay = act.assay, slot = sl, 
                               logfc.threshold = logfc.thresh, test.use = test.to.use, min.pct = min.pct, min.diff.pct = mdp, 
                               verbose = F, only.pos = only.pos)
    tmp$diff.pct <- tmp$pct.1 - tmp$pct.2
    cluster.degs[[cluster.name]] <- tmp[base::order(tmp$avg_log2FC, decreasing = T),]
    
    p.tmp <- tmp[base::which(tmp$p_val_adj < p.thresh),]
    cluster.p.degs[[cluster.name]] <- p.tmp[base::order(p.tmp$avg_log2FC, decreasing = T),]
    
    genes.to.plot <- base::c(genes.to.plot, utils::head(x = base::row.names(x = cluster.p.degs[[cluster.name]]), n = top.genes))
  }
  genes.to.plot <- base::unique(x = genes.to.plot)
  
  base::saveRDS(object = cluster.degs, file = base::file.path(save.rds.dir, base::sprintf(fmt = "%s_%s_DEGs%s.RDS", fn.pre, test.to.use, fn.suf)))
  openxlsx::write.xlsx(x = cluster.degs, file = base::file.path(save.rst.dir, base::sprintf(fmt = "%s_%s_DEGs%s.xlsx", fn.pre, test.to.use, fn.suf)), 
                       startCol = 1, startRow = 1, colNames = T, rowNames = T, overwrite = T)
  openxlsx::write.xlsx(x = cluster.p.degs, 
                       file = base::file.path(save.rst.dir, base::sprintf(fmt = "%s_%s_DEGs%s_p%s.xlsx", fn.pre, test.to.use, fn.suf, p.thresh)), 
                       startCol = 1, startRow = 1, colNames = T, rowNames = T, overwrite = T)
  
  if (do.dot) {
    d <- 18
    f <- 18
    ## The values in DotPlot are extracted from the @data slot, averaged, and then passed to scale.
    ## These are then Min-Maxed based on the col.min and col.max parameter values.
    ## The DotPlot shows scaled values (which can be both positive and negative), 
    ## so that we can visualize highly and lowly expressed genes on the same color scale.
    dp <- Seurat::DotPlot(object = seurat.obj, assay = act.assay, features = genes.to.plot, dot.scale = d, cluster.idents = F, scale = T) + 
      Seurat::RotatedAxis() + 
      viridis::scale_color_viridis(option = "D", direction = -1) + 
      Seurat::FontSize(x.text = f, y.text = f, x.title = f, y.title = f) + 
      ggplot2::theme(axis.line = ggplot2::element_blank(), 
                     panel.grid = ggplot2::element_line(colour = "gray", size = 0.5, linetype = "dashed"), 
                     panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1))
    w <- 0.5 * base::length(x = genes.to.plot) + 3
    if (7 > w) {
      w <- 7
    }
    h <- 0.7 * base::length(x = clusters)
    if (7 > h) {
      h <- 7
    }
    grDevices::pdf(file = base::file.path(save.fig.dir, base::sprintf(fmt = "dotplot_%s_%s_top%smarkers%s_p%s.pdf", 
                                                                      fn.pre, test.to.use, top.genes, fn.suf, p.thresh)), 
                   width = w, height = h, onefile = T, paper = "special")
    base::print(x = dp)
    grDevices::dev.off()
  }
  
  ### Plot heatmap
  if (do.heatmap) {
    seurat.obj <- Seurat::ScaleData(object = seurat.obj, assay = act.assay, features = genes.to.plot, verbose = F)
    hm <- NULL
    w <- NULL
    if (0 < cells.per.cluster) {
      hm <- Seurat::DoHeatmap(object = subset(x = seurat.obj, downsample = cells.per.cluster, seed = 1), features = genes.to.plot, 
                              group.by = "ident", slot = "scale.data", assay = act.assay, label = T, size = 20)
      w <- 0.06 * cells.per.cluster * base::length(x = clusters) + 3
    } else {
      hm <- Seurat::DoHeatmap(object = seurat.obj, features = genes.to.plot, group.by = "ident", slot = "scale.data", 
                              assay = act.assay, label = T, size = 20)
      w <- 0.06 * base::ncol(x = seurat.obj) + 3
    }
    h <- 0.15 * base::length(x = genes.to.plot)
    grDevices::pdf(file = base::file.path(save.fig.dir, base::sprintf(fmt = "heatmap_%s_%s_top%smarkers%s_p%s.pdf", 
                                                                      fn.pre, test.to.use, top.genes, fn.suf, p.thresh)), 
                   width = w, height = h, onefile = T, paper = "special")
    base::print(x = hm)
    grDevices::dev.off()
  }
  
  ### Add gene annotations
  annotations <- FetchGeneAnnotations(save.rds.dir = save.rds.dir, save.rst.dir = save.rst.dir, organism = organism)
  for (i in base::seq(from = 1, to = base::length(x = cluster.p.degs), by = 1)) {
    ## Combine DEGs with gene descriptions
    cluster.p.degs[[i]] <- cluster.p.degs[[i]] %>% 
      tibble::rownames_to_column(var="gene") %>% 
      dplyr::left_join(y = base::unique(annotations[, base::c("gene_name", "description")]), by = c("gene" = "gene_name")) 
  }
  openxlsx::write.xlsx(x = cluster.p.degs, file = base::file.path(save.rst.dir, base::sprintf(fmt = "%s_%s_annotated_DEGs%s_p%s.xlsx", 
                                                                                              fn.pre, test.to.use, fn.suf, p.thresh)), 
                       startCol = 1, startRow = 1, colNames = T, rowNames = F, overwrite = T)
  
  base::return(cluster.degs)
}

CalcCosgMarkers <- function(seurat.obj, known.clusters, save.rds.dir, save.rst.dir, save.fig.dir, fn.pre, var.as.cluster = NULL, 
                            clusters.to.inc = "all", act.assay = "RNA", markers.to.calc = 100, fn.suf = "", top.markers = 20, 
                            cells.per.cluster = 100) {
  if (!base::is.null(x = var.as.cluster)) {
    Seurat::Idents(object = seurat.obj) <- var.as.cluster
  }
  base::cat("\n")
  base::print(x = base::table(Seurat::Idents(object = seurat.obj)))
  clusters <- base::levels(x = seurat.obj)
  clusters <- known.clusters[base::sort(x = base::match(x = clusters, table = known.clusters), decreasing = F)]
  Seurat::Idents(object = seurat.obj) <- base::factor(x = Seurat::Idents(object = seurat.obj), levels = clusters)
  base::cat("\n")
  base::print(x = base::table(Seurat::Idents(object = seurat.obj)))
  
  cosg.markers <- COSG::cosg(object = seurat.obj, groups = clusters.to.inc, assay = act.assay, 
                             slot = "data", mu = 1, n_genes_user = markers.to.calc)
  base::saveRDS(object = cosg.markers, file = base::file.path(save.rds.dir, base::sprintf(fmt = "%s_COSG_markers%s.RDS", fn.pre, fn.suf)))
  openxlsx::write.xlsx(x = cosg.markers, file = base::file.path(save.rst.dir, base::sprintf(fmt = "%s_COSG_markers%s.xlsx", fn.pre, fn.suf)), 
                       startCol = 1, startRow = 1, colNames = T, rowNames = F, overwrite = T)
  
  marker.names <- cosg.markers$names
  base::cat("\n")
  base::print(x = utils::head(x = marker.names))
  
  markers.to.plot <- NULL
  for (cluster in base::colnames(x = marker.names)) {
    cluster.markers <- marker.names[1:top.markers, cluster]
    markers.to.plot <- base::c(markers.to.plot, cluster.markers)
  }
  markers.to.plot <- base::unique(x = markers.to.plot)
  
  ### Dot plot
  d <- 18
  f <- 18
  ## The values in DotPlot are extracted from the @data slot, averaged, and then passed to scale.
  ## These are then Min-Maxed based on the col.min and col.max parameter values.
  ## The DotPlot shows scaled values (which can be both positive and negative), 
  ## so that we can visualize highly and lowly expressed genes on the same color scale.
  dp <- Seurat::DotPlot(object = seurat.obj, assay = act.assay, features = markers.to.plot, dot.scale = d, cluster.idents = F, scale = T) + 
    Seurat::RotatedAxis() + 
    viridis::scale_color_viridis(option = "D", direction = -1) + 
    Seurat::FontSize(x.text = f, y.text = f, x.title = f, y.title = f) + 
    ggplot2::theme(axis.line = ggplot2::element_blank(), 
                   panel.grid = ggplot2::element_line(colour = "gray", size = 0.5, linetype = "dashed"), 
                   panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1))
  w <- 0.5 * base::length(x = markers.to.plot) + 3
  if (7 > w) {
    w <- 7
  }
  h <- 0.7 * base::ncol(x = marker.names)
  if (7 > h) {
    h <- 7
  }
  grDevices::pdf(file = base::file.path(save.fig.dir, base::sprintf(fmt = "dotplot_%s_COSG_top%smarkers%s.pdf", fn.pre, top.markers, fn.suf)), 
                 width = w, height = h, onefile = T, paper = "special")
  base::print(x = dp)
  grDevices::dev.off()
  
  ### Clustered dot plot
  dp <- Seurat::DotPlot(object = seurat.obj, assay = act.assay, features = markers.to.plot, dot.scale = d, cluster.idents = T, scale = T) + 
    Seurat::RotatedAxis() + 
    viridis::scale_color_viridis(option = "D", direction = -1) + 
    Seurat::FontSize(x.text = f, y.text = f, x.title = f, y.title = f) + 
    ggplot2::theme(axis.line = ggplot2::element_blank(), 
                   panel.grid = ggplot2::element_line(colour = "gray", size = 0.5, linetype = "dashed"), 
                   panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1))
  grDevices::pdf(file = base::file.path(save.fig.dir, base::sprintf(fmt = "dotplot_clustered_%s_COSG_top%smarkers%s.pdf", 
                                                                    fn.pre, top.markers, fn.suf)), 
                 width = w, height = h, onefile = T, paper = "special")
  base::print(x = dp)
  grDevices::dev.off()
  
  ### Heatmap
  seurat.obj <- Seurat::ScaleData(object = seurat.obj, assay = act.assay, features = markers.to.plot, verbose = F)
  hm <- NULL
  w <- NULL
  if (0 < cells.per.cluster) {
    hm <- Seurat::DoHeatmap(object = subset(x = seurat.obj, downsample = cells.per.cluster, seed = 1), features = markers.to.plot, 
                            group.by = "ident", slot = "scale.data", assay = act.assay, label = T, size = 20)
    w <- 0.06 * cells.per.cluster * base::ncol(x = marker.names) + 3
  } else {
    hm <- Seurat::DoHeatmap(object = seurat.obj, features = markers.to.plot, group.by = "ident", slot = "scale.data", 
                            assay = act.assay, label = T, size = 20)
    w <- 0.06 * base::ncol(x = seurat.obj) + 3
  }
  h <- 0.15 * base::length(x = markers.to.plot)
  grDevices::pdf(file = base::file.path(save.fig.dir, base::sprintf(fmt = "heatmap_%s_COSG_top%smarkers%s.pdf", fn.pre, top.markers, fn.suf)), 
                 width = w, height = h, onefile = T, paper = "special")
  base::print(x = hm)
  grDevices::dev.off()
  
  base::return(cosg.markers)
}

CalcSampleInClusterDaps <- function(seurat.obj, act.assay, var.as.cluster, var.as.sample, save.rds.dir, save.rst.dir, known.clusters, 
                                    sort.clusters, fn.pre, rename.clusters = NULL, samples.to.inc = NULL, against.sample = NULL, 
                                    logfc.thresh = 0.25, test.to.use = "LR", min.pct = 0.05, only.pos = F, latent.vars = NULL, p.thresh = 0.05, 
                                    fn.suf = "", min.cells.per.sample = 10) {
  if (base::is.null(x = samples.to.inc)) {
    Seurat::Idents(object = seurat.obj) <- var.as.sample
    samples.to.inc <- base::levels(x = seurat.obj)
  }
  if ((!base::is.null(x = against.sample)) && (against.sample %in% samples.to.inc)) {
    samples.to.inc <- samples.to.inc[!(samples.to.inc %in% base::c(against.sample))]
  }
  base::cat("\n")
  base::print(x = samples.to.inc)
  Seurat::Idents(object = seurat.obj) <- var.as.cluster
  base::cat("\n")
  base::print(x = base::table(Seurat::Idents(object = seurat.obj)))
  clusters <- base::levels(x = seurat.obj)
  if (sort.clusters) {
    clusters <- known.clusters[base::sort(x = base::match(x = clusters, table = known.clusters), decreasing = F)]
  } else {
    clusters <- known.clusters
  }
  Seurat::Idents(object = seurat.obj) <- base::factor(x = Seurat::Idents(object = seurat.obj), levels = clusters)
  base::cat("\n")
  base::print(x = base::table(Seurat::Idents(object = seurat.obj)))
  split.sobj <- Seurat::SplitObject(object = seurat.obj, split.by = var.as.cluster)
  
  sl <- "data"
  mdp <- -Inf
  
  daps <- base::list()
  p.daps <- base::list()
  for (cluster in clusters) {
    cluster.name <- cluster
    if (!base::is.null(x = rename.clusters)) {
      cluster.name <- base::sprintf(fmt = "%s %s", rename.clusters, cluster)
    }
    base::cat(base::sprintf(fmt = "\nProcessing %s ...\n", cluster.name))
    seurat.obj <- split.sobj[[cluster]]
    base::cat("\n")
    base::print(x = base::table(Seurat::Idents(object = seurat.obj)))
    Seurat::Idents(object = seurat.obj) <- var.as.sample
    cells.per.sample <- base::table(Seurat::Idents(object = seurat.obj))
    base::cat("\n")
    base::print(x = cells.per.sample)
    valid.samples <- base::names(x = cells.per.sample[cells.per.sample >= min.cells.per.sample])
    base::cat("\n")
    base::print(x = valid.samples)
    if ((!base::is.null(x = against.sample)) && (!(against.sample %in% valid.samples))) {
      base::cat(base::sprintf(fmt = "\nWarning: less than %s cells from %s group, not enough for DAP calculation.\n", min.cells.per.sample, against.sample))
      base::cat(base::sprintf(fmt = "%s skipped.\n", cluster.name))
      next
    }
    for (sample in samples.to.inc) {
      curr.name <- base::sprintf(fmt = "%s in %s", sample, cluster.name)
      if (!(sample %in% valid.samples)) {
        base::cat(base::sprintf(fmt = "\nWarning: less than %s cells from %s group, not enough for DAP calculation.\n", min.cells.per.sample, curr.name))
        base::cat(base::sprintf(fmt = "%s skipped.\n", curr.name))
        next
      }
      base::cat(base::sprintf(fmt = "\nProcessing %s ...\n", curr.name))
      tmp <- Seurat::FindMarkers(object = seurat.obj, ident.1 = sample, ident.2 = against.sample, assay = act.assay, slot = sl, 
                                 logfc.threshold = logfc.thresh, test.use = test.to.use, min.pct = min.pct, min.diff.pct = mdp, 
                                 verbose = F, only.pos = only.pos, latent.vars = latent.vars)
      tmp$diff.pct <- tmp$pct.1 - tmp$pct.2
      daps[[curr.name]] <- tmp[base::order(tmp$avg_log2FC, decreasing = T),]
      
      p.tmp <- tmp[base::which(tmp$p_val_adj < p.thresh),]
      p.daps[[curr.name]] <- p.tmp[base::order(p.tmp$avg_log2FC, decreasing = T),]
    }
  }
  base::saveRDS(object = daps, file = base::file.path(save.rds.dir, base::sprintf(fmt = "%s_%s_DAPs%s.RDS", fn.pre, test.to.use, fn.suf)))
  openxlsx::write.xlsx(x = daps, file = base::file.path(save.rst.dir, base::sprintf(fmt = "%s_%s_DAPs%s.xlsx", fn.pre, test.to.use, fn.suf)), 
                       startCol = 1, startRow = 1, colNames = T, rowNames = T, overwrite = T)
  openxlsx::write.xlsx(x = p.daps, file = base::file.path(save.rst.dir, base::sprintf(fmt = "%s_%s_DAPs%s_p%s.xlsx", fn.pre, test.to.use, fn.suf, p.thresh)), 
                       startCol = 1, startRow = 1, colNames = T, rowNames = T, overwrite = T)
  
  base::return(daps)
}

CalcSampleInClusterDegs <- function(seurat.obj, act.assay, var.as.cluster, var.as.sample, save.rds.dir, save.rst.dir, known.clusters, 
                                    sort.clusters, rename.clusters, fn.pre, samples.to.inc = NULL, against.sample = NULL, logfc.thresh = 0.25, 
                                    test.to.use = "MAST", min.pct = 0.1, only.pos = F, p.thresh = 0.05, fn.suf = "", organism = "Mus musculus") {
  if (base::is.null(x = samples.to.inc)) {
    Seurat::Idents(object = seurat.obj) <- var.as.sample
    samples.to.inc <- base::levels(x = seurat.obj)
    base::cat("\n")
    base::print(x = samples.to.inc)
  }
  if ((!base::is.null(x = against.sample)) && (against.sample %in% samples.to.inc)) {
    samples.to.inc <- samples.to.inc[!(samples.to.inc %in% base::c(against.sample))]
    base::cat("\n")
    base::print(x = samples.to.inc)
  }
  Seurat::Idents(object = seurat.obj) <- var.as.cluster
  base::cat("\n")
  base::print(x = base::table(Seurat::Idents(object = seurat.obj)))
  clusters <- base::levels(x = seurat.obj)
  if (sort.clusters) {
    clusters <- known.clusters[base::sort(x = base::match(x = clusters, table = known.clusters), decreasing = F)]
  } else {
    clusters <- known.clusters
  }
  Seurat::Idents(object = seurat.obj) <- base::factor(x = Seurat::Idents(object = seurat.obj), levels = clusters)
  base::cat("\n")
  base::print(x = base::table(Seurat::Idents(object = seurat.obj)))
  split.sobj <- Seurat::SplitObject(object = seurat.obj, split.by = var.as.cluster)
  
  sl <- "data"
  mdp <- -Inf
  
  degs <- base::list()
  p.degs <- base::list()
  for (cluster in clusters) {
    cluster.name <- cluster
    if (rename.clusters) {
      cluster.name <- base::sprintf(fmt = "Cluster %s", cluster)
    }
    base::cat(base::sprintf(fmt = "\nProcessing %s ...\n", cluster.name))
    seurat.obj <- NULL
    base::gc()
    seurat.obj <- split.sobj[[cluster]]
    base::cat("\n")
    base::print(x = base::table(Seurat::Idents(object = seurat.obj)))
    Seurat::Idents(object = seurat.obj) <- var.as.sample
    base::cat("\n")
    base::print(x = base::table(Seurat::Idents(object = seurat.obj)))
    
    for (sample in samples.to.inc) {
      curr.name <- base::sprintf(fmt = "%s in %s", sample, cluster.name)
      base::cat(base::sprintf(fmt = "\nProcessing %s ...\n", curr.name))
      tmp <- Seurat::FindMarkers(object = seurat.obj, ident.1 = sample, ident.2 = against.sample, assay = act.assay, slot = sl, 
                                 logfc.threshold = logfc.thresh, test.use = test.to.use, min.pct = min.pct, min.diff.pct = mdp, 
                                 verbose = F, only.pos = only.pos)
      tmp$diff.pct <- tmp$pct.1 - tmp$pct.2
      degs[[curr.name]] <- tmp[base::order(tmp$avg_log2FC, decreasing = T),]
      
      p.tmp <- tmp[base::which(tmp$p_val_adj < p.thresh),]
      p.degs[[curr.name]] <- p.tmp[base::order(p.tmp$avg_log2FC, decreasing = T),]
    }
  }
  base::saveRDS(object = degs, file = base::file.path(save.rds.dir, base::sprintf(fmt = "%s_%s_DEGs%s.RDS", fn.pre, test.to.use, fn.suf)))
  openxlsx::write.xlsx(x = degs, file = base::file.path(save.rst.dir, base::sprintf(fmt = "%s_%s_DEGs%s.xlsx", fn.pre, test.to.use, fn.suf)), 
                       startCol = 1, startRow = 1, colNames = T, rowNames = T, overwrite = T)
  openxlsx::write.xlsx(x = p.degs, 
                       file = base::file.path(save.rst.dir, base::sprintf(fmt = "%s_%s_DEGs%s_p%s.xlsx", fn.pre, test.to.use, fn.suf, p.thresh)), 
                       startCol = 1, startRow = 1, colNames = T, rowNames = T, overwrite = T)
  
  ### Add gene annotations
  annotations <- FetchGeneAnnotations(save.rds.dir = save.rds.dir, save.rst.dir = save.rst.dir, organism = organism)
  for (i in base::seq(from = 1, to = base::length(x = p.degs), by = 1)) {
    ## Combine markers with gene descriptions
    p.degs[[i]] <- p.degs[[i]] %>% 
      tibble::rownames_to_column(var = "gene") %>% 
      dplyr::left_join(y = base::unique(annotations[, base::c("gene_name", "description")]), by = c("gene" = "gene_name"))
  }
  openxlsx::write.xlsx(x = p.degs, file = base::file.path(save.rst.dir, base::sprintf(fmt = "%s_%s_annotated_DEGs%s_p%s.xlsx", 
                                                                                      fn.pre, test.to.use, fn.suf, p.thresh)), 
                       startCol = 1, startRow = 1, colNames = T, rowNames = F, overwrite = T)
  
  base::return(degs)
}

DimPlotCells <- function(seurat.obj, save.fig.dir, red.to.plot = "umap", fn.suf = "", var.as.cell = "cell.type", sort.cell.types = F, 
                         known.cell.types = NULL, plot.label = T, no.legend = T, split.cells = F, do.split = F, var.to.split = "sample", 
                         w1 = NULL, w2 = NULL, w3 = NULL) {
  Seurat::Idents(object = seurat.obj) <- var.as.cell
  cell.types <- base::levels(x = seurat.obj)
  if (sort.cell.types) {
    cell.types <- known.cell.types[base::sort(x = base::match(x = cell.types, table = known.cell.types), decreasing = F)]
    Seurat::Idents(object = seurat.obj) <- base::factor(x = Seurat::Idents(object = seurat.obj), levels = cell.types)
  }
  
  l <- 5
  f <- 18
  h <- 7
  
  if (base::is.null(x = w1)) {
    w1 <- 7
    if (!no.legend) {
      w1 <- w1 + 2
    }
  }
  if (base::is.null(x = w2)) {
    w2 <- 7 * base::length(x = base::levels(x = seurat.obj))
  }
  if (base::is.null(x = w3)) {
    w3 <- 7 * base::length(x = base::unique(x = seurat.obj@meta.data[[var.to.split]]))
  }
  
  dp <- Seurat::DimPlot(object = seurat.obj, reduction = red.to.plot, label = plot.label, label.size = l, repel = T) + 
    Seurat::FontSize(x.text = f, y.text = f, x.title = f, y.title = f)
  if (no.legend) {
    dp <- dp + Seurat::NoLegend()
  }
  grDevices::pdf(file = base::file.path(save.fig.dir, base::sprintf(fmt = "%s_%s%s.pdf", red.to.plot, var.as.cell, fn.suf)), 
                 width = w1, height = h, onefile = T, paper = "special")
  base::print(x = dp)
  grDevices::dev.off()
  
  if (split.cells) {
    sdp <- Seurat::DimPlot(object = seurat.obj, reduction = red.to.plot, label = plot.label, label.size = l, repel = T, 
                           split.by = var.as.cell) + 
      Seurat::FontSize(x.text = f, y.text = f, x.title = f, y.title = f) + 
      Seurat::NoLegend()
    grDevices::pdf(file = base::file.path(save.fig.dir, base::sprintf(fmt = "%s_split%s%s.pdf", red.to.plot, var.as.cell, fn.suf)),
                   width = w2, height = h, onefile = T, paper = "special")
    base::print(x = sdp)
    grDevices::dev.off()
  }
  
  if (do.split) {
    sdp <- Seurat::DimPlot(object = seurat.obj, reduction = red.to.plot, label = plot.label, label.size = l, repel = T, 
                           split.by = var.to.split) + 
      Seurat::FontSize(x.text = f, y.text = f, x.title = f, y.title = f) + 
      Seurat::NoLegend()
    grDevices::pdf(file = base::file.path(save.fig.dir, base::sprintf(fmt = "%s_%s_split%s%s.pdf", 
                                                                      red.to.plot, var.as.cell, var.to.split, fn.suf)), 
                   width = w3, height = h, onefile = T, paper = "special")
    base::print(x = sdp)
    grDevices::dev.off()
  }
}

DotPlotCells <- function(seurat.obj, save.fig.dir, act.assay, genes.to.plot, fn.suf = "", var.as.cell = "cell.type", sort.cell.types = F, 
                         known.cell.types = NULL, rotate.axis = T, do.split = F, split.obj = F, var.to.split = "sample", colors.to.use = NULL) {
  Seurat::Idents(object = seurat.obj) <- var.as.cell
  cell.types <- base::levels(x = seurat.obj)
  if (sort.cell.types) {
    cell.types <- known.cell.types[base::sort(x = base::match(x = cell.types, table = known.cell.types), decreasing = F)]
    Seurat::Idents(object = seurat.obj) <- base::factor(x = Seurat::Idents(object = seurat.obj), levels = cell.types)
  }
  
  num.celltypes <- base::length(x = cell.types)
  h <- 0.7 * num.celltypes
  add.h <- 2
  if (!rotate.axis) {
    h <- h + add.h
  }
  if (7 > h) {
    h <- 7
  }
  w <- 0.5 * base::length(x = genes.to.plot) + 2
  if (7 > w) {
    w <- 7
  }
  c <- base::c("lightgrey", "blue")
  d <- 18
  f <- 18
  ## The values in DotPlot are extracted from the @data slot, averaged, and then passed to scale.
  ## These are then Min-Maxed based on the col.min and col.max parameter values.
  ## The DotPlot shows scaled values (which can be both positive and negative), 
  ## so that we can visualize highly and lowly expressed genes on the same color scale.
  dp <- Seurat::DotPlot(object = seurat.obj, assay = act.assay, features = genes.to.plot, 
                        cols = c, dot.scale = d, cluster.idents = F, scale = T) + 
    Seurat::FontSize(x.text = f, y.text = f, x.title = f, y.title = f)
  if (rotate.axis) {
    dp <- dp + Seurat::RotatedAxis()
  }
  grDevices::pdf(file = base::file.path(save.fig.dir, base::sprintf(fmt = "dotplot%s.pdf", fn.suf)), 
                 width = w, height = h, onefile = T, paper = "special")
  base::print(x = dp)
  grDevices::dev.off()
  
  if (do.split) {
    vars <- base::unique(x = seurat.obj@meta.data[[var.to.split]])
    num.vars <- base::length(x = vars)
    sdp <- NULL
    if (split.obj) {
      split.sobj <- Seurat::SplitObject(object = seurat.obj, split.by = var.to.split)
      seurat.obj <- NULL
      base::gc()
      dot.plots <- base::list()
      for (var in vars) {
        dp <- Seurat::DotPlot(object = split.sobj[[var]], assay = act.assay, features = genes.to.plot, 
                              cols = c, dot.scale = d, cluster.idents = F, scale = T) + 
          ggplot2::labs(title = var) + 
          Seurat::FontSize(x.text = f, y.text = f, x.title = f, y.title = f, main = f) + 
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
        if (rotate.axis) {
          dp <- dp + Seurat::RotatedAxis()
        }
        dot.plots[[var]] <- dp
      }
      nc <- 1
      nr <- num.vars
      h <- 0.7 * num.celltypes
      if (!rotate.axis) {
        h <- h + add.h
      }
      h <- h * nr
      if (7 > h) {
        h <- 7
      }
      sdp <- cowplot::plot_grid(plotlist = dot.plots, nrow = nr, ncol = nc)
      
    } else {
      if (base::is.null(x = colors.to.use)) {
        colors.to.use <- (scales::hue_pal())(n = num.vars)
      }
      h <- 0.7 * num.celltypes
      h <- h * num.vars
      if (!rotate.axis) {
        h <- h + add.h
      }
      if (7 > h) {
        h <- 7
      }
      sdp <- Seurat::DotPlot(object = seurat.obj, assay = act.assay, features = genes.to.plot, cols = colors.to.use, 
                             dot.scale = d, split.by = var.to.split, cluster.idents = F, scale = T) + 
        Seurat::FontSize(x.text = f, y.text = f, x.title = f, y.title = f)
      if (rotate.axis) {
        sdp <- sdp + Seurat::RotatedAxis()
      }
    }
    grDevices::pdf(file = base::file.path(save.fig.dir, base::sprintf(fmt = "dotplot_split%s%s.pdf", var.to.split, fn.suf)), 
                   width = w, height = h, onefile = T, paper = "special")
    base::print(x = sdp)
    grDevices::dev.off()
  }
}

FetchGeneAnnotations <- function(save.rds.dir, save.rst.dir, organism = "Mus musculus") {
  organism.name <- base::gsub(pattern = " ", replacement = "_", x = organism)
  annotation.csv <- base::file.path(save.rst.dir, base::sprintf(fmt = "gene_annotations_%s.csv", organism.name))
  annotations <- NULL
  if (base::file.exists(annotation.csv)) {
    annotations <- utils::read.csv(file = annotation.csv, header = T, sep = ",", stringsAsFactors = F, fileEncoding = "UTF-8-BOM")
  } else {
    ## Connect to AnnotationHub
    ah <- AnnotationHub::AnnotationHub()
    ## Access the Ensembl database for organism
    ahDb <- AnnotationHub::query(ah, pattern = base::c(organism, "EnsDb"), ignore.case = TRUE)
    ## Acquire the latest annotation files
    id <- ahDb %>% 
      S4Vectors::mcols() %>% 
      base::rownames() %>% 
      utils::tail(n = 1)
    ## Download the appropriate Ensembldb database
    edb <- ah[[id]]
    ## Extract gene-level information from database
    annotations <- ensembldb::genes(edb, return.type = "data.frame")
    ## Select annotations of interest
    annotations <- annotations %>% 
      dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
    base::saveRDS(object = annotations, file = base::file.path(save.rds.dir, base::sprintf(fmt = "gene_annotations_%s.RDS", organism.name)))
    utils::write.csv(x = annotations, file = annotation.csv, row.names = F)
  }
  base::return(annotations)
}

Mouse2HumanGenes <- function(genes.to.map, save.rst.dir = NULL, invert = F) {
  base::require(package = "biomaRt")
  
  human <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse <- biomaRt::useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  tmp <- NULL
  fn.pre <- NULL
  if (invert) {
    ## human to mouse
    tmp <- biomaRt::getLDS(attributes = base::c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), filters = "hgnc_symbol", 
                           values = genes.to.map, mart = human, attributesL = base::c("mgi_symbol", "chromosome_name", "start_position", "end_position"), 
                           martL = mouse, uniqueRows = T)
    tmp$index <- base::match(x = tmp$HGNC.symbol, table = genes.to.map)
    fn.pre <- "human_to_mouse"
    
  } else {
    ## mouse to human
    tmp <- biomaRt::getLDS(attributes = base::c("mgi_symbol", "chromosome_name", "start_position", "end_position"), filters = "mgi_symbol", 
                           values = genes.to.map, mart = mouse, attributesL = base::c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), 
                           martL = human, uniqueRows = T)
    tmp$index <- base::match(x = tmp$MGI.symbol, table = genes.to.map)
    fn.pre <- "mouse_to_human"
  }
  out.genes <- tmp[base::order(tmp$index, decreasing = F),]
  out.genes$index <- NULL
  
  if (!base::is.null(x = save.rst.dir)) {
    openxlsx::write.xlsx(x = out.genes, file = base::file.path(save.rst.dir, base::sprintf(fmt = "%s_genes.xlsx", fn.pre)), 
                         startCol = 1, startRow = 1, colNames = T, rowNames = F, overwrite = T)
  }
  
  base::return(out.genes)
}

PlotClusters <- function(seurat.obj, assay.clustered, resols, save.dir, red.to.plot = "umap", do.red = T, 
                         split.clusters = T, do.split = F, var.to.split = "sample", do.bar = F, var.to.bar = "sample", 
                         do.dot = T, assay.to.plot = "RNA", markers.to.plot = default.known.markers, fig.suf = "") {
  l <- 5
  f <- 18
  d <- 18
  w1 <- 7 * base::length(x = base::unique(x = seurat.obj@meta.data[[var.to.split]]))
  all.bar.vars <- base::unique(x = seurat.obj@meta.data[[var.to.bar]])
  num.bar.vars <- base::length(x = all.bar.vars)
  w2 <- num.bar.vars
  if (7 > w2) {
    w2 <- 7
  }
  w3 <- (0.5 * base::length(x = markers.to.plot)) + 2
  if (7 > w3) {
    w3 <- 7
  }
  
  for (res.to.use in resols) {
    clustering <- base::sprintf(fmt = "%s_snn_res.%s", assay.clustered, res.to.use)
    Seurat::Idents(object = seurat.obj) <- clustering
    num.clusters <- base::length(x = base::levels(x = seurat.obj))
    plot.title <- base::sprintf(fmt = "%s: %s clusters", clustering, num.clusters)
    clusters <- SortClusters(seurat.obj = seurat.obj)
    Seurat::Idents(object = seurat.obj) <- base::factor(x = Seurat::Idents(object = seurat.obj), levels = clusters)
    
    if (do.red) {
      rp <- Seurat::DimPlot(object = seurat.obj, reduction = red.to.plot, label = T, label.size = l, repel = T) + 
        ggplot2::labs(title = plot.title) + 
        Seurat::FontSize(x.text = f, y.text = f, x.title = f, y.title = f, main = f)
      grDevices::pdf(file = base::file.path(save.dir, 
                                            base::sprintf(fmt = "%splot1_%s_res%s%s.pdf", 
                                                          red.to.plot, assay.clustered, res.to.use, fig.suf)), 
                     width = 9, height = 7, onefile = T, paper = "special")
      base::print(x = rp)
      grDevices::dev.off()
    }
    if (split.clusters) {
      srp <- Seurat::DimPlot(object = seurat.obj, reduction = red.to.plot, label = T, label.size = l, repel = T, 
                             split.by = clustering) + 
        ggplot2::labs(title = plot.title) + 
        Seurat::NoLegend() + 
        Seurat::FontSize(x.text = f, y.text = f, x.title = f, y.title = f, main = f)
      grDevices::pdf(file = base::file.path(save.dir,
                                            base::sprintf(fmt = "%splot1_splitclusters_%s_res%s%s.pdf",
                                                          red.to.plot, assay.clustered, res.to.use, fig.suf)),
                     width = (7 * num.clusters), height = 7, onefile = T, paper = "special")
      base::print(x = srp)
      grDevices::dev.off()
    }
    if (do.split) {
      srp <- Seurat::DimPlot(object = seurat.obj, reduction = red.to.plot, label = T, label.size = l, repel = T, 
                             split.by = var.to.split) + 
        ggplot2::labs(title = plot.title) + 
        Seurat::NoLegend() + 
        Seurat::FontSize(x.text = f, y.text = f, x.title = f, y.title = f, main = f)
      grDevices::pdf(file = base::file.path(save.dir,
                                            base::sprintf(fmt = "%splot1_split%s_%s_res%s%s.pdf", red.to.plot, 
                                                          var.to.split, assay.clustered, res.to.use, fig.suf)),
                     width = w1, height = 7, onefile = T, paper = "special")
      base::print(x = srp)
      grDevices::dev.off()
    }
    
    if (do.bar) {
      cell.counts <- base::table(Seurat::Idents(object = seurat.obj), seurat.obj@meta.data[[var.to.bar]])
      ## for each cell type, number and fraction of cells from each sample
      bar.plots <- base::list()
      for (cluster in clusters) {
        cc <- cell.counts[cluster,]
        df <- base::data.frame(var = base::names(x = cc), count = cc)
        t <- base::sprintf(fmt = "Cluster %s: %s cells", cluster, base::sum(cc))
        bar.plots[[t]] <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = var, y = count, fill = var)) + 
          ggplot2::geom_bar(stat = "identity", show.legend = F) + 
          ggplot2::geom_text(mapping = ggplot2::aes(label = count), vjust = -0.2, size = 5) + 
          ggplot2::labs(x = var.to.bar, y = "count", title = t) + 
          ggplot2::theme_classic() + 
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) + 
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
          ggplot2::theme(text = ggplot2::element_text(size = f))
      }
      grDevices::pdf(file = base::file.path(save.dir, 
                                            base::sprintf(fmt = "barplots1_%spercluster_%s_res%s%s.pdf", 
                                                          var.to.bar, assay.clustered, res.to.use, fig.suf)), 
                     width = w2, height = (7 * num.clusters), onefile = T, paper = "special")
      base::print(x = cowplot::plot_grid(plotlist = bar.plots, nrow = num.clusters, ncol = 1))
      grDevices::dev.off()
      
      ## for each sample, number and fraction of cells from each cell type
      bar.plots <- base::list()
      for (curr.bar.var in all.bar.vars) {
        cc <- cell.counts[,curr.bar.var]
        df <- base::data.frame(cluster = base::names(x = cc), count = cc)
        t <- base::sprintf(fmt = "%s: %s cells", curr.bar.var, base::sum(cc))
        bar.plots[[t]] <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = cluster, y = count, fill = cluster)) + 
          ggplot2::geom_bar(stat = "identity", show.legend = F) + 
          ggplot2::geom_text(mapping = ggplot2::aes(label = count), vjust = -0.2, size = 5) + 
          ggplot2::labs(x = "cluster", y = "count", title = t) + 
          ggplot2::theme_classic() + 
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) + 
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
          ggplot2::theme(text = ggplot2::element_text(size = f))
      }
      w <- num.clusters
      if (7 > w) {
        w <- 7
      }
      grDevices::pdf(file = base::file.path(save.dir, 
                                            base::sprintf(fmt = "barplots1_clusterper%s_%s_res%s%s.pdf", 
                                                          var.to.bar, assay.clustered, res.to.use, fig.suf)), 
                     width = w, height = (7 * num.bar.vars), onefile = T, paper = "special")
      base::print(x = cowplot::plot_grid(plotlist = bar.plots, nrow = num.bar.vars, ncol = 1))
      grDevices::dev.off()
    }
    
    if (do.dot) {
      ## The values in DotPlot are extracted from the @data slot, averaged, and then passed to scale.
      ## These are then Min-Maxed based on the col.min and col.max parameter values.
      ## The DotPlot shows scaled values (which can be both positive and negative), 
      ## so that we can visualize highly and lowly expressed genes on the same color scale.
      dp <- Seurat::DotPlot(object = seurat.obj, assay = assay.to.plot, features = markers.to.plot, dot.scale = d, 
                            cluster.idents = F, scale = T) + 
        Seurat::RotatedAxis() + 
        ggplot2::labs(title = plot.title) + 
        viridis::scale_color_viridis(option = "D", direction = -1) + 
        Seurat::FontSize(x.text = f, y.text = f, x.title = f, y.title = f, main = f) + 
        ggplot2::theme(axis.line = ggplot2::element_blank(), 
                       panel.grid = ggplot2::element_line(colour = "gray", size = 0.5, linetype = "dashed"), 
                       panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1))
      h <- 0.7 * num.clusters
      if (7 > h) {
        h <- 7
      }
      grDevices::pdf(file = base::file.path(save.dir, 
                                            base::sprintf(fmt = "dotplot1_%s_res%s%s.pdf", 
                                                          assay.to.plot, res.to.use, fig.suf)), 
                     width = w3, height = h, onefile = T, paper = "special")
      base::print(x = dp)
      grDevices::dev.off()
    }
  }
}

PlotGenes <- function(seurat.obj, act.assay, genes.to.plot, save.fig.dir, dims.to.plot = base::c(1, 2), 
                      red.to.plot = "umap", var.to.split = NULL, fn.suf = "", col.n = 2, plot.label = T) {
  ### feature plots
  Seurat::DefaultAssay(object = seurat.obj) <- act.assay
  c <- base::c("grey", "red")
  l <- 5
  nc <- col.n
  nr <- base::ceiling(x = (base::length(x = genes.to.plot) / nc))
  if (!base::is.null(x = var.to.split)) {
    nc <- base::length(x = base::unique(x = seurat.obj@meta.data[[var.to.split]]))
    nr <- base::length(x = genes.to.plot)
  }
  w <- 7 * nc
  h <- 7 * nr
  
  fp <- Seurat::FeaturePlot(object = seurat.obj, features = genes.to.plot, dims = dims.to.plot, cols = c, order = T, 
                            reduction = red.to.plot, split.by = var.to.split, slot = "data", label = plot.label, 
                            label.size = l, repel = T, ncol = nc, combine = T)
  grDevices::pdf(file = base::file.path(save.fig.dir, base::sprintf(fmt = "%s_genes%s.pdf", red.to.plot, fn.suf)), 
                 width = w, height = h, onefile = T, paper = "special")
  base::print(x = fp)
  grDevices::dev.off()
}

PlotIntegration <- function(seurat.obj, act.assay, features.to.plot, save.dir, fn.suf = "") {
  ## violin plot
  w <- 0.5 * base::length(x = base::unique(x = seurat.obj$sample))
  if (7 > w) {
    w <- 7
  }
  h <- 7
  p <- 0
  grDevices::pdf(file = base::file.path(save.dir, base::sprintf("violin_embeddings%s.pdf", fn.suf)), 
                 width = w, height = h, onefile = T, paper = "special")
  vp <- Seurat::VlnPlot(object = seurat.obj, features = features.to.plot, pt.size = p, assay = act.assay, 
                        group.by = "sample", slot = "data") + Seurat::NoLegend()
  base::print(x = vp)
  grDevices::dev.off()
}

PlotMonocle2 <- function(cds, feature.to.plot, save.fig.dir, feature.name, fn.suf, plot.together = T, plot.split = T) {
  if (plot.together) {
    ctp <- monocle::plot_cell_trajectory(cds = cds, color_by = feature.to.plot, cell_size = 0.3) + 
      ggplot2::theme(legend.position = "right")
    grDevices::pdf(file = base::file.path(save.fig.dir, base::sprintf(fmt = "trajectory_monocle2_%s%s.pdf", feature.name, fn.suf)), 
                   width = 9, height = 7, onefile = T, paper = "special")
    base::print(x = ctp)
    grDevices::dev.off()
  }
  
  if (plot.split) {
    feature.n <- base::length(x = base::unique(x = Biobase::pData(object = cds)[[feature.to.plot]]))
    base::cat("\n")
    base::print(x = feature.n)
    nc <- base::min(10, feature.n)
    base::print(x = nc)
    nr <- base::ceiling(x = feature.n / nc)
    base::print(x = nr)
    ctp <- monocle::plot_cell_trajectory(cds = cds, color_by = feature.to.plot, cell_size = 0.4) + 
      ggplot2::facet_wrap(facets = base::c(feature.to.plot), ncol = nc, dir = "h") + 
      ggplot2::theme(legend.position = "none")
    grDevices::pdf(file = base::file.path(save.fig.dir, base::sprintf(fmt = "trajectory_monocle2_split%s%s.pdf", feature.name, fn.suf)), 
                   width = 5 * nc, height = 5 * nr, onefile = T, paper = "special")
    base::print(x = ctp)
    grDevices::dev.off()
  }
}

PlotNumberOfCells <- function(seurat.obj, var.to.count, save.fig.dir, x.label = NULL, fn.suf = "") {
  metadata <- seurat.obj@meta.data
  f1 <- 5
  f2 <- 18
  p <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes_string(x = var.to.count, fill = var.to.count)) + 
    ggplot2::geom_bar(show.legend = F) + 
    ggplot2::geom_text(mapping = ggplot2::aes(label = ggplot2::after_stat(x = count)), stat = "count", 
                       vjust = -0.2, size = f1) + 
    ggplot2::ggtitle(label = base::sprintf(fmt = "%s cells in total", base::length(x = metadata[[var.to.count]]))) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), 
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"), 
      text = ggplot2::element_text(size = f2)
    )
  if (base::is.null(x = x.label)) {
    p <- p + ggplot2::theme(axis.title.x = ggplot2::element_blank())
  } else {
    p <- p + ggplot2::xlab(label = x.label)
  }
  
  w <- 0.5 * base::length(x = base::unique(x = metadata[[var.to.count]]))
  if (w < 7) {
    w <- 7
  }
  h <- 7
  grDevices::pdf(file = base::file.path(save.fig.dir, base::sprintf(fmt = "number_of_cells%s.pdf", fn.suf)), 
                 width = w, height = h, onefile = T, paper = "special")
  base::print(x = p)
  grDevices::dev.off()
}

PlotPredictionScore <- function(seurat.obj, var.predicted, save.fig.dir) {
  metadata <- seurat.obj@meta.data
  f <- 18
  p <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes_string(x = var.predicted, 
                                                                      y = base::sprintf(fmt = "%s.score", var.predicted), 
                                                                      color = var.predicted)) + 
    ggplot2::geom_boxplot(show.legend = F) + 
    ggplot2::geom_hline(yintercept = 0.5) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
    ggplot2::theme(text = ggplot2::element_text(size = f))
  w <- 0.5 * base::length(x = base::unique(x = metadata[[var.predicted]]))
  if (14 > w) {
    w <- 14
  }
  h <- 7
  grDevices::pdf(file = base::file.path(save.fig.dir, base::sprintf(fmt = "%s_score.pdf", var.predicted)), 
                 width = w, height = h, onefile = T, paper = "special")
  base::print(x = p)
  grDevices::dev.off()
}

PlotPrincipalComponents <- function(seurat.obj, save.dir, dim.to.plot, red.to.use, fn.suf = "", do.elbow = T) {
  Seurat::Idents(object = seurat.obj) <- "sample"
  
  ### Dimensional reduction plot
  plots <- base::list()
  for (d in 1:(dim.to.plot - 1)) {
    plots[[d]] <- Seurat::DimPlot(object = seurat.obj, dims = c(d, d + 1), reduction = red.to.use) + Seurat::NoLegend()
  }
  
  nc <- 10
  nr <- base::ceiling(x = (base::length(x = plots) / nc))
  w <- 5 * nc
  h <- 5 * nr
  grDevices::jpeg(file = base::file.path(save.dir, base::sprintf(fmt = "PC%s.jpg", fn.suf)), width = w, height = h, 
                  units = "in", pointsize = 12, quality = 100, res = 300)
  base::print(x = cowplot::plot_grid(plotlist = plots, nrow = nr, ncol = nc))
  grDevices::dev.off()
  
  if (do.elbow) {
    ### Elbow plots
    grDevices::pdf(file = base::file.path(save.dir, base::sprintf(fmt = "elbow1%s.pdf", fn.suf)), width = 14, height = 7, 
                   onefile = T, paper = "special")
    base::print(x = Seurat::ElbowPlot(object = seurat.obj, ndims = dim.to.plot, reduction = red.to.use))
    grDevices::dev.off()
    
    ## Determine percent of variation associated with each PC
    pct <- seurat.obj[[red.to.use]]@stdev / base::sum(seurat.obj[[red.to.use]]@stdev) * 100
    ## Calculate cumulative percents for each PC
    cumu <- base::cumsum(pct)
    ## Determine which PC exhibits cumulative percent greater than 90% 
    ## and % variation associated with the PC as less than 5
    co1 <- base::which((cumu > 90) & (pct < 5))[1]
    ## Determine the difference between variation of PC and subsequent PC
    ## last point where change of % of variation is more than 0.1%.
    co2 <- base::sort(base::which((pct[1:(base::length(pct) - 1)] - pct[2:base::length(pct)]) > 0.1), 
                      decreasing = T)[1] + 1 
    # Minimum of the two calculations
    pcs <- base::min(co1, co2)
    
    base::sink(file = base::file.path(save.dir, base::sprintf(fmt = "PC%s.txt", fn.suf)), append = F, split = F)
    base::cat(base::sprintf(fmt = "%s %s\n", co1, co2))
    base::cat(base::sprintf(fmt = "%s\n", pcs))
    base::sink()
    
    ## Create a dataframe with values
    plot_df <- base::data.frame(pct = pct, 
                                cumu = cumu, 
                                rank = 1:base::length(pct)) 
    ## Elbow plot to visualize
    p <- ggplot2::ggplot(data = plot_df, mapping = ggplot2::aes(x = cumu, y = pct, label = rank, color = (rank > pcs))) + 
      ggplot2::geom_text() + 
      ggplot2::geom_vline(xintercept = 90, color = "grey") + 
      ggplot2::geom_hline(yintercept = base::min(pct[pct > 5]), color = "grey") + 
      ggplot2::theme_bw()
    grDevices::pdf(file = base::file.path(save.dir, base::sprintf(fmt = "elbow2%s.pdf", fn.suf)), width = 14, height = 7, 
                   onefile = T, paper = "special")
    base::print(x = p)
    grDevices::dev.off()
    
    base::return(pcs)
  }
}

PlotQuality <- function(seurat.obj, save.fig.dir, red.to.plot = "umap", feature.to.split = NULL, fn.suf = "", 
                        x.label = "cluster") {
  features.to.plot <- base::c("nCount_ATAC", "pct_reads_in_peaks", "blacklist_ratio", "nucleosome_signal", "TSS.enrichment")
  c <- base::c("grey", "red")
  l <- 5
  nc <- base::length(x = features.to.plot)
  nr <- 1
  if (!base::is.null(x = feature.to.split)) {
    nc <- base::length(x = base::unique(x = seurat.obj@meta.data[[feature.to.split]]))
    nr <- base::length(x = features.to.plot)
  }
  w <- 7 * nc
  h <- 7 * nr
  fp <- Seurat::FeaturePlot(object = seurat.obj, features = features.to.plot, cols = c, order = T, reduction = red.to.plot, 
                            split.by = feature.to.split, slot = "data", label = T, label.size = l, repel = T, ncol = nc, combine = T)
  grDevices::jpeg(file = base::file.path(save.fig.dir, base::sprintf(fmt = "%s_quality%s.jpg", red.to.plot, fn.suf)), 
                  width = w, height = h, units = "in", pointsize = 12, quality = 100, res = 300)
  base::print(x = fp)
  grDevices::dev.off()
  
  seurat.obj$act.clusters <- base::as.character(x = Seurat::Idents(object = seurat.obj))
  metadata <- seurat.obj@meta.data
  f <- 18
  plots <- base::list()
  
  ## Visualize the total number of fragments in peaks.
  plots$prf <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes(x = act.clusters, y = nCount_ATAC, color = act.clusters)) + 
    ggplot2::geom_boxplot(show.legend = F) + 
    ggplot2::labs(x = x.label) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), 
      text = ggplot2::element_text(size = f)
    )
  
  ## Visualize the fraction of fragments in peaks.
  plots$prp <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes(x = act.clusters, y = pct_reads_in_peaks, color = act.clusters)) + 
    ggplot2::geom_boxplot(show.legend = F) + 
    ggplot2::labs(x = x.label) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), 
      text = ggplot2::element_text(size = f)
    )
  
  ## Visualize the ratio of reads in genomic blacklist regions.
  plots$br <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes(x = act.clusters, y = blacklist_ratio, color = act.clusters)) + 
    ggplot2::geom_boxplot(show.legend = F) + 
    ggplot2::labs(x = x.label) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), 
      text = ggplot2::element_text(size = f)
    )
  
  ## Visualize the nucleosome signal score per cell.
  plots$ns <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes(x = act.clusters, y = nucleosome_signal, color = act.clusters)) + 
    ggplot2::geom_boxplot(show.legend = F) + 
    ggplot2::labs(x = x.label) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), 
      text = ggplot2::element_text(size = f)
    )
  
  ## Visualize the TSS enrichment score per cell.
  plots$tsse <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes(x = act.clusters, y = TSS.enrichment, color = act.clusters)) + 
    ggplot2::geom_boxplot(show.legend = F) + 
    ggplot2::labs(x = x.label) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), 
      text = ggplot2::element_text(size = f)
    )
  
  nc <- 1
  nr <- base::length(x = plots)
  w <- 0.5 * base::length(x = base::levels(x = seurat.obj))
  if (14 > w) {
    w <- 14
  }
  h <- 7 * nr
  grDevices::pdf(file = base::file.path(save.fig.dir, base::sprintf(fmt = "boxplot_quality%s.pdf", fn.suf)), 
                 width = w, height = h, onefile = T, paper = "special")
  base::print(x = cowplot::plot_grid(plotlist = plots, nrow = nr, ncol = nc))
  grDevices::dev.off()
}

PlotQuality.Rna <- function(seurat.obj, save.fig.dir, features.to.add = NULL, red.to.plot = "umap", 
                            feature.to.split = NULL, fn.suf = "", x.label = "cluster") {
  features.to.plot <- base::c("nUMI", "nGene", "mitoRatio", features.to.add)
  c <- base::c("grey", "red")
  l <- 5
  nc <- base::length(x = features.to.plot)
  nr <- 1
  if (!base::is.null(x = feature.to.split)) {
    nc <- base::length(x = base::unique(x = seurat.obj@meta.data[[feature.to.split]]))
    nr <- base::length(x = features.to.plot)
  }
  w <- 7 * nc
  h <- 7 * nr
  
  fp <- Seurat::FeaturePlot(object = seurat.obj, features = features.to.plot, cols = c, order = T, 
                            reduction = red.to.plot, split.by = feature.to.split, slot = "data", label = T, 
                            label.size = l, repel = T, ncol = nc, combine = T)
  grDevices::jpeg(file = base::file.path(save.fig.dir, base::sprintf(fmt = "%s_quality%s.jpg", red.to.plot, fn.suf)), 
                  width = w, height = h, units = "in", pointsize = 12, quality = 100, res = 300)
  base::print(x = fp)
  grDevices::dev.off()
  
  seurat.obj$act.clusters <- base::as.character(x = Seurat::Idents(object = seurat.obj))
  metadata <- seurat.obj@meta.data
  f <- 18
  plots <- base::list()
  
  ## Visualize the number of UMIs/transcripts per cell
  plots$nUMI <- ggplot2::ggplot(data = metadata, 
                                mapping = ggplot2::aes(x = act.clusters, y = nUMI, color = act.clusters)) + 
    ggplot2::geom_boxplot(show.legend = F) + 
    ggplot2::labs(x = x.label) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), 
      text = ggplot2::element_text(size = f)
    )
  
  ## Visualize the number of genes detected per cell
  plots$nGene <- ggplot2::ggplot(data = metadata, 
                                 mapping = ggplot2::aes(x = act.clusters, y = nGene, color = act.clusters)) + 
    ggplot2::geom_boxplot(show.legend = F) + 
    ggplot2::labs(x = x.label) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), 
      text = ggplot2::element_text(size = f)
    )
  
  ## Visualize the percentage of mitochondrial genes detected per cell
  plots$mitoRatio <- ggplot2::ggplot(data = metadata, 
                                     mapping = ggplot2::aes(x = act.clusters, y = mitoRatio, color = act.clusters)) + 
    ggplot2::geom_boxplot(show.legend = F) + 
    ggplot2::labs(x = x.label) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), 
      text = ggplot2::element_text(size = f)
    )
  
  ## Visualizing the number of genes detected per transcript
  plots$nGenePerUMI <- ggplot2::ggplot(data = metadata, 
                                       mapping = ggplot2::aes(x = act.clusters, y = nGenePerUMI, color = act.clusters)) + 
    ggplot2::geom_boxplot(show.legend = F) + 
    ggplot2::labs(x = x.label) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), 
      text = ggplot2::element_text(size = f)
    )
  
  nc <- 1
  nr <- base::length(x = plots)
  w <- 0.5 * base::length(x = base::levels(x = seurat.obj))
  if (14 > w) {
    w <- 14
  }
  h <- 7 * nr
  grDevices::pdf(file = base::file.path(save.fig.dir, base::sprintf(fmt = "boxplot_quality%s.pdf", fn.suf)), 
                 width = w, height = h, onefile = T, paper = "special")
  base::print(x = cowplot::plot_grid(plotlist = plots, nrow = nr, ncol = nc))
  grDevices::dev.off()
}

PlotQualityMetrics <- function(seurat.obj, save.fig.dir, var.to.split = "sample", unsplit.box = T, fn.suf = "", 
                               all.metrics = T) {
  seurat.obj$all.cells <- "all cells"
  metadata <- seurat.obj@meta.data
  features.to.plot <- base::c("pct_reads_in_peaks", "nCount_ATAC")
  if (all.metrics) {
    features.to.plot <- base::c(features.to.plot, "TSS.enrichment", "blacklist_ratio", 
                                "nucleosome_signal")
  }
  f <- 18
  p <- 0.4
  plots <- base::list()
  
  ## Visualize the number of cells per sample.
  plots$noc <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes_string(x = var.to.split, fill = var.to.split)) + 
    ggplot2::geom_bar(show.legend = F) + 
    ggplot2::geom_text(mapping = ggplot2::aes(label = ggplot2::after_stat(x = count), vjust = -0.2), stat = "count") + 
    ggplot2::theme_classic() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) + 
    ggplot2::ggtitle(label = base::sprintf(fmt = "%s Cellular Barcodes", base::length(x = metadata[[var.to.split]]))) + 
    ggplot2::theme(text = ggplot2::element_text(size = f))
  
  ## Visualize the total number of fragments in peaks.
  plots$prf <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes_string(x = var.to.split, 
                                                                                  y = "nCount_ATAC", 
                                                                                  color = var.to.split)) + 
    ggplot2::geom_boxplot(show.legend = F) + 
    ggplot2::geom_hline(yintercept = 3000) + 
    ggplot2::geom_hline(yintercept = 20000) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
    ggplot2::theme(text = ggplot2::element_text(size = f))
  
  ## Visualize the fraction of fragments in peaks.
  plots$prp <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes_string(x = var.to.split, 
                                                                                  y = "pct_reads_in_peaks", 
                                                                                  color = var.to.split)) + 
    ggplot2::geom_boxplot(show.legend = F) + 
    ggplot2::geom_hline(yintercept = 15) + 
    ggplot2::geom_hline(yintercept = 40) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
    ggplot2::theme(text = ggplot2::element_text(size = f))
  
  ## Visualize the ratio of reads in mitochondrial DNA regions.
  # plots$mr <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes_string(x = var.to.split, 
  #                                                                            y = "mito_ratio", 
  #                                                                            color = var.to.split)) + 
  #   ggplot2::geom_boxplot(show.legend = F) + 
  #   ggplot2::geom_hline(yintercept = 0.5) + 
  #   ggplot2::geom_hline(yintercept = 0.2) + 
  #   ggplot2::theme_classic() + 
  #   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
  #   ggplot2::theme(text = ggplot2::element_text(size = f))
  
  if (all.metrics) {
    ## Visualize the ratio of reads in genomic blacklist regions.
    plots$br <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes_string(x = var.to.split, 
                                                                               y = "blacklist_ratio", 
                                                                               color = var.to.split)) + 
      ggplot2::geom_boxplot(show.legend = F) + 
      ggplot2::geom_hline(yintercept = 0.025) + 
      ggplot2::geom_hline(yintercept = 0.05) + 
      ggplot2::theme_classic() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
      ggplot2::theme(text = ggplot2::element_text(size = f))
    
    ## Visualize the nucleosome signal score per cell.
    plots$ns <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes_string(x = var.to.split, 
                                                                               y = "nucleosome_signal", 
                                                                               color = var.to.split)) + 
      ggplot2::geom_boxplot(show.legend = F) + 
      ggplot2::geom_hline(yintercept = 4) + 
      ggplot2::theme_classic() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
      ggplot2::theme(text = ggplot2::element_text(size = f))
    
    ## Visualize the frequency that fragments of different lengths are present for different groups of cells.
    seurat.obj$nucleosome_group <- base::ifelse(test = (seurat.obj$nucleosome_signal > 4), 
                                                yes = "NS > 4", no = "NS < 4")
    plots$fh <- Signac::FragmentHistogram(object = seurat.obj, region = "chr1-1-10000000", 
                                          group.by = "nucleosome_group")
    
    ## Visualize the TSS enrichment score per cell.
    plots$tsse <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes_string(x = var.to.split, 
                                                                                 y = "TSS.enrichment", 
                                                                                 color = var.to.split)) + 
      ggplot2::geom_boxplot(show.legend = F) + 
      ggplot2::geom_hline(yintercept = 2) + 
      ggplot2::theme_classic() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
      ggplot2::theme(text = ggplot2::element_text(size = f))
    
    ## Visualize the normalized TSS enrichment score at each position relative to the TSS.
    seurat.obj$high.tss <- base::ifelse(test = (seurat.obj$TSS.enrichment > 2), yes = "High", no = "Low")
    plots$tssp <- Signac::TSSPlot(object = seurat.obj, group.by = "high.tss") + Seurat::NoLegend()
  }
  
  if (unsplit.box) {
    ## Visualize the total number of fragments in peaks.
    plots$unsplit.prf <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes(x = all.cells, 
                                                                                 y = nCount_ATAC, 
                                                                                 color = all.cells)) + 
      ggplot2::geom_boxplot(show.legend = F) + 
      ggplot2::geom_hline(yintercept = 3000) + 
      ggplot2::geom_hline(yintercept = 20000) + 
      ggplot2::theme_classic() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
      ggplot2::theme(text = ggplot2::element_text(size = f))
    
    ## Visualize the fraction of fragments in peaks.
    plots$unsplit.prp <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes(x = all.cells, 
                                                                                 y = pct_reads_in_peaks, 
                                                                                 color = all.cells)) + 
      ggplot2::geom_boxplot(show.legend = F) + 
      ggplot2::geom_hline(yintercept = 15) + 
      ggplot2::geom_hline(yintercept = 40) + 
      ggplot2::theme_classic() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
      ggplot2::theme(text = ggplot2::element_text(size = f))
    
    ## Visualize the ratio of reads in mitochondrial DNA regions.
    # plots$unsplit.mr <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes(x = all.cells, 
    #                                                                             y = mito_ratio, 
    #                                                                             color = all.cells)) + 
    #   ggplot2::geom_boxplot(show.legend = F) + 
    #   ggplot2::geom_hline(yintercept = 0.5) + 
    #   ggplot2::geom_hline(yintercept = 0.2) + 
    #   ggplot2::theme_classic() + 
    #   ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
    #   ggplot2::theme(text = ggplot2::element_text(size = f))
    
    if (all.metrics) {
      ## Visualize the ratio of reads in genomic blacklist regions.
      plots$unsplit.br <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes(x = all.cells, 
                                                                                  y = blacklist_ratio, 
                                                                                  color = all.cells)) + 
        ggplot2::geom_boxplot(show.legend = F) + 
        ggplot2::geom_hline(yintercept = 0.025) + 
        ggplot2::geom_hline(yintercept = 0.05) + 
        ggplot2::theme_classic() + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
        ggplot2::theme(text = ggplot2::element_text(size = f))
      
      ## Visualize the nucleosome signal score per cell.
      plots$unsplit.ns <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes(x = all.cells, 
                                                                                  y = nucleosome_signal, 
                                                                                  color = all.cells)) + 
        ggplot2::geom_boxplot(show.legend = F) + 
        ggplot2::geom_hline(yintercept = 4) + 
        ggplot2::theme_classic() + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
        ggplot2::theme(text = ggplot2::element_text(size = f))
      
      ## Visualize the TSS enrichment score per cell.
      plots$unsplit.tsse <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes(x = all.cells, 
                                                                                    y = TSS.enrichment, 
                                                                                    color = all.cells)) + 
        ggplot2::geom_boxplot(show.legend = F) + 
        ggplot2::geom_hline(yintercept = 2) + 
        ggplot2::theme_classic() + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
        ggplot2::theme(text = ggplot2::element_text(size = f))
    }
  }
  
  plots$vp <- Seurat::VlnPlot(object = seurat.obj, features = features.to.plot, pt.size = 0, 
                              group.by = "all.cells", ncol = base::length(x = features.to.plot))
  vp <- Seurat::VlnPlot(object = seurat.obj, features = features.to.plot, pt.size = p, 
                        group.by = "all.cells", ncol = base::length(x = features.to.plot))
  
  nc <- 1
  nr <- base::length(x = plots)
  w <- 0.5 * base::length(x = base::unique(x = metadata[[var.to.split]]))
  if (14 > w) {
    w <- 14
  }
  h <- 7 * nr
  grDevices::pdf(file = base::file.path(save.fig.dir, base::sprintf("QC%s.pdf", fn.suf)), width = w, height = h,
                 onefile = T, paper = "special")
  base::print(x = cowplot::plot_grid(plotlist = plots, nrow = nr, ncol = nc))
  grDevices::dev.off()
  grDevices::jpeg(file = base::file.path(save.fig.dir, base::sprintf("QC%s.jpg", fn.suf)), width = w, height = 7,
                  units = "in", pointsize = 12, quality = 100, res = 300)
  base::print(x = vp)
  grDevices::dev.off()
}

PlotQualityMetrics.Rna <- function(seurat.obj, save.fig.dir, var.to.split = "sample", features.to.add = NULL, do.correlation = T, 
                                   do.density = T, do.box.t = T, fn.suf = "", low.umi = 500, high.umi = 1000, low.gene = 250, 
                                   high.gene = 2500, low.mito = 0.2, high.mito = 0.5, low.complex = 0.25) {
  seurat.obj$all.cells <- "all cells"
  metadata <- seurat.obj@meta.data
  split.vars <- base::sort(x = base::unique(x = metadata[[var.to.split]]), decreasing = F)
  features.to.plot <- base::c("nUMI", "nGene", "mitoRatio", features.to.add)
  f <- 18
  p <- 0.4
  plots <- base::list()
  
  ## Visualize the cell count per sample
  plots$nCell <- ggplot2::ggplot(data = metadata, mapping = ggplot2::aes_string(x = var.to.split, fill = var.to.split)) + 
    ggplot2::geom_bar(show.legend = F) + 
    ggplot2::geom_text(mapping = ggplot2::aes(label = ggplot2::after_stat(x = count), vjust = -0.2), stat = "count") + 
    ggplot2::theme_classic() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) + 
    ggplot2::ggtitle(label = base::sprintf(fmt = "%s Cellular Barcodes", base::length(x = metadata[[var.to.split]]))) + 
    ggplot2::theme(text = ggplot2::element_text(size = f))
  
  if (do.correlation) {
    ## Visualize the correlation between genes detected and UMI counts 
    ## and determine whether strong presence of cells with low numbers of genes/UMIs
    plots$correlation <- ggplot2::ggplot(data = metadata, 
                                         mapping = ggplot2::aes(x = nUMI, y = nGene, color = mitoRatio)) + 
      ggplot2::geom_point() + 
      ggplot2::scale_colour_gradient(low = "gray90", high = "black") + 
      ggplot2::stat_smooth(method = lm) + 
      ggplot2::scale_x_log10() + 
      ggplot2::scale_y_log10() + 
      ggplot2::theme_classic() + 
      ggplot2::geom_vline(xintercept = low.umi) + 
      ggplot2::geom_vline(xintercept = high.umi) + 
      ggplot2::geom_hline(yintercept = low.gene) + 
      ggplot2::geom_hline(yintercept = high.gene) + 
      ggplot2::facet_wrap(stats::as.formula(object = base::paste("~", var.to.split, sep = ""))) + 
      ggplot2::theme(axis.title = ggplot2::element_text(size = f))
  }
  
  if (do.density) {
    ## Visualize the distribution of UMI count per cell
    plots$nUMI.density <- ggplot2::ggplot(data = metadata, 
                                          mapping = ggplot2::aes_string(color = var.to.split, x = "nUMI", 
                                                                        fill = var.to.split)) + 
      ggplot2::geom_density(alpha = 0.2) + 
      ggplot2::scale_x_log10() + 
      ggplot2::theme_classic() + 
      ggplot2::labs(y = "cell density") + 
      ggplot2::geom_vline(xintercept = low.umi) + 
      ggplot2::geom_vline(xintercept = high.umi) + 
      ggplot2::theme(text = ggplot2::element_text(size = f))
    
    ## Visualize the distribution of genes detected per cell
    plots$nGene.density <- ggplot2::ggplot(data = metadata, 
                                           mapping = ggplot2::aes_string(color = var.to.split, x = "nGene", 
                                                                         fill = var.to.split)) + 
      ggplot2::geom_density(alpha = 0.2) + 
      ggplot2::theme_classic() + 
      ggplot2::scale_x_log10() + 
      ggplot2::labs(y = "cell density") + 
      ggplot2::geom_vline(xintercept = low.gene) + 
      ggplot2::geom_vline(xintercept = high.gene) + 
      ggplot2::theme(text = ggplot2::element_text(size = f))
    
    ## Visualize the distribution of mitochondrial gene expression detected per cell
    plots$mitoRatio.density <- ggplot2::ggplot(data = metadata, 
                                               mapping = ggplot2::aes_string(color = var.to.split, x = "mitoRatio", 
                                                                             fill = var.to.split)) + 
      ggplot2::geom_density(alpha = 0.2) + 
      ggplot2::scale_x_log10() + 
      ggplot2::theme_classic() + 
      ggplot2::labs(y = "cell density") + 
      ggplot2::geom_vline(xintercept = low.mito) + 
      ggplot2::geom_vline(xintercept = high.mito) + 
      ggplot2::theme(text = ggplot2::element_text(size = f))
    
    ## Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
    plots$nGenePerUMI.density <- ggplot2::ggplot(data = metadata, 
                                                 mapping = ggplot2::aes_string(color = var.to.split, x = "nGenePerUMI", 
                                                                               fill = var.to.split)) +
      ggplot2::geom_density(alpha = 0.2) +
      ggplot2::theme_classic() +
      ggplot2::labs(y = "cell density") +
      ggplot2::geom_vline(xintercept = low.complex) +
      ggplot2::theme(text = ggplot2::element_text(size = f))
  }
  
  ## Visualize the distribution of UMI count per cell
  plots$nUMI.box <- ggplot2::ggplot(data = metadata, 
                                    mapping = ggplot2::aes_string(x = var.to.split, y = "nUMI", color = var.to.split)) + 
    ggplot2::geom_boxplot(show.legend = F) + 
    ggplot2::geom_hline(yintercept = 500) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
    ggplot2::theme(text = ggplot2::element_text(size = f))
  
  ## Visualize the distribution of genes detected per cell
  plots$nGene.box <- ggplot2::ggplot(data = metadata, 
                                     mapping = ggplot2::aes_string(x = var.to.split, y = "nGene", color = var.to.split)) + 
    ggplot2::geom_boxplot(show.legend = F) + 
    ggplot2::geom_hline(yintercept = 1500) + 
    ggplot2::geom_hline(yintercept = 1000) + 
    ggplot2::geom_hline(yintercept = 500) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
    ggplot2::theme(text = ggplot2::element_text(size = f))
  
  ## Visualize the distribution of mitochondrial gene expression detected per cell
  plots$mitoRatio.box <- ggplot2::ggplot(data = metadata, 
                                         mapping = ggplot2::aes_string(x = var.to.split, y = "mitoRatio", 
                                                                       color = var.to.split)) + 
    ggplot2::geom_boxplot(show.legend = F) + 
    ggplot2::geom_hline(yintercept = 0.05) + 
    ggplot2::geom_hline(yintercept = 0.01) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
    ggplot2::theme(text = ggplot2::element_text(size = f))
  
  ## Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
  plots$nGenePerUMI.box <- ggplot2::ggplot(data = metadata, 
                                           mapping = ggplot2::aes_string(x = var.to.split, y = "nGenePerUMI", 
                                                                         color = var.to.split)) + 
    ggplot2::geom_boxplot(show.legend = F) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
    ggplot2::theme(text = ggplot2::element_text(size = f))
  
  if (do.box.t) {
    ## Visualize the distribution of UMI count per cell
    plots$nUMI.box.t <- ggplot2::ggplot(data = metadata, 
                                        mapping = ggplot2::aes(x = all.cells, y = nUMI, color = all.cells)) + 
      ggplot2::geom_boxplot(show.legend = F) + 
      ggplot2::geom_hline(yintercept = 500) + 
      ggplot2::theme_classic() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
      ggplot2::theme(text = ggplot2::element_text(size = f))
    
    ## Visualize the distribution of genes detected per cell
    plots$nGene.box.t <- ggplot2::ggplot(data = metadata, 
                                         mapping = ggplot2::aes(x = all.cells, y = nGene, color = all.cells)) + 
      ggplot2::geom_boxplot(show.legend = F) + 
      ggplot2::geom_hline(yintercept = 1500) + 
      ggplot2::geom_hline(yintercept = 1000) + 
      ggplot2::geom_hline(yintercept = 500) + 
      ggplot2::theme_classic() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
      ggplot2::theme(text = ggplot2::element_text(size = f))
    
    ## Visualize the distribution of mitochondrial gene expression detected per cell
    plots$mitoRatio.box.t <- ggplot2::ggplot(data = metadata, 
                                             mapping = ggplot2::aes(x = all.cells, y = mitoRatio, color = all.cells)) + 
      ggplot2::geom_boxplot(show.legend = F) + 
      ggplot2::geom_hline(yintercept = 0.05) + 
      ggplot2::geom_hline(yintercept = 0.01) + 
      ggplot2::theme_classic() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
      ggplot2::theme(text = ggplot2::element_text(size = f))
    
    ## Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
    plots$nGenePerUMI.box.t <- ggplot2::ggplot(data = metadata, 
                                               mapping = ggplot2::aes(x = all.cells, y = nGenePerUMI, 
                                                                      color = all.cells)) + 
      ggplot2::geom_boxplot(show.legend = F) + 
      ggplot2::theme_classic() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) + 
      ggplot2::theme(text = ggplot2::element_text(size = f))
  }
  
  plots$nodeless.violin <- Seurat::VlnPlot(object = seurat.obj, features = features.to.plot, pt.size = 0, 
                                           group.by = "all.cells", ncol = base::length(x = features.to.plot))
  vp <- Seurat::VlnPlot(object = seurat.obj, features = features.to.plot, pt.size = p, 
                        group.by = "all.cells", ncol = base::length(x = features.to.plot))
  
  nc <- 1
  nr <- base::length(x = plots)
  w <- 0.5 * base::length(x = split.vars)
  if (14 > w) {
    w <- 14
  }
  h <- 7 * nr
  grDevices::pdf(file = base::file.path(save.fig.dir, base::sprintf("QC%s.pdf", fn.suf)), width = w, height = h,
                 onefile = T, paper = "special")
  base::print(x = cowplot::plot_grid(plotlist = plots, nrow = nr, ncol = nc))
  grDevices::dev.off()
  grDevices::jpeg(file = base::file.path(save.fig.dir, base::sprintf("QC%s.jpg", fn.suf)), width = w, height = 7,
                  units = "in", pointsize = 12, quality = 100, res = 300)
  base::print(x = vp)
  grDevices::dev.off()
}

QueryByDegree <- function(data, groups, colors.by.degree, mode, ...) {
  upset.data = ComplexUpset::upset_data(data = data, intersect = groups, mode = mode, ...)
  intersections = base::unique(x = (upset.data$plot_intersections_subset))
  base::lapply(
    X = intersections,
    FUN = function(x) {
      members = base::strsplit(x = x, split = "-", fixed = T)[[1]]
      if (!utils::hasName(x = colors.by.degree, name = base::as.character(x = base::length(x = members)))) {
        base::stop(base::sprintf(fmt = "missing specification of colors for degree %s.", base::length(x = members)))
      }
      args = base::c(
        base::list(intersect = members, only_components = "intersections_matrix"), 
        colors.by.degree[[base::as.character(x = base::length(x = members))]]
      )
      if (0 < upset.data$sizes[[mode]][x]) {
        args = base::c(
          base::list(intersect = members), 
          colors.by.degree[[base::as.character(x = base::length(x = members))]]
        )
      }
      base::do.call(what = ComplexUpset::upset_query, args = args)
    }
  )
}

Rat2MouseGenes <- function(genes.to.map, save.rst.dir = NULL, invert = F) {
  base::require(package = "biomaRt")
  
  rat <- biomaRt::useMart(biomart = "ensembl", dataset = "rnorvegicus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse <- biomaRt::useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  tmp <- NULL
  fn.pre <- NULL
  from.symbol <- NULL
  to.symbol <- NULL
  if (invert) {
    ## mouse to rat
    tmp <- biomaRt::getLDS(attributes = base::c("mgi_symbol", "chromosome_name", "start_position", "end_position"), 
                           filters = "mgi_symbol", 
                           values = genes.to.map, mart = mouse, 
                           attributesL = base::c("rgd_symbol", "chromosome_name", "start_position", "end_position"), 
                           martL = rat, uniqueRows = T)
    tmp$index <- base::match(x = tmp$MGI.symbol, table = genes.to.map)
    fn.pre <- "mouse_to_rat"
    from.symbol <- "MGI.symbol"
    to.symbol <- "RGD.symbol"
  } else {
    ## rat to mouse
    tmp <- biomaRt::getLDS(attributes = base::c("rgd_symbol", "chromosome_name", "start_position", "end_position"), 
                           filters = "rgd_symbol", 
                           values = genes.to.map, mart = rat, 
                           attributesL = base::c("mgi_symbol", "chromosome_name", "start_position", "end_position"), 
                           martL = mouse, uniqueRows = T)
    tmp$index <- base::match(x = tmp$RGD.symbol, table = genes.to.map)
    fn.pre <- "rat_to_mouse"
    from.symbol <- "RGD.symbol"
    to.symbol <- "MGI.symbol"
  }
  out.genes <- tmp[base::order(tmp$index, decreasing = F),]
  out.genes$index <- NULL
  
  if (!base::is.null(x = save.rst.dir)) {
    openxlsx::write.xlsx(x = out.genes, file = base::file.path(save.rst.dir, base::sprintf(fmt = "%s_genes.xlsx", fn.pre)), 
                         startCol = 1, startRow = 1, colNames = T, rowNames = F, overwrite = T)
  }
  
  out.genes <- out.genes[,base::c(from.symbol, to.symbol)]
  base::print(x = utils::head(x = out.genes, n = 3L))
  base::print(x = utils::tail(x = out.genes, n = 3L))
  
  from.genes <- NULL
  to.genes <- NULL
  for (i in 1:base::nrow(x = out.genes)) {
    curr.from.gene <- out.genes[i, from.symbol]
    curr.to.gene <- out.genes[i, to.symbol]
    if ((base::is.na(x = curr.from.gene)) || ("" == curr.from.gene) || (curr.from.gene %in% from.genes) || 
        (base::is.na(x = curr.to.gene)) || ("" == curr.to.gene) || (curr.to.gene %in% to.genes)) {
      next
    }
    from.genes <- base::c(from.genes, curr.from.gene)
    to.genes <- base::c(to.genes, curr.to.gene)
  }
  base::names(x = to.genes) <- from.genes
  
  base::return(to.genes)
}

ScaleRows <- function(x) { # function "scale_rows" from https://github.com/raivokolde/pheatmap/blob/master/R/pheatmap.r
  m = base::apply(X = x, MARGIN = 1, FUN = mean, na.rm = T)
  s = base::apply(X = x, MARGIN = 1, FUN = sd, na.rm = T)
  base::return((x - m) / s)
}

SelectClusters <- function(seurat.obj, clusters.to.select, save.rds.dir, fn.suf, invert = F, save.seurat.obj = T, 
                           select.from = NULL, fn.pre = "after_QC", do.profile = F, save.fig.dir = NULL, features.to.add = NULL) {
  base::cat("\n")
  base::print(x = base::table(Seurat::Idents(object = seurat.obj)))
  select.sobj <- subset(x = seurat.obj, idents = clusters.to.select, invert = invert)
  base::cat("\n")
  base::print(x = base::table(Seurat::Idents(object = select.sobj)))
  
  cells.to.select <- base::as.character(x = select.sobj$cellID)
  cells.clustered <- base::as.character(x = Seurat::Idents(object = select.sobj))
  base::names(x = cells.clustered) <- cells.to.select
  base::cat(base::sprintf(fmt = "\nNumber of select cells: %s\n", base::length(x = cells.clustered)))
  base::print(x = utils::head(x = cells.clustered, n = 3L))
  base::print(x = utils::tail(x = cells.clustered, n = 3L))
  select.sobj <- NULL
  base::gc()
  base::saveRDS(object = cells.clustered, file = base::file.path(save.rds.dir, base::sprintf(fmt = "select_cells%s.RDS", fn.suf)))
  
  if (save.seurat.obj) {
    if (base::is.null(x = select.from)) {
      select.from <- base::file.path(save.rds.dir, "after_QC.RDS")
    }
    sobj <- base::readRDS(file = select.from)
    select.sobj <- subset(x = sobj, subset = (cellID %in% cells.to.select))
    base::cat(base::sprintf(fmt = "\nNumber of select cells: %s\n", base::ncol(x = select.sobj)))
    base::print(x = utils::head(x = select.sobj[[]], n = 3L))
    base::print(x = utils::tail(x = select.sobj[[]], n = 3L))
    sobj <- NULL
    base::gc()
    base::saveRDS(object = select.sobj, file = base::file.path(save.rds.dir, base::sprintf(fmt = "%s%s.RDS", fn.pre, fn.suf)))
    
    if (do.profile) {
      PlotQualityMetrics.Rna(seurat.obj = select.sobj, save.fig.dir = save.fig.dir, var.to.split = "sample", 
                             features.to.add = features.to.add, do.correlation = F, do.density = F, do.box.t = T, 
                             fn.suf = base::sprintf(fmt = "_profile%s", fn.suf))
    }
  }
}

SortClusters <- function(seurat.obj, var.as.cluster = NULL) {
  if (!base::is.null(x = var.as.cluster)) {
    Seurat::Idents(object = seurat.obj) <- var.as.cluster
  }
  clusters <- base::levels(x = seurat.obj)
  clusters <- base::as.numeric(x = clusters)
  clusters <- base::sort(x = clusters, decreasing = F)
  clusters <- base::as.character(x = clusters)
  base::return(clusters)
}

TestPredictionScore <- function(seurat.obj, var.predicted, scores.to.test, red.to.plot, save.fig.dir, fn.suf) {
  var.as.score <- base::sprintf(fmt = "%s.score", var.predicted)
  base::cat("\n")
  base::print(x = var.as.score)
  Seurat::Idents(object = seurat.obj) <- var.predicted
  seurat.obj$act.score <- seurat.obj[[var.as.score]]
  base::cat("\n")
  base::print(x = seurat.obj)
  base::cat("\n")
  base::print(x = utils::head(x = seurat.obj[[]], n = 3L))
  base::print(x = utils::tail(x = seurat.obj[[]], n = 3L))
  
  l <- 5
  f <- 18
  h <- 7
  w1 <- 7
  for (score.thresh in scores.to.test) {
    base::cat(base::sprintf(fmt = "\nProcessing prediction score >= %s ...\n", score.thresh))
    
    sele.sobj <- subset(x = seurat.obj, subset = (act.score >= score.thresh))
    Seurat::Idents(object = sele.sobj) <- var.predicted
    base::cat("\n")
    base::print(x = sele.sobj)
    base::cat("\n")
    base::print(x = utils::head(x = sele.sobj[[]], n = 3L))
    base::print(x = utils::tail(x = sele.sobj[[]], n = 3L))
    base::cat("\n")
    base::print(x = base::min(sele.sobj$act.score))
    
    plot.title <- base::sprintf(fmt = "score >= %s: %s cells", score.thresh, base::ncol(x = sele.sobj))
    dp <- Seurat::DimPlot(object = sele.sobj, reduction = red.to.plot, label = T, label.size = l, repel = T) + 
      ggplot2::labs(title = plot.title) + 
      Seurat::FontSize(x.text = f, y.text = f, x.title = f, y.title = f, main = f) + 
      Seurat::NoLegend()
    grDevices::pdf(file = base::file.path(save.fig.dir, 
                                          base::sprintf(fmt = "%s_%s%s_scoreGE%s.pdf", 
                                                        red.to.plot, var.predicted, fn.suf, score.thresh)), 
                   width = w1, height = h, onefile = T, paper = "special")
    base::print(x = dp)
    grDevices::dev.off()
    
    sdp <- Seurat::DimPlot(object = sele.sobj, reduction = red.to.plot, label = T, label.size = l, repel = T, 
                           split.by = var.predicted) + 
      ggplot2::labs(title = plot.title) + 
      Seurat::FontSize(x.text = f, y.text = f, x.title = f, y.title = f, main = f) + 
      Seurat::NoLegend()
    w2 <- 7 * base::length(x = base::levels(x = sele.sobj))
    grDevices::pdf(file = base::file.path(save.fig.dir, 
                                          base::sprintf(fmt = "%s_split%s%s_scoreGE%s.pdf", 
                                                        red.to.plot, var.predicted, fn.suf, score.thresh)),
                   width = w2, height = h, onefile = T, paper = "special")
    base::print(x = sdp)
    grDevices::dev.off()
  }
}

Xlsx2List <- function(xlsx.file, startRow = 1, colNames = T, rowNames = T, rows = NULL, cols = NULL, sep.names = " ") {
  clusters <- openxlsx::getSheetNames(file = xlsx.file)
  cluster.markers <- base::list()
  for (cluster in clusters) {
    cluster.markers[[cluster]] <- openxlsx::read.xlsx(xlsxFile = xlsx.file, sheet = cluster, startRow = startRow, colNames = colNames, 
                                                      rowNames = rowNames, rows = rows, cols = cols, sep.names = sep.names, namedRegion = NULL)
  }
  base::return(cluster.markers)
}
