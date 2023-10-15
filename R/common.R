# common.R
base::cat("\nRunning \"common.R\" ...\n")

### Global settings ------
base::options(max.print = 99999, warn = 1)
base::suppressPackageStartupMessages({
  base::library(package = "plot1cell")
  base::library(package = "AnnotationHub")
  base::library(package = "Biobase")
  base::library(package = "BiocGenerics")
  base::library(package = "biomaRt")
  base::library(package = "CellChat")
  base::library(package = "circlize")
  base::library(package = "ComplexHeatmap")
  base::library(package = "ComplexUpset")
  base::library(package = "COSG")
  base::library(package = "data.table")
  base::library(package = "decoupleR")
  base::library(package = "DESeq2")
  base::library(package = "destiny")
  base::library(package = "devtools")
  base::library(package = "dplyr")
  base::library(package = "DoubletFinder")
  base::library(package = "EnsDb.Mmusculus.v79")
  base::library(package = "ensembldb")
  base::library(package = "factoextra")
  base::library(package = "future")
  base::library(package = "GenomeInfoDb")
  base::library(package = "GenomicRanges")
  base::library(package = "GGally")
  base::library(package = "ggplot2")
  base::library(package = "ggpubr")
  base::library(package = "ggvenn")
  base::library(package = "ggVennDiagram")
  base::library(package = "GOSemSim")
  base::library(package = "grDevices")
  base::library(package = "grid")
  base::library(package = "hdWGCNA")
  base::library(package = "igraph")
  base::library(package = "IRanges")
  base::library(package = "Matrix")
  base::library(package = "MetaNeighbor")
  base::library(package = "methods")
  base::library(package = "monocle")
  base::library(package = "nichenetr")
  base::library(package = "org.Hs.eg.db")
  base::library(package = "org.Mm.eg.db")
  base::library(package = "patchwork")
  base::library(package = "pheatmap")
  base::library(package = "reshape2")
  base::library(package = "reticulate")
  base::library(package = "rGREAT")
  base::library(package = "rjson")
  base::library(package = "rrvgo")
  base::library(package = "S4Vectors")
  base::library(package = "scales")
  base::library(package = "sceasy")
  base::library(package = "scITD")
  base::library(package = "scmap")
  base::library(package = "Seurat")
  base::library(package = "SeuratDisk")
  base::library(package = "SeuratWrappers")
  base::library(package = "Signac")
  base::library(package = "SoupX")
  base::library(package = "stats")
  base::library(package = "stringr")
  base::library(package = "tibble")
  base::library(package = "tidyr")
  base::library(package = "tidyverse")
  base::library(package = "utils")
  base::library(package = "uuid")
  base::library(package = "velocyto.R")
  base::library(package = "VGAM")
  base::library(package = "viridis")
  base::library(package = "WGCNA")
})
base::set.seed(seed = 1)

ensdb.mm.v79.ucsc.rds <- base::file.path(datasets.dir, "EnsDb.Mmusculus.v79.UCSC.RDS")
macs2.to.use <- "/home/jzhou88/installed/miniconda3/envs/r4/bin/macs2"

### Datasets ------
## 29 datasets
susztak.lab.rna <- base::c("FA1", "FA2", "UUO1", "UUO2",
                           "longIRI1d1", "longIRI1d2", "longIRI3d1", "longIRI3d2", "longIRI14d1", "longIRI14d2",
                           "shortIRI1d1", "shortIRI1d2", "shortIRI3d1", "shortIRI3d2", "shortIRI14d1", "shortIRI14d2",
                           "ESRRA",
                           "norm1", "norm2", "norm3", "norm4", "norm5", "MMctrl7",
                           "Notch3", "PGC1a-1", "PGC1a-2",
                           "G2CA-E", "G2CA-C",
                           "MMctrl8")
## 6 datasets
susztak.lab.human.dkd.snRNA <- base::c("HK2596", "HK2844", "HK2862", "HK2886", "HK2898", "HK2899")
## 7 datasets
gse151658 <- base::c("GSE151658-LPS0hr", "GSE151658-LPS1hr", "GSE151658-LPS4hr", "GSE151658-LPS16hr", "GSE151658-LPS27hr", 
                     "GSE151658-LPS36hr", "GSE151658-LPS48hr")
## 11 datasets
gse184652 <- base::c("GSE184652-dbdbPBSrep1", "GSE184652-dbdbPBSrep2", "GSE184652-dbdbPBSrep3", 
                     "GSE184652-dbdbPBSrep4", "GSE184652-dbdbPBSrep5", "GSE184652-dbdbPBSrep6", 
                     "GSE184652-dbmPBSrep1", "GSE184652-dbmPBSrep2", "GSE184652-dbmPBSrep3", 
                     "GSE184652-dbmPBSrep4", "GSE184652-dbmPBSrep5")
gse184652.ids <- base::c(
  "GSE184652-dbdbPBSrep1" = "N3024", 
  "GSE184652-dbdbPBSrep2" = "B3025", 
  "GSE184652-dbdbPBSrep3" = "A3026", 
  "GSE184652-dbdbPBSrep4" = "B3027", 
  "GSE184652-dbdbPBSrep5" = "C3028", 
  "GSE184652-dbdbPBSrep6" = "A3029", 
  "GSE184652-dbmPBSrep1" = "E3019", 
  "GSE184652-dbmPBSrep2" = "A3020", 
  "GSE184652-dbmPBSrep3" = "F3021", 
  "GSE184652-dbmPBSrep4" = "G3022", 
  "GSE184652-dbmPBSrep5" = "H3023"
)

## RNA disease datasets
rna.disease.samples <- base::c("FA1", "FA2", "Notch3", "Notch4", "PGC1a-1", "PGC1a-2", "UUO1", "UUO2", 
                               "longIRI1d1", "longIRI1d2", "longIRI3d1", "longIRI3d2", "longIRI14d1", "longIRI14d2", 
                               "shortIRI1d1", "shortIRI1d2", "shortIRI3d1", "shortIRI3d2", "shortIRI14d1", "shortIRI14d2",
                               "ESRRA", "G2CA-E",
                               "GSE151658-LPS1hr", "GSE151658-LPS4hr", "GSE151658-LPS16hr", "GSE151658-LPS27hr", 
                               "GSE151658-LPS36hr", "GSE151658-LPS48hr", 
                               "HK2596", "HK2844", "HK2862", "HK2886", 
                               "GSE184652-dbdbPBSrep1", "N3024", "GSE184652-dbdbPBSrep2", "B3025", "GSE184652-dbdbPBSrep3", "A3026", 
                               "GSE184652-dbdbPBSrep4", "B3027", "GSE184652-dbdbPBSrep5", "C3028", "GSE184652-dbdbPBSrep6", "A3029")
## RNA control datasets
rna.control.samples <- base::c("G2CA-C",
                               "norm1", "norm2", "norm3", "norm4", "norm5", "MMctrl7", "MMctrl8", 
                               "GSE151658-LPS0hr",
                               "HK2898", "HK2899", 
                               "GSE184652-dbmPBSrep1", "E3019", "GSE184652-dbmPBSrep2", "A3020", "GSE184652-dbmPBSrep3", "F3021", 
                               "GSE184652-dbmPBSrep4", "G3022", "GSE184652-dbmPBSrep5", "H3023")

## 30 datasets
select.b6.mouse.scRNA <- base::c("FA1", "FA2", "UUO1", "UUO2", 
                                 "longIRI1d1", "longIRI1d2", "longIRI3d1", "longIRI3d2", "longIRI14d1", "longIRI14d2", 
                                 "shortIRI1d1", "shortIRI1d2", "shortIRI3d1", "shortIRI3d2", "shortIRI14d1", "shortIRI14d2", 
                                 "ESRRA", 
                                 "norm1", "norm2", "norm3", "norm4", "norm5", "MMctrl7", 
                                 gse151658)
## 6 datasets
select.fvb.mouse.scRNA <- base::c("Notch3", "PGC1a-1", "PGC1a-2", "G2CA-E", "G2CA-C", "MMctrl8")
## 9 datasets
select.mouse.snRNA <- base::c("GSE184652-dbdbPBSrep1", "GSE184652-dbdbPBSrep2", "GSE184652-dbdbPBSrep3", "GSE184652-dbdbPBSrep6", 
                              "GSE184652-dbmPBSrep1", "GSE184652-dbmPBSrep2", "GSE184652-dbmPBSrep3", "GSE184652-dbmPBSrep4", 
                              "GSE184652-dbmPBSrep5")
## 6 datasets
select.human.snRNA <- base::c("HK2596", "HK2844", "HK2862", "HK2886", "HK2898", "HK2899")

### Markers ------
default.known.markers <- base::c(
  "Nrp1", "Kdr", # GEC, Endo
  "Ehd3", # GEC
  "Igfbp3", # Endo
  "Nphs1", "Nphs2", # Podo
  "Slc27a2", "Lrp2", # PT, "Aqp1" is also a good marker
  "Slc5a2", "Snhg11", "Slc5a12", # PT S1
  "Slc22a6", "Slc13a3", # PT S2
  "Slc7a13", "Atp11a", "Slc22a30", # PT S3
  "Slc4a11", "Ptgds", "Cp", # DLOH
  "Slc12a1", "Ppp1r1b", # ALOH
  "Slc12a3", "Pvalb", # DCT, CNT
  "Trpv5", # CNT
  "Aqp2", "Hsd11b2", # PC
  "Insrr", "Rhbg", # Trans
  "Atp6v1g3", "Atp6v0d2", # IC
  "Slc4a1", "Aqp6", # A IC
  "Slc26a4", "Hmx2", # B IC
  "Col1a1", "Col3a1", "Vim", "Fn1", # Fib
  "Plac8", "S100a4", # Mono
  "F13a1", "Chil3", # Mono Ly6Chi
  "Ace", "Treml4", # Mono Ly6Clo
  "Clec10a", # DC
  "Cd209a", # DC 11b+
  "Xcr1", # DC 11b-
  "Siglech", "Ccr9", # pDC
  "Fcer1a", "Mcpt8", "Csrp3", "Cd200r3", "Cyp11a1", "Ms4a2", # Baso
  "C1qa", "C1qb", "C1qc", # Macro
  "S100a8", "S100a9", "Hp", # Granul
  "Cd79a", "Cd79b", "Igkc", # B lymph
  "Ebf1", "Ighm", # B1
  "Igha", "Jchain", "Derl3", # B2
  "Ltb", "Cxcr6", # T lymph
  "Lef1", # T naive
  "Cd28", # T mem
  "Icos", "Rora", "Actn2", "Ly6g5b", # T gd
  "Cd3g", "Cd3d", "Cd8a", # CD8 effector
  "Ncr1", "Ccl5", "Nkg7", # NK
  "Gzma", # NK2
  "Mki67", "Cdca3", "Stmn1", "Lockd", # Novel
  "Top2a", # Novel2
  "Pdpn", "Wt1", "Mafb", "Synpo", "Cdkn1c", "Ptpro", # podocyte (DOI: 10.1681/asn.2020020220)
  "Cldn1", "Pax8", # + podocyte markers, parietal epithelial cell (PEC)
  "Pdgfrb", "Gata3", "Des", "Itga8", # mesangial (DOI: 10.1681/asn.2020020220)
  "Plvap", "Prkca", "Art3", "Nt5e", # mesangial (DOI: 10.1681/asn.2020020220, Fig. 1)
  "Flt1", "Tie1", "Pecam1", "Emcn", # endothelial (DOI: 10.1681/asn.2020020220)
  "Acta2", "Myh11", "Tagln", # smooth muscle cell (SMC) (DOI: 10.1681/asn.2020020220)
  "Ren1", "Akr1b7", "Rgs5", "Rergl", "Map3k7cl", # SMC (DOI: 10.1681/asn.2020020220, Fig. 1)
  "Ptprc", "Lyz1", "Csf1r", "Itgam", "Ms4a1", # immune (DOI: 10.1681/asn.2020020220)
  "Fxyd2", "Slc14a2", "Aqp1", "Umod" # tubular epithelial cell (TEC) (DOI: 10.1681/asn.2020020220)
)

### Local settings ------
picked.name <- proj.name # change here to analyze different data
if (!(picked.name %in% base::names(x = all.mice))) {
  base::stop(base::sprintf("could not identify the picked name \"%s\".", picked.name))
}

curr.mice <- all.mice[[picked.name]]
curr.out.dir <- base::file.path(cwd, time.stamp, picked.name)
rds.dir <- base::file.path(curr.out.dir, "RDS")
rst.dir <- base::file.path(curr.out.dir, "results")
rst.deg.dir <- base::file.path(rst.dir, "DEG")
rst.mac.dir <- base::file.path(rst.dir, "macs2")
rst.bul.dir <- base::file.path(rst.dir, "bulk")
rst.pub.dir <- base::file.path(rst.dir, "publication")
rst.clu.dir <- base::file.path(rst.dir, "clusters")
rst.tra.dir <- base::file.path(rst.dir, "trajectory")
rst.mo2.dir <- base::file.path(rst.tra.dir, "monocle2")
rst.scv.dir <- base::file.path(rst.tra.dir, "scvelo")
rst.dap.dir <- base::file.path(rst.dir, "DAP")
rst.hom.dir <- base::file.path(rst.dap.dir, "HOMER")
rst.gre.dir <- base::file.path(rst.dap.dir, "GREAT")
rst.wgc.dir <- base::file.path(rst.dir, "WGCNA")
rst.qc.dir <- base::file.path(rst.dir, "QC")
fig.dir <- base::file.path(curr.out.dir, "figures")
fig.qc.dir <- base::file.path(fig.dir, "QC")
fig.pc.dir <- base::file.path(fig.dir, "PC")
fig.int.dir <- base::file.path(fig.dir, "integration")
fig.clu.dir <- base::file.path(fig.dir, "clusters")
fig.deg.dir <- base::file.path(fig.dir, "DEG")
fig.pub.dir <- base::file.path(fig.dir, "publication")
fig.bul.dir <- base::file.path(fig.dir, "bulk")
fig.tra.dir <- base::file.path(fig.dir, "trajectory")
fig.des.dir <- base::file.path(fig.tra.dir, "destiny")
fig.vel.dir <- base::file.path(fig.tra.dir, "velocyto")
fig.mo2.dir <- base::file.path(fig.tra.dir, "monocle2")
fig.scv.dir <- base::file.path(fig.tra.dir, "scvelo")
fig.pat.dir <- base::file.path(fig.dir, "pathway")
fig.wgc.dir <- base::file.path(fig.dir, "WGCNA")
base::dir.create(path = curr.out.dir, showWarnings = T, recursive = T, mode = "0755")
base::dir.create(path = rds.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = rst.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = rst.deg.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = rst.mac.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = rst.bul.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = rst.pub.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = rst.clu.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = rst.tra.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = rst.mo2.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = rst.scv.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = rst.dap.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = rst.hom.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = rst.gre.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = rst.wgc.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = rst.qc.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.qc.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.pc.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.int.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.clu.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.deg.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.pub.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.bul.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.tra.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.des.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.vel.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.mo2.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.scv.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.pat.dir, showWarnings = T, recursive = F, mode = "0755")
base::dir.create(path = fig.wgc.dir, showWarnings = T, recursive = F, mode = "0755")

base::source(file = base::file.path(src.dir, "auxiliary_functions.R"))
