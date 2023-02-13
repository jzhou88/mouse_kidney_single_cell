# bulk_analysis.R
base::cat("\nRunning \"bulk_analysis.R\" ...\n")

iri.bulk <- base::c("longIRI14d", "shortIRI14d")
bulk.models <- base::c("FAN", "UUO", iri.bulk, "Cisplatin", "APOL1")
genes.results.files <- base::list()
sample.names <- base::list()
sample.conditions <- base::list()
model.contrasts <- base::list()

### FAN
genes.results.files[["FAN"]] <- base::c(
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/FAN/rsem/L4C1p02CTL_S25_R1_001.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/FAN/rsem/L4C1p04CTL_S26_R1_001.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/FAN/rsem/L4C1p05CTL_S27_R1_001.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/FAN/rsem/NI19FAN_S30_R1_001.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/FAN/rsem/PS2p45FAN_S28_R1_001.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/FAN/rsem/PS2p46FAN_S29_R1_001.genes.results"
)
sample.names[["FAN"]] <- base::c("L4C1p02CTL", "L4C1p04CTL", "L4C1p05CTL", "NI19FAN", "PS2p45FAN", "PS2p46FAN")
sample.conditions[["FAN"]] <- base::factor(x = base::c("Control", "Control", "Control", "FAN", "FAN", "FAN"))
model.contrasts[["FAN"]] <- base::c("condition", "FAN", "Control")

### UUO
genes.results.files[["UUO"]] <- base::c(
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/UUO/rsem/PS63SHAM_ACAGTG_L007_R1_001.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/UUO/rsem/PS64SHAM_GCCAAT_L007_R1_001.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/UUO/rsem/PS8SHAM_TGACCA_L007_R1_001.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/UUO/rsem/WT3UUO_ATCACG_L007_R1_001.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/UUO/rsem/WT5UUO_CGATGT_L007_R1_001.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/UUO/rsem/WT6UUO_TTAGGC_L007_R1_001.genes.results"
)
sample.names[["UUO"]] <- base::c("PS63SHAM", "PS64SHAM", "PS8SHAM", "WT3UUO", "WT5UUO", "WT6UUO")
sample.conditions[["UUO"]] <- base::factor(x = base::c("UUO", "UUO", "UUO", "Control", "Control", "Control"))
model.contrasts[["UUO"]] <- base::c("condition", "UUO", "Control")

### longIRI14d: IRI9, IRI12
genes.results.files[["longIRI14d"]] <- base::c(
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/paired-end/IRI/rsem/IRI9_CRRA200010445-1a_H55HFDSXY_L1.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/paired-end/IRI/rsem/IRI12_CRRA200010448-1a_H55HFDSXY_L1.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/paired-end/Cisplatin/rsem/A1.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/paired-end/Cisplatin/rsem/A2.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/paired-end/Cisplatin/rsem/A3.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/paired-end/Cisplatin/rsem/A4.genes.results"
)
sample.names[["longIRI14d"]] <- base::c("IRI9", "IRI12", "A1", "A2", "A3", "A4")
sample.conditions[["longIRI14d"]] <- base::factor(x = base::c("longIRI14d", "longIRI14d", "Control", "Control", "Control", "Control"))
model.contrasts[["longIRI14d"]] <- base::c("condition", "longIRI14d", "Control")

### shortIRI14d: IRI8, IRI10
genes.results.files[["shortIRI14d"]] <- base::c(
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/paired-end/IRI/rsem/IRI8_CRRA200010444-1a_H55HFDSXY_L1.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/paired-end/IRI/rsem/IRI10_CRRA200010446-1a_H55HFDSXY_L1.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/paired-end/Cisplatin/rsem/A1.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/paired-end/Cisplatin/rsem/A2.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/paired-end/Cisplatin/rsem/A3.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/paired-end/Cisplatin/rsem/A4.genes.results"
)
sample.names[["shortIRI14d"]] <- base::c("IRI8", "IRI10", "A1", "A2", "A3", "A4")
sample.conditions[["shortIRI14d"]] <- base::factor(x = base::c("shortIRI14d", "shortIRI14d", "Control", "Control", "Control", "Control"))
model.contrasts[["shortIRI14d"]] <- base::c("condition", "shortIRI14d", "Control")

### Cisplatin
genes.results.files[["Cisplatin"]] <- base::c(
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/paired-end/Cisplatin/rsem/A1.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/paired-end/Cisplatin/rsem/A2.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/paired-end/Cisplatin/rsem/A3.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/paired-end/Cisplatin/rsem/A4.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/paired-end/Cisplatin/rsem/B1.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/paired-end/Cisplatin/rsem/B2.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/paired-end/Cisplatin/rsem/B3.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/paired-end/Cisplatin/rsem/B4.genes.results"
)
sample.names[["Cisplatin"]] <- base::c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4")
sample.conditions[["Cisplatin"]] <- base::factor(x = base::c("Control", "Control", "Control", "Control", 
                                                             "Cisplatin", "Cisplatin", "Cisplatin", "Cisplatin"))
model.contrasts[["Cisplatin"]] <- base::c("condition", "Cisplatin", "Control")

### APOL1
genes.results.files[["APOL1"]] <- base::c(
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/APOL1/rsem/C1.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/APOL1/rsem/C2.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/APOL1/rsem/C3.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/APOL1/rsem/WT1.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/APOL1/rsem/WT2.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/APOL1/rsem/WT3.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/APOL1/rsem/G1.1.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/APOL1/rsem/G1.2.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/APOL1/rsem/G2.1.genes.results", 
  "/home/jzhou88/projects/KidneyDisease/2022-01-26/bulkRNAseq/single-end/APOL1/rsem/G2.2.genes.results"
)
sample.names[["APOL1"]] <- base::c("G0.1", "G0.2", "G0.3", "WT1", "WT2", "WT3", "G1.1", "G1.2", "G2.1", "G2.2")
sample.conditions[["APOL1"]] <- base::factor(x = base::c("Control", "Control", "Control", "Control", "Control", "Control", "APOL1", "APOL1", "APOL1", "APOL1"))
model.contrasts[["APOL1"]] <- base::c("condition", "APOL1", "Control")

count.list <- base::list()
tpm.list <- base::list()
count.matrix <- NULL
tpm.matrix <- NULL
for (curr.model in bulk.models) {
  base::cat(base::sprintf(fmt = "\nProcessing %s ...\n", curr.model))
  curr.files <- genes.results.files[[curr.model]]
  curr.names <- sample.names[[curr.model]]
  base::cat("\n")
  base::print(x = curr.files)
  base::cat("\n")
  base::print(x = curr.names)
  
  count.data <- base::data.frame(data.table::fread(curr.files[1]))[base::c(1,5)]
  tpm.data <- base::data.frame(data.table::fread(curr.files[1]))[base::c(1,6)]
  base::cat("\n")
  base::print(x = utils::head(x = count.data))
  base::cat("\n")
  base::print(x = utils::head(x = tpm.data))
  ## Loop and read the 5th and 6th columns in the remaining files
  for(i in 2:base::length(x = curr.files)) {
    count.data <- base::cbind(count.data, base::data.frame(data.table::fread(curr.files[i]))[5])
    tpm.data <- base::cbind(tpm.data, base::data.frame(data.table::fread(curr.files[i]))[6])
    base::cat("\n")
    base::print(x = utils::head(x = count.data))
    base::cat("\n")
    base::print(x = utils::head(x = tpm.data))
  }
  
  base::colnames(x = count.data) <- base::c("GeneID", curr.names)
  base::colnames(x = tpm.data) <- base::c("GeneID", curr.names)
  base::cat("\n")
  base::print(x = utils::head(x = count.data))
  base::cat("\n")
  base::print(x = utils::head(x = tpm.data))
  base::rownames(x = count.data) <- count.data$GeneID
  base::rownames(x = tpm.data) <- tpm.data$GeneID
  base::cat("\n")
  base::print(x = utils::head(x = count.data))
  base::cat("\n")
  base::print(x = utils::head(x = tpm.data))
  count.data <- count.data[,base::c(2:base::ncol(x = count.data))]
  tpm.data <- tpm.data[,base::c(2:base::ncol(x = tpm.data))]
  base::cat("\n")
  base::print(x = utils::head(x = count.data))
  base::cat("\n")
  base::print(x = utils::head(x = tpm.data))
  
  ## Save original counts and TPMs.
  count.list[[curr.model]] <- count.data
  tpm.list[[curr.model]] <- tpm.data
  
  if (base::is.null(x = count.matrix)) {
    if (curr.model %in% iri.bulk) {
      count.matrix <- count.data[,base::c(1,2)]
    } else {
      count.matrix <- count.data
    }
  } else {
    if (curr.model %in% iri.bulk) {
      count.matrix <- base::cbind(count.matrix, count.data[,base::c(1,2)])
    } else {
      count.matrix <- base::cbind(count.matrix, count.data)
    }
  }
  if (base::is.null(x = tpm.matrix)) {
    if (curr.model %in% iri.bulk) {
      tpm.matrix <- tpm.data[,base::c(1,2)]
    } else {
      tpm.matrix <- tpm.data
    }
  } else {
    if (curr.model %in% iri.bulk) {
      tpm.matrix <- base::cbind(tpm.matrix, tpm.data[,base::c(1,2)])
    } else {
      tpm.matrix <- base::cbind(tpm.matrix, tpm.data)
    }
  }
  base::cat("\n")
  base::print(x = utils::head(x = count.matrix))
  base::cat("\n")
  base::print(x = utils::head(x = tpm.matrix))
}
base::saveRDS(object = count.list, file = base::file.path(rds.dir, "bulkRNAseq_count.RDS"))
base::saveRDS(object = tpm.list, file = base::file.path(rds.dir, "bulkRNAseq_TPM.RDS"))
base::saveRDS(object = count.matrix, file = base::file.path(rds.dir, "bulkRNAseq_countmatrix.RDS"))
base::saveRDS(object = tpm.matrix, file = base::file.path(rds.dir, "bulkRNAseq_tpmmatrix.RDS"))
openxlsx::write.xlsx(x = count.list, file = base::file.path(rst.bul.dir, "bulkRNAseq_count.xlsx"), 
                     startCol = 1, startRow = 1, colNames = T, rowNames = T, overwrite = T)
openxlsx::write.xlsx(x = tpm.list, file = base::file.path(rst.bul.dir, "bulkRNAseq_TPM.xlsx"), 
                     startCol = 1, startRow = 1, colNames = T, rowNames = T, overwrite = T)
openxlsx::write.xlsx(x = count.matrix, file = base::file.path(rst.bul.dir, "bulkRNAseq_countmatrix.xlsx"), 
                     startCol = 1, startRow = 1, colNames = T, rowNames = T, overwrite = T)
openxlsx::write.xlsx(x = tpm.matrix, file = base::file.path(rst.bul.dir, "bulkRNAseq_tpmmatrix.xlsx"), 
                     startCol = 1, startRow = 1, colNames = T, rowNames = T, overwrite = T)

tpm.thresh <- 1
sample.thresh <- 2
p.thresh <- 0.05

dds.list <- base::list()
res.list <- base::list()
res.tab.list <- base::list()
for (curr.model in bulk.models) {
  base::cat(base::sprintf(fmt = "\nProcessing %s ...\n", curr.model))
  count.data <- count.list[[curr.model]]
  tpm.data <- tpm.list[[curr.model]]
  condition <- sample.conditions[[curr.model]]
  contrast <- model.contrasts[[curr.model]]
  base::cat("\n")
  base::print(x = base::nrow(x = count.data))
  base::print(x = utils::head(x = count.data))
  base::cat("\n")
  base::print(x = base::nrow(x = tpm.data))
  base::print(x = utils::head(x = tpm.data))
  base::cat("\n")
  base::print(x = condition)
  base::cat("\n")
  base::print(x = contrast)
  
  ## Filter out genes with low expression.
  keep <- (base::rowSums(x = (tpm.data >= tpm.thresh)) >= sample.thresh)
  count.data <- count.data[keep,]
  base::cat("\n")
  base::print(x = base::nrow(x = count.data))
  base::print(x = utils::head(x = count.data))
  
  ## Round to integer.
  count.data <- base::round(x = count.data)
  
  sample.data <- base::data.frame(condition, row.names = base::colnames(x = count.data))
  base::cat("\n")
  base::print(x = sample.data)
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = count.data, colData = sample.data, design = ~ condition)
  base::cat("\n")
  base::print(x = dds)
  dds <- DESeq2::DESeq(object = dds)
  base::cat("\n")
  base::print(x = dds)
  
  dds.list[[curr.model]] <- dds
  
  res <- DESeq2::results(object = dds, contrast = contrast, alpha = p.thresh)
  base::cat("\n")
  base::print(x = DESeq2::summary(object = res))
  res$GeneSymbol <- stringr::str_split(string = base::rownames(x = res), pattern = "_", simplify = T)[,2]
  base::cat("\n")
  base::print(x = utils::head(x = res))
  
  res.list[[curr.model]] <- res
  
  res.tab <- base::as.data.frame(x = res)
  res.tab <- res.tab[base::order(res.tab$log2FoldChange, decreasing = T),]
  
  res.tab.list[[curr.model]] <- res.tab
}
base::saveRDS(object = dds.list, file = base::file.path(rds.dir, "bulkRNAseq_DESeqDataSet.RDS"))
base::saveRDS(object = res.list, file = base::file.path(rds.dir, "bulkRNAseq_DESeqResult.RDS"))
openxlsx::write.xlsx(x = res.tab.list, file = base::file.path(rst.bul.dir, "bulkRNAseq_DESeqResult.xlsx"), 
                     startCol = 1, startRow = 1, colNames = T, rowNames = T, overwrite = T)
