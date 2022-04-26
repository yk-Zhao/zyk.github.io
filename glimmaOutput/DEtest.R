
library(edgeR)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Glimma")
library(Glimma)

dge <- readRDS(system.file("RNAseq123/dge.rds", package = "Glimma"))

glimmaVolcano(efit, dge = dge)

design <- readRDS(
  system.file("RNAseq123/design.rds", package = "Glimma"))
contr.matrix <- readRDS(
  system.file("RNAseq123/contr.matrix.rds", package = "Glimma"))

v <- voom(dge, design)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)

glimmaVolcano(efit, dge = dge)
glimmaMA(efit, dge = dge)

