rm(list=ls())

suppressMessages(suppressWarnings(library(writexl)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(DESeq2)))
suppressMessages(suppressWarnings(library(RColorBrewer)))
suppressMessages(suppressWarnings(library(mclust)))
suppressMessages(suppressWarnings(library(cluster)))
suppressMessages(suppressWarnings(library(dendextend)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(pheatmap)))
suppressMessages(suppressWarnings(library(grid)))
suppressMessages(suppressWarnings(library(clusterProfiler)))
suppressMessages(suppressWarnings(library(org.Hs.eg.db)))
suppressMessages(suppressWarnings(library(grid)))
suppressMessages(suppressWarnings(library(gridExtra)))
suppressMessages(suppressWarnings(library(DT)))                 



set.seed(42) 

cell_names <- c("Astro", "Excite", "Inhibit", "Micro", "Oligo", "OPC")

rename_list <- list("Astro", "Astro", "Excite", "Inhibit", "Micro", "Micro", "Oligo", "OPC")
names(rename_list) <- c("Astro", "Astrocyte", "Excite", "Inhibit", "Microglia", "Micro", "Oligo", "OPC")


all_canonical_genes <- c("AQP4", "SLC1A2", "SLC4A4", "GFAP", "ALDH1L1", "GJA1",
                         "GAD1", "GAD2", "GRIK2", "KCNC2", "CDH9", "KIT", "SST", "VIP", "SV2C",
                         "SLC17A7", "SATB2", "CAMK2A", "CBLN2", "IL1RAPL2", "TRPM3", "SEMA3E", "SYT1", "RBFOX3", "CUX2", "RORB", "TLE4",
                         "TYROBP", "APBB1IP", "P2RY12", "CSF1R", "CX3CR1", "PTPRC", "DOCK8", "HLA.DRA",
                         "OPALIN", "MOBP", "MBP", "PLP1", "ST18", "MOG",
                         "PDGFRA", "PCDH15", "VCAN", "CSPG4")


all_canonical_gene_list <- list(c("AQP4", "SLC1A2", "SLC4A4", "GFAP", "ALDH1L1", "GJA1"),
                                c("GAD1", "GAD2", "GRIK2", "KCNC2", "CDH9", "KIT", "SST", "VIP", "SV2C"),
                                c("SLC17A7", "SATB2", "CAMK2A", "CBLN2", "IL1RAPL2", "TRPM3", "SEMA3E", "SYT1", "RBFOX3", "CUX2", "RORB", "TLE4"),
                                c("TYROBP", "APBB1IP", "P2RY12", "CSF1R", "CX3CR1", "PTPRC", "DOCK8", "HLA.DRA"),
                                c("OPALIN", "MOBP", "MBP", "PLP1", "ST18", "MOG"),
                                c("PDGFRA", "PCDH15", "VCAN", "CSPG4"))

names(all_canonical_gene_list) <- c("Astro", "Inhibit", "Excite", "Micro", "Oligo", "OPC")

