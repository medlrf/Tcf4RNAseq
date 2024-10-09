###### 10xGenomics WTA 3' scRNA-seq data analysis reproduction ######
# Create: 2024-10-08
# Author: Ruofan Li, reproduction from Marius Stephan
# Last edit: 2024-10-08

##### Load required libraries #####
library(Seurat) # Single-cell RNA-seq analysis
library(tidyverse) # Data manipulation
library(patchwork) # Combining multiple ggplot2 plots together
library(cowplot) # Combining multiple ggplot2 plots together  (alternative to patchwork)

# install.packages("BiocManager")
# BiocManager::install("fgsea")
library(fgsea) # Fast Gene Set Enrichment Analysis
# BiocManager::install("Nebulosa")
library(Nebulosa) # Gene set enrichment analysis

BiocManager::install(c("scds", "SingleCellExperiment", "SCMarker", "One2One"))
library(scds) # Single-cell differential splicing analysis
library(SingleCellExperiment) # Single-cell data representation
#library(SCMarker) # Single-cell marker genes identification
#library(One2One) # Single-cell data integration
##### Setting the random seed for reproducibility #####
set.seed(1234) 

##### User-defined functions #####
prepSCdata <- function(filename, projectname, samplename){
  print("Importing Data ...")
  raw <- Read10X_h5(filename)
  d1 <- CreateSeuratObject(counts = raw, project = projectname, min.cells = 3, min.features = 200)
  
      # !!!  Make sure that the gene names are unique
      rownames(d1) <- toupper(rownames(d1))
      #d1 <- subset(d1, cells = sample(Cells(d1), 5000))  
  
  # Doublet identification with SCDS
  print("Identifying potential doublets ...")
  sce <- SingleCellExperiment(list(counts = raw)) # SCDS expects SingleCellExperiment object
  sce = cxds_bcds_hybrid(sce, estNdbl = T)
  
  # Add Doublet calls to meta.data table
  d1@meta.data$doubletHybridScore <- sce$hybrid_score
  d1@meta.data$doubletCall <- sce$hybrid_call
  
  print(paste0("There were ", sum(d1$doubletCall), " out of ", length(d1$doubletCall), " cells flagged as doublets!"))
  
  # Add mito gene percentage
  print("Finding mitochondrial transcript level ...")
  d1[["percent.mt"]] <- PercentageFeatureSet(d1, pattern = "^mt-")
  
  # # Cell cycle scoring
  # print("Scoring cell cycle ...")
  # s.genes <- cc.genes$s.genes
  # g2m.genes <- cc.genes$g2m.genes
  # d1 <- CellCycleScoring(d1,s.genes,g2m.genes)
      
  # Adding sample name
  d1@meta.data$group <- samplename
  
  print("Done!")
  return(d1)
  
}

##### Load data #####
# Tcf4_vHIP_TG_SD_filtered <- Read10X_h5("/Users/medlrf/Documents/Github/Tcf4RNAseq/vHIP Count Matrices/Tcf4_vHIP_TG_SD_filtered.h5")
# Tcf4_vHIP_WT_HC_filtered <- Read10X_h5("/Users/medlrf/Documents/Github/Tcf4RNAseq/vHIP Count Matrices/Tcf4_vHIP_WT_HC_filtered.h5")

d1 = prepSCdata("/Users/medlrf/Documents/Github/Tcf4RNAseq/mPFC Count Matrices/TG_HC_counts.h5", "Tcf4 SD mPFC", "TG_HC")
dList <- list(
  d1 = prepSCdata("/Users/medlrf/Documents/Github/Tcf4RNAseq/mPFC Count Matrices/TG_HC_counts.h5", "Tcf4 SD mPFC", "TG_HC")
)

# QC plots
for (d in dList){
  print(VlnPlot(d, features = c("nFeature_RNA", "nCount_RNA","percent.mt", "doubletHybridScore"), group.by = "doubletCall", ncol = 4, pt.size = 0)+geom_hline(yintercept = 0.9))
}


# # Filter and normalize Expression Data
# dList <- lapply(X = dList, FUN = function(x) {
#   x <- SCTransform(x, vars.to.regress = c("percent.mt","G2M.Score","S.Score"), vst.flavor = "v2")
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
# })

# Integrate data
IntFeats <- SelectIntegrationFeatures(object.list = dList, nfeatures = 2000)
dList <- PrepSCTIntegration(object.list = dList, anchor.features = IntFeats)
IntAnchors <- FindIntegrationAnchors(object.list = dList, normalization.method = "SCT", anchor.features = IntFeats)
dm <- IntegrateData(anchorset = IntAnchors, normalization.method = "SCT")
rm(IntAnchors,dList)

# Run PCA
dm <- RunPCA(dm)
ElbowPlot(dm,ndims = 30)
