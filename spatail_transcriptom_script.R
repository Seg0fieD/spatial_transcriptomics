#===============================================================================================
#Library =====================================================
suppressPackageStartupMessages({	
	library(Seurat)
	library(SeuratData)
	library(SeuratWrappers)
	library(ggplot2)
	library(dplyr)
	library(patchwork)
	library(SCDC)
	library(gridExtra)
	library(grid)
	library(CellChat)
	library(Biobase)
	library(glmGamPoi)
	library(xbioc)
	library(BiocNeighbors)
	library(Matrix)
	library(spatstat)
	options(stringsAsFactors = FALSE)
  })	

#===============================================================================================
#===============================================================================================
# Set working Directory ======================================
setwd("~/Spatial_transcriptome/")

# Read in the image file ======================================
image_1 <- Read10X_Image("~/Spatial_transcriptome/Section_1/spatial")
image_2 <- Read10X_Image("~/Spatial_transcriptome/Section_2/spatial")


# Load-10X--spatial Data ========================================
spdata1 <- Load10X_Spatial(data.dir = "~/Spatial_transcriptome/Section_1/",
                          filename = "V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.h5",
                          slice= "slice1", image = image_1)

spdata2 <- Load10X_Spatial(data.dir = "~/Spatial_transcriptome/Section_2",
                           filename = "V1_Mouse_Brain_Sagittal_Posterior_Section_2_filtered_feature_bc_matrix.h5",
                           slice= "slice2", image = image_2)
# Set the indents 
spdata1 <- SetIdent(spdata1, value = "slice1")
spdata2 <- SetIdent(spdata2, value = "slice2")

spdata1.tab <- data.frame(spdata1@meta.data)
spdata2.tab <- data.frame(spdata2@meta.data)


spdata1@images$slice1
spdata2@images$slice2

#===============================================================================================
#===============================================================================================
# Adding Mitochondrial Percentage==========================================
spdata1[["percent.mt"]] <- PercentageFeatureSet(spdata1, pattern = "^mt.")
spdata2[["percent.mt"]] <- PercentageFeatureSet(spdata2, pattern = "^mt.")

spdata1[[]]
spdata2[[]]

#===============================================================================================
#===============================================================================================
# Visualization Pre-Filtration ==================================
## = = = = = Mitochondrial Percentage= = = = = 
### -- -- -- -- -- -- section-1 -- -- -- -- -- --
png(file = "png/Prefilter-percent-mito-section-1.png", width = 800, height = 500)
wrap_plots(VlnPlot(spdata1, features = "percent.mt", pt.size = 0.1, )+ NoLegend(),
          SpatialFeaturePlot(spdata1, features = "percent.mt")+ theme(legend.position = "right"))
dev.off()

### -- -- -- -- -- -- section-2 -- -- -- -- -- --
png(file = "png/Prefilter-percent-mito-section-2.png", width = 800, height = 500)
wrap_plots(VlnPlot(spdata2, features = "percent.mt", pt.size = 0.1)+ NoLegend(),
            SpatialFeaturePlot(spdata2, features = "percent.mt")+ theme(legend.position = "right"))
dev.off()

## = = = = = nFeature_Spatial = = = = = 
### -- -- -- -- -- -- section-1 -- -- -- -- -- --
png(file = "png/Prefilter-nFeature_Spatial-section-1.png", width = 800, height = 500)
wrap_plots(VlnPlot(spdata1, features = "nFeature_Spatial", pt.size = 0.1, )+ NoLegend(),
           SpatialFeaturePlot(spdata1, features = "nFeature_Spatial")+ theme(legend.position = "right"))
dev.off()

### -- -- -- -- -- -- section-2-- -- -- -- -- --
png(file = "png/Prefilter-nFeature_Spatial-section-2.png", width = 800, height = 500)
wrap_plots(VlnPlot(spdata2, features = "nFeature_Spatial", pt.size = 0.1)+ NoLegend(),
           SpatialFeaturePlot(spdata2, features = "nFeature_Spatial")+ theme(legend.position = "right"))
dev.off()

## = = = = = nCount_Spatial = = = = =
### -- -- -- -- -- -- section-1 -- -- -- -- -- --
png(file = "png/Prefilter-nCount_Spatial-section-1.png", width = 800, height = 500)
wrap_plots(VlnPlot(spdata1, features = "nCount_Spatial", pt.size = 0.1, )+ NoLegend(),
           SpatialFeaturePlot(spdata1, features = "nCount_Spatial")+ theme(legend.position = "right"))
dev.off()

### -- -- -- -- -- -- section-2 -- -- -- -- -- --
png(file = "png/Prefilter-nCount_Spatial-section-2.png", width = 800, height = 500)
wrap_plots(VlnPlot(spdata2, features = "nCount_Spatial", pt.size = 0.1)+ NoLegend(),
           SpatialFeaturePlot(spdata2, features = "nCount_Spatial")+ theme(legend.position = "right"))
dev.off()

# Pick up random gene to visualize 
png(file =  "png/Hpca_&_Ttr-gene-section1.png", width = 1000, height = 600) 
SpatialFeaturePlot(spdata1, features = c("Hpca","Ttr"), pt.size.factor = 1.6, alpha = c(0.1, 1.6))
dev.off()

png(file =  "png/Hpca_&_Ttr-gene-section2.png", width = 1000, height = 600) 
SpatialFeaturePlot(spdata2, features = c("Hpca","Ttr"), pt.size.factor = 1.6, alpha = c(0.1, 1.6))
dev.off()


#===============================================================================================
#===============================================================================================
# Filtering the data=============================================
spdata1 <- subset(spdata1, nFeature_Spatial > 200 & nFeature_Spatial < 8000 
                  & nCount_Spatial < 40000 & percent.mt > 5 & percent.mt < 25)

spdata2 <- subset(spdata2, nFeature_Spatial > 200 & nFeature_Spatial < 8000 
                  & nCount_Spatial < 40000 & percent.mt > 5 & percent.mt < 25)


#===============================================================================================
#===============================================================================================
# Visualization Post-filtration===================================
png(file="png/Voilin-plot-for-section1-POST-filtration.png", width = 800, height = 500)
VlnPlot(spdata1, features = c("nCount_Spatial","nFeature_Spatial","percent.mt"), pt.size = 0.1,)+ NoLegend()
dev.off()

png(file =  "png/Voilin-plot-for-section2-POST-filtration.png", width = 800, height = 500) 
VlnPlot(spdata2, features = c("nCount_Spatial","nFeature_Spatial","percent.mt"), pt.size = 0.1)+ NoLegend()
dev.off()


#===============================================================================================
#===============================================================================================
# SC Transformation======================================

spdata1 <- SCTransform(spdata1, assay = "Spatial", verbose = FALSE)
spdata2 <- SCTransform(spdata2, assay = "Spatial", verbose = FALSE)

#===============================================================================================
#===============================================================================================
# Dimensionality Reduction=====================================
##   +PCA =============================================
spdata1 <- RunPCA(spdata1, assay = "SCT", verbose = FALSE, ndims = 50)
spdata2 <- RunPCA(spdata2, assay = "SCT", verbose = FALSE, ndims = 50)

##   +Elbow Plot ==============
png(file =  "png/ Elbow plot-section-1.png", width = 1000, height = 600) 
ElbowPlot(spdata1, ndims = 50 )
dev.off()

png(file =  "png/ Elbow plot-section-2.png", width = 1000, height = 600) 
ElbowPlot(spdata2, ndims = 50 )
dev.off()

##   +Finding Neighbors================================
spdata1 <- FindNeighbors(spdata1, reduction = "pca", dims = 1:30)
spdata2 <- FindNeighbors(spdata2, reduction = "pca", dims = 1:30)


#===============================================================================================
#===============================================================================================
# Clustering===================================================
##   +Find_clusters==================
spdata1 <- FindClusters(spdata1, verbose = FALSE, resolution =  0.5)
spdata2 <- FindClusters(spdata2, verbose = FALSE, resolution =  0.5)

table(spdata1$SCT_snn_res.0.5)
table(spdata2$SCT_snn_res.0.5)

##   +UMAP===========================
spdata1 <- RunUMAP(spdata1, reduction = 'pca', dims = 1:30)
spdata2 <- RunUMAP(spdata2, reduction = 'pca', dims = 1:30)

##   +cluster_Visualization=========== 
### -- -- -- -- -- -- section-1 -- -- -- -- -- --
png(file =  "png/UMAP-resolution_0.5-section-1.png", width = 850, height = 600)
DimPlot(spdata1, reduction = 'umap', label = TRUE)+DarkTheme()
dev.off()

png(file =  "png/Sptial-Dim-plot-section-1.png", width = 900, height = 600)
SpatialDimPlot(spdata1, label = TRUE, label.size = 10, label.box = TRUE, label.color = "black")
dev.off()

png(file =  "png/Cell-highlight-section-1.png", width = 1500, height = 1000)
SpatialDimPlot(spdata1, cells.highlight = CellsByIdentities(object = spdata1, idents = c(2,1,4,3,5,8)), 
               facet.highlight = TRUE, ncol = 3)
dev.off()

### -- -- -- -- -- -- section-2 -- -- -- -- -- --

png(file =  "png/UMAP-resolution_0.5-section-2.png", width = 850, height = 600)
DimPlot(spdata2, reduction = 'umap', label = TRUE)+DarkTheme()
dev.off()

png(file =  "png/Sptial-Dim-plot-section-2.png", width = 900, height = 600)
SpatialDimPlot(spdata2, label = TRUE, label.size = 10, label.box = TRUE, label.color = "black")
dev.off()

png(file =  "png/Cell-highlight-section-2.png", width = 1500, height = 1000)
SpatialDimPlot(spdata2, cells.highlight = CellsByIdentities(object = spdata2, idents = c(2,1,4,3,5,8)), 
               facet.highlight = TRUE, ncol = 3)
dev.off()


#===============================================================================================
#===============================================================================================
# Differential Expression Analysis==============================
##   +Comparison Gene-Expression between clusters====
### - - - - - - - Marker from section-1 clusters- - - - - - - -
sp1mrk1 <- FindMarkers(spdata1, ident.1 = 0, min.pct = 0.25)
sp1mrk2 <- FindMarkers(spdata1, ident.1 = 3, min.pct = 0.25)
sp1mrk3 <- FindMarkers(spdata1, ident.1 = 6, min.pct = 0.25)

Mrk.t1 <- data.frame(head(sp1mrk1))
Mrk.t2 <- data.frame(head(sp1mrk2))
Mrk.t3 <- data.frame(head(sp1mrk3))


### - - - - - - - Marker from section-2 clusters- - - - - - - -
sp2mrk1 <- FindMarkers(spdata2, ident.1 = 0, min.pct = 0.25)
sp2mrk2 <- FindMarkers(spdata2, ident.1 = 1, min.pct = 0.25)
sp2mrk3 <- FindMarkers(spdata2, ident.1 = 3, min.pct = 0.25)

Mrk.t4 <- data.frame(head(sp2mrk1))
Mrk.t5 <- data.frame(head(sp2mrk2))
Mrk.t6 <- data.frame(head(sp2mrk3))


##   +Find All marker=====================
### -- -- -- -- -- -- section-1 -- -- -- -- -- -- 
sp1_all_mrk <- FindAllMarkers(spdata1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

png(file =  "png/ Section-1-All-Marker.png", width = 1500, height = 900) 
SpatialFeaturePlot(object = spdata1, features = rownames(sp1_all_mrk)[1:6], alpha =c(0.1,1), ncol = 3)
dev.off()

### -- -- -- -- -- -- section-2 -- -- -- -- -- -- 
sp2_all_mrk <- FindAllMarkers(spdata2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

png(file =  "png/ Section-2-All-Marker.png", width = 1500, height = 900) 
SpatialFeaturePlot(object = spdata2, features = rownames(sp2_all_mrk)[1:6], alpha =c(0.1,1), ncol = 3)
dev.off()

##   +Find spatially Variable Features============
spdata1 <- FindSpatiallyVariableFeatures(spdata1, assay = "SCT", 
                                        features = VariableFeatures(spdata1)[1:100], 
                                        selection.method = "markvariogram")

spdata2 <- FindSpatiallyVariableFeatures(spdata2, assay = "SCT", 
                                        features = VariableFeatures(spdata2)[1:100], 
                                        selection.method = "markvariogram")

##   +Visualization of Top Markers================
###  -- -- -- -- -- -- section-1 -- -- -- -- -- -- 
section_1.Top_Features <- SpatiallyVariableFeatures(spdata1, selection.method = "markvariogram")
feature.tab1 <- data.frame(section_1.Top_Features )
top3feature_sec1  <- head(SpatiallyVariableFeatures(spdata1, selection.method = "markvariogram"),3)

png(file =  "png/Top-3-spatially-varible-features-section-1.png", width = 1200, height = 700) 
SpatialFeaturePlot(spdata1, features = top3feature_sec1, ncol = 3, alpha = c(0.1, 1))
dev.off()

###  -- -- -- -- -- -- section-2 -- -- -- -- -- -- 
section_2.Top_Features <- SpatiallyVariableFeatures(spdata2, selection.method = "markvariogram")
feature.tab2 <- data.frame(section_2.Top_Features)
top3feature_sec2  <- head(SpatiallyVariableFeatures(spdata2, selection.method = "markvariogram"),3)

png(file =  "png/ Top-3-spatially-varible-features-section-2.png", width = 1200, height = 700) 
SpatialFeaturePlot(spdata2, features = top3feature_sec1, ncol = 3, alpha = c(0.1, 1))
dev.off()

# Save spatial data (section1 & section2)=====================================
saveRDS(spdata1, file = "~/Spatial_transcriptome/saveloca/sperate/sp-section-1.rds")
saveRDS(spdata2, file = "~/Spatial_transcriptome/saveloca/sperate/sp-section-2.rds")

#===============================================================================================
#===============================================================================================
# Data Merging==================================================================
## Column addition for grouping section separately
spdata1 <- readRDS("~/Spatial_transcriptome/saveloca/sperate/sp-section-1.rds")
spdata2 <- readRDS("~/Spatial_transcriptome/saveloca/sperate/sp-section-2.rds")
spdata1$Slice <- c('slice1')
spdata2$Slice <- c('slice2')
#===============================================================================================
#===============================================================================================
##    +Merging NO/without Batch-correction=============
### Merging= = = = = = 
merg.sections.nBC <- merge(spdata1, spdata2)

merg.nBC.tab <- data.frame(merg.sections.nBC@meta.data)

### Normalization- SCTransform  = = = = =
merg.sections@assays
merg.sections.nBC <- SCTransform(merg.sections.nBC, assay = "SCT", verbose = FALSE)

### Dimensionality reduction  = = = = = = =
merg.sections.nBC <- RunPCA(merg.sections.nBC, assay = "SCT", verbose = FALSE, ndims = 50)

### Elbow-plot = = = = = = = = = = = = = = 
png(file =  "png/mgNB/Elbow plot-Merg-data--nBC.png", width = 800, height = 500) 
ElbowPlot(merg.sections.nBC, ndims = 50)
dev.off()
### Find Neighbours  = = = = = = = = = = = = = =
merg.sections.nBC <- FindNeighbors(merg.sections.nBC, assay = 'SCT',reduction = 'pca', dims = 1:30)

### Find clusters = = = = = = = = = = = = = =
merg.sections.nBC <- FindClusters(merg.sections.nBC, verbose = FALSE, resolution = 0.5)

### UMAP = = = = = = = = = = = = = =
merg.sections.nBC <- RunUMAP(merg.sections.nBC, reduction = 'pca', dims =  1:30)

### Visualization UMAP and  spatial-dim-plot = = = = = = = =
png(file =  "png/mgNB/UMAP-merg-sections-group-by-SCT_snn_res.0.5--nBC.png", width = 900, height = 600)
DimPlot(merg.sections.nBC, reduction = 'umap', group.by = 'SCT_snn_res.0.5')+DarkTheme()
dev.off()

png(file =  "png/mgNB/UMAP-merg-sections-resolution_0.5-group-by-slice-nBC.png", width = 900, height = 600)
DimPlot(merg.sections.nBC, reduction = 'umap', group.by = c('Slice'),  label = TRUE)+DarkTheme()
dev.off()


png(file =  "png/mgNB/Sptial-Dim-plot-merg-sections--nBC.png", width = 1200, height = 750)
SpatialDimPlot(merg.sections.nBC, label = TRUE, label.size = 7, label.box = TRUE, label.color = "black")
dev.off()

# Saving the data without Batch correction = = = = = = = = = 
saveRDS(merg.sections.nBC, file = "~/Spatial_transcriptome/saveloca/merg.sections_nBC.rds")
# read 
merg.sections.nBC <- readRDS('~/Spatial_transcriptome/saveloca/merg.sections_nBC.rds')

#===============================================================================================
#===============================================================================================
##    +Merging with Batch-Correction method========================================
### Merging data
merg_wBCr <- merge(spdata1, spdata2)

### Splitting by slice column
merg_wBCr.list <- SplitObject(merg_wBCr, split.by = "Slice")

### Normalizing and finding variable features
merg_wBCr.list <- lapply(X = merg_wBCr.list, FUN = SCTransform, assay = 'SCT')

### integration of features that are variable across the data sets
feature.wBCr <- SelectIntegrationFeatures(object.list = merg_wBCr.list, nfeatures = 3000)

### Prepare SCT integration
merg_wBCr.list <- PrepSCTIntegration(object.list = merg_wBCr.list, anchor.features = feature.wBCr)

### Find integration anchors
anchor.wBCr <- FindIntegrationAnchors(object.list = merg_wBCr.list,
                                         normalization.method = "SCT",
                                         anchor.features = feature.wBCr)

### Integrated data assay creation
combined.wBCr <- IntegrateData(anchorset = anchor.wBCr, normalization.method = "SCT")

### performing integrated analysis
DefaultAssay(combined.wBCr) <- 'integrated'
combined.wBCr@active.assay
combi.tab <- data.frame(combined.wBCr@meta.data)

### Dimensionality Reduction- PCA 
combined.wBCr <- RunPCA(combined.wBCr, verbose =  FALSE ,ndims = 50)

### Elbow-plot 
png(file =  "png/wBC/Elbow plot-Merg-data--wBC.png", width = 800, height = 500) 
ElbowPlot(combined.wBCr, ndims = 50)
dev.off()

### Find Neighbours  = = = = = = = = = = = = = =
combined.wBCr <- FindNeighbors(combined.wBCr, reduction = 'pca', dims = 1:30)

### Find clusters = = = = = = = = = = = = = =
combined.wBCr <- FindClusters(combined.wBCr, verbose = FALSE, resolution = 0.5)

### UMAP = = = = = = = = = = = = = =
combined.wBCr <- RunUMAP(combined.wBCr, reduction = 'pca', dims =  1:30)

data.frame(combined.wBCr@meta.data)

### Visualization UMAP = = = = = = = =
png(file =  "png/wBC/UMAP-merg-sections-group-by-slices--wBC.png", width = 900, height = 600)
DimPlot(combined.wBCr, reduction = 'umap', group.by = c('Slice'), cols = c('yellow', 'red'))+DarkTheme()
dev.off()

png(file =  "png/wBC/UMAP-merg-sections-group-by-SCT_snn_res.0.5--wBC.png", width = 900, height = 600)
DimPlot(combined.wBCr, reduction = 'umap', group.by = 'SCT_snn_res.0.5')+DarkTheme()
dev.off()

png(file =  "png/wBC/UMAP-merg-sections-group-by-integrated_snn_res.0.5--wBC.png", width = 900, height = 600)
DimPlot(combined.wBCr, reduction = 'umap', group.by = 'integrated_snn_res.0.5')+DarkTheme()
dev.off()

# saving the integrated batch corrected data
saveRDS(combined.wBCr, file = "~/Spatial_transcriptome/saveloca/combined-data-with-BatchCorr.rds")

# Read 
combined.wBCr <- readRDS('~/Spatial_transcriptome/saveloca/combined-data-with-BatchCorr.rds')




#===============================================================================================
#===============================================================================================
# Cell-Type Identification======================================================
##   +Processing the Reference data============
### Read the reference data = = = = =
refrnce <- readRDS("~/Spatial_transcriptome/referncedata/allen_cortex.rds")

ref.tab <- data.frame(refrnce@meta.data)
### SC Transform  = = = = = = 
refrnce <- SCTransform(refrnce, assay = 'RNA', verbose =  FALSE)

### Dimensionality Reduction- PCA = = = = = = 
refrnce <- RunPCA(refrnce, ndims = 50, verbose = FALSE)

### Elbow plot= = = = = = =
png(file =  "png/ref/Elbow-plot-reference-data.png", width = 800, height = 500) 
ElbowPlot(refrnce, ndims = 50)
dev.off()

### Find Neighbors = = = = = = = =
refrnce <- FindNeighbors(refrnce, assay = 'RNA', reduction = 'pca', dims = 1:30)

### Find Clusters = = = = = = = =  
refrnce <- FindClusters(refrnce, verbose = FALSE, resolution = 0.5)

### Run UMAP = = = = = = = =  
refrnce <- RunUMAP(refrnce, reduction = 'pca', dims = 1:30)

### Visualize Reference = = = = = 
png(file =  "png/ref/UMAP-reference-subclass.png", width = 1400, height = 900)
DimPlot(refrnce, reduction = 'umap', group.by = 'subclass')+DarkTheme()
dev.off()


## Saving the processed reference data = = = = = 
saveRDS(refrnce, file = "~/Spatial_transcriptome/saveloca/reference.rds")

##   +Automatic Annotation using Data Integration=========
### Reference data set 
ref.cortex <- readRDS('~/Spatial_transcriptome/saveloca/reference.rds')

ref.tab <- data.frame(ref.cortex@meta.data)
ref.cortex@assays

## label projection processing with  batch corrected data as query
### - - -setting  batch corrected data  as query set 
query.wBC <- readRDS('~/Spatial_transcriptome/saveloca/combined-data-with-BatchCorr.rds')

anchor.wbc <- FindTransferAnchors(reference = ref.cortex, query = query.wBC, dims = 1:30)
prediction.wbc <- TransferData(anchorset = anchor.wbc , refdata = ref.cortex$subclass, dims = 1:30)
query.wBC <- AddMetaData(object = query.wBC, metadata = prediction.wbc)
query.wBC <- RunUMAP(query.wBC, reduction = 'pca', dims = 1:30)

query.tab2 <- data.frame(query.wBC@meta.data)
table(query.wBC$predicted.id)

### Visualization
png(file =  "png/cellid/Annotated-ID-UMAP-with-batch-correction.png", width = 1400, height = 900) 
DimPlot(query.wBC, group.by = c('predicted.id'), reduction = 'umap', 
        label = TRUE, label.size = 7, label.box = TRUE, 
        label.color = 'black', repel =  TRUE)+DarkTheme()
dev.off()

### Save 
saveRDS(query.wBC,'~/Spatial_transcriptome/saveloca/celltypeAnno/cellAnnotate-ID-wBC.rds')

#===============================================================================================
#===============================================================================================
##   +Manual Annotation = = = = = = = = = = = = = = = = = = ====
query.wBC.mrk <- FindAllMarkers(query.wBC, only.pos = TRUE, min.pct = 0.25, 
                                logfc.threshold = 0.25) 

head(query.wBC.mrk)


## Visualization 
png(file =  "png/mAnn/Spatial-feature-plot-QueryData-with-batch-corr.png", width = 1700, height = 1000) 
SpatialFeaturePlot(object = query.wBC, features = c('Klk6'), alpha =c(0.1,1),)
dev.off()

png(file =  "png/mAnn/UMAP-plot-QueryData-with-batch-corr.png", width = 1400, height = 900) 
FeaturePlot(object = query.wBC, features = c('Trf')) + DarkTheme() 
dev.off()

#===============================================================================================
#===============================================================================================
# Deconvolution===============================================================
## ## Integration for deconvolution  = = = = = = 
## creating a list of original data 
sp.merg.list <- list(slice1 = spdata1, slice2 = spdata2 )

## SCTransform
sp.merg.list <- lapply(sp.merg.list, SCTransform, assay = 'Spatial', method = 'poisson')

## Maxsize
options(future.globals.maxSize = 2000 * 1024^2)

### Creating the Feature list
featurs.list <- SelectIntegrationFeatures(object.list = sp.merg.list, nfeatures = 3000, verbose = FALSE)

### Prepare SCT integration
sp.merg.list <- PrepSCTIntegration(object.list = sp.merg.list, anchor.features = featurs.list)

### Find integration anchors
anchor.int <- FindIntegrationAnchors(object.list = sp.merg.list,
                                     normalization.method = "SCT",
                                     anchor.features = featurs.list,
                                     verbose = FALSE)

### Integrated data assay creation
combi.sp.list <- IntegrateData(anchorset = anchor.int, normalization.method = "SCT",
                               verbose = FALSE)



rm(anchor.int, sp.merg.list)
gc()


# Standard Processing
combi.sp.list <- RunPCA(combi.sp.list, verbose = FALSE)
combi.sp.list <- FindNeighbors(combi.sp.list, dims = 1:30)
combi.sp.list <- FindClusters(combi.sp.list, verbose = FALSE)
combi.sp.list <- RunUMAP(combi.sp.list, dims = 1:30)

DimPlot(combi.sp.list, reduction = 'umap', group.by = c('ident', 'orig.ident'))

SpatialDimPlot(combi.sp.list)

saveRDS(combi.sp.list, '~/Spatial_transcriptome/saveloca/deconvo/integrated-ST-wBC-deconvo.rds')


## prepare the reference data for Deconvolution= = = = = = = = = = = = = = = 
### Read the processed reference data file
ref.dcnvo <- readRDS('~/Spatial_transcriptome/saveloca/reference.rds')
integrat  <- readRDS('~/Spatial_transcriptome/saveloca/deconvo/integrated-ST-wBC-deconvo.rds') 


ref.dcnv.tab <-data.frame(ref.dcnvo@meta.data)
table(ref.dcnvo$subclass)
ref.dcnvo@active.assay = "RNA"

### to downsize the cell to 250 per subclass of the sample
### setting subclass as active indent
Idents(ref.dcnvo) <- ref.dcnvo$subclass
ref.dcnvo <- subset(ref.dcnvo, cells = WhichCells(ref.dcnvo, downsample = 250))

### SC transform 
ref.dcnvo <- SCTransform(ref.dcnvo, ncells = 3000, 
                         verbose = FALSE, method = 'poisson') %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)


png(file =  "png/deconvo/UMAP-plot-deconvolution-with-annotation-p.png", width = 1400, height = 900) 
DimPlot(ref.dcnvo, reduction = 'umap', label = TRUE, label.size = 7,  
        label.box = TRUE,  label.color = 'black', repel =  TRUE)+DarkTheme()
dev.off()



table(integrat$Slice)
table(integrat$seurat_clusters)

integrat@images$slice1
integrat@images$slice2.1

table(integrat$orig.ident)
### Subsetting the integrated spatial data 
ST.wBC <- subset(integrat, idents = c(0, 1, 2, 3, 4, 5))

integrat <- subset(integrat, slice1_imagerow > 400 | slice1_imagecol < 150, invert = TRUE)
integrat <- subset(integrat, slice1_imagerow > 275 & slice1_imagecol > 370, invert = TRUE)
integrat <- subset(integrat, slice1_imagerow > 250 & slice1_imagecol > 440, invert = TRUE)

integrat <- subset(integrat, slice2_imagerow > 400 | slice2_imagecol < 150, invert = TRUE)
integrat <- subset(integrat, slice2_imagerow > 275 & slice2_imagecol > 370, invert = TRUE)
integrat <- subset(integrat, slice2_imagerow > 250 & slice2_imagecol > 440, invert = TRUE)


## Select genes for deconvolution  = = = = = = = = = = = 
deconvo_marker <- FindAllMarkers(ref.dcnvo, only.pos = TRUE, logfc.threshold = 0.1,
                                 test.use = 'wilcox', min.pct = 0.05, min.diff.pct =  0.1, 
                                 max.cells.per.ident = 250) 

## Filtering the gene that are present in the ST data
deconvo_marker <- deconvo_marker[deconvo_marker$gene  %in%  rownames(integrat), ]


### selecting 20 differential expressed genes per cell subclass
deconvo_marker$pct.diff <- deconvo_marker$pct.1 - deconvo_marker$pct.2

deconvo_marker$log.pct.diff <- log2((deconvo_marker$pct.1 * 99 + 1) / (deconvo_marker$pct.2 * 99 + 1))

deconvo_marker %>% 
  group_by(cluster) %>% 
  top_n(-100, p_val) %>% 
  top_n(50, pct.diff) %>% 
  top_n(20, log.pct.diff) -> top20_DEG


mrk.deco.uq <- unique(as.character(top20_DEG$gene))



## Create an Expression-set-Object = = = = = = = = = = =
exprs.set.SC <- ExpressionSet(assayData = as.matrix(ref.dcnvo@assays$RNA@counts[mrk.deco.uq,]),
                              phenoData = AnnotatedDataFrame(ref.dcnvo@meta.data))

exprs.set.ST <- ExpressionSet(assayData = as.matrix(integrat@assays$Spatial@counts[mrk.deco.uq,]),
                              phenoData = AnnotatedDataFrame(integrat@meta.data))
  


## Deconvolution= = = = = = = = = = = = = = = = = = = 

Deconvolution <- SCDC::SCDC_prop(bulk.eset = exprs.set.ST,sc.eset = exprs.set.SC, 
                                 ct.varname = 'subclass',
                                  ct.sub = as.character(unique(exprs.set.SC$subclass)))

head(Deconvolution$prop.est.mvw)

### Adding the deconvolution results and add it to Seurat object as new assay- - - 


integrat@assays[["SCDC"]] <- CreateAssayObject(data = t(Deconvolution$prop.est.mvw))

if (length(integrat@assays$SCDC@key) == 0){
   integrat@assays$SCDC@key = 'scdc_'
   }

DefaultAssay(integrat) <- 'SCDC'

## Spatial Feature plot
png(file =  "png/deconvo/SpatialFeaturePlot-Deconvolution.png", width = 1400, height = 900) 
SpatialFeaturePlot(integrat, features =  c('L5 PT','Macrophage'), pt.size.factor = 1.6, crop = TRUE)
dev.off()

integrat <- FindSpatiallyVariableFeatures(integrat, assay = "SCDC", selection.method = "markvariogram",
                                         features = rownames(integrat), r.metric = 5, slot = 'data')

top.cluster.deconvo <- head(SpatiallyVariableFeatures(integrat),3)

png(file =  "png/deconvo/Top-Cluster-Feature-Deconvolution.png", width = 1400, height = 900) 
SpatialPlot(object = integrat, features = top.cluster.deconvo, ncol = 2)
dev.off()

#===============================================================================================
#===============================================================================================





























