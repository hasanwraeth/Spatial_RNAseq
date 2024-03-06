#Library---------------
library(dplyr)
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(SeuratData)
library(ggplot2)
library(cowplot)
library(patchwork)
library(metap)
library(presto)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(velocyto.R)
library(pagoda2)
library(SeuratWrappers)
library(monocle3)
library(CellChat)
library(viridis)

#10X Visium------------------
#InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")

plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))

# plot <- SpatialFeaturePlot(brain, features = c("Ttr")) + theme(legend.text = element_text(size = 0),
#                                                                legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"))
# jpeg(filename = "../spatial_vignette_ttr.jpg", height = 700, width = 1200, quality = 50)
# print(plot)
# dev.off()

p1 <- SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features = "Ttr", alpha = c(0.1, 1))
p1 + p2

brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2

SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3,
                                                                                     5, 8)), facet.highlight = TRUE, ncol = 3)
