library(Seurat)

rownames(counts) = counts$tiles
counts = counts[,2:229]

counts = as.matrix(count)

seu = CreateSeuratObject(counts,scale.data = counts)
seu@assays$RNA@scale.data = counts

seu = FindVariableFeatures(seu)

seu <- RunPCA(seu,features = rownames(seu))

seu <- FindNeighbors(seu, reduction = "pca", dims = 1:50)
seu <- FindClusters(seu, resolution = 0.8, by = 0.1)
counts <- RunUMAP(counts, reduction = "pca", dims = 1:50)
DimPlot(counts,reduction = 'umap')
