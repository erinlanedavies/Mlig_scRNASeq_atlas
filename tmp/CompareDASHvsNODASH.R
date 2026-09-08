

### Comparing DASH versus non-DASH

# SCT
DefaultAssay(seurat_integrated_final) <- "SCT"
seurat_integrated_final <- PrepSCTFindMarkers(seurat_integrated_final)
DEGS_Ad_non_DASH_vs_DASH_SCT <- FindMarkers(seurat_integrated_final, ident.1 = "2_Adult_Unsorted", ident.2 = "2_Adult_Unsorted_DASH")


# RNA
DefaultAssay(seurat_integrated_final) <- "RNA"
seurat_integrated_final <- NormalizeData(seurat_integrated_final, normalization.method = "LogNormalize", scale.factor = 10000)
DEGS_Ad_non_DASH_vs_DASH_RNA <- FindMarkers(seurat_integrated_final, ident.1 = "2_Adult_Unsorted", ident.2 = "2_Adult_Unsorted_DASH")



