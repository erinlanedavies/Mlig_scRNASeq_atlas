
## Generating GO enrichment analysis for Sam's talk

library(topGO)

seurat_integrated_final <- PrepSCTFindMarkers(seurat_integrated_final)

Idents(seurat_integrated_final) <- "integrated_snn_res.1.2"

marker.genes <- FindAllMarkers(seurat_integrated_final, assay = "SCT")

all.genes <- rownames(seurat_integrated_final@assays$SCT)

geneID2GOs <- read.table("/data/Mlig_scRNA_Seq/Francisco/data/geneid2go.map", sep = "\t", header = FALSE)
gene2GO <- strsplit(geneID2GOs$V2, ", ",fixed=T)
names(gene2GO) <- geneID2GOs$V1

results_enrichment <- data.frame(GO.ID=character(),
                                 Term=character(),
                                 Annotated=integer(),
                                 Significant=integer(),
                                 Expected =double(),
                                 Fisher=double(),
                                 Fisher_FDR=double(),
                                 ClusterID=character()
)

for (cluster in unique(Idents(seurat_integrated_final))) {
  print(cluster)
  genes <- marker.genes$gene[marker.genes$cluster == cluster & marker.genes$p_val_adj < 0.00001 & marker.genes$avg_log2FC > 2]
  # define geneList as 1 if gene is in expressed.genes, 0 otherwise
  geneList <- ifelse(all.genes %in% genes, 1, 0)
  names(geneList) <- all.genes
  # Create topGOdata object
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                geneSelectionFun = function(x)(x == 1),
                annot = annFUN.gene2GO,
                nodeSize = 3,
                gene2GO = gene2GO)
  
  # Test for enrichment using Fisher's Exact Test
  resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  results <- GenTable(GOdata, Fisher = resultFisher, topNodes = length(GOdata@graph@nodes), numChar = 1000)
  results$Fisher <- ifelse(results$Fisher == "< 1e-30", 1e-30, as.numeric(results$Fisher))
  results$Fisher_FDR <- p.adjust(results$Fisher)
  results$ClusterID <- cluster
  print(results[results$Fisher_FDR < 0.2,])
  results_enrichment <- rbind(results_enrichment, results)
}

openxlsx::write.xlsx(x = results_enrichment, file = "../results/GO_enrichment_analysis_Sam_talk.xlsx")

showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 20, useInfo = 'all')
