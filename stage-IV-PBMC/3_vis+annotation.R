rm(list = ls())

srPostPro = function(srsc, plotPf, res = 600, ndim = 15, rdsSave = FALSE, plotFeat = NULL, metaPlot = NULL, metaFPlot = NULL, adtPlot = FALSE, compCols = NULL,
                     jsFlag = FALSE, tsneFlag = FALSE, logFCThr = 0.25, qcFlag = FALSE, bm1 = NULL, bm2 = NULL) {
  cat("Start to QC and pre-processing Seurat object...\n")
  
  if (qcFlag) {
    cat("\tVisualize the basic QC parameters for integrated object...\n")
    qcgg <- VlnPlot(srsc, features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb", "percent_plat"), ncol = 3)
    ggsave(plot = qcgg, filename = paste(plotPf, "QCPlot1.png", sep = ""), dpi = res, width = 18, height = 18)
  }
  ##############
  srsc <- ScaleData(srsc, assay = "integrated")
  srsc <- ScaleData(srsc, assay = "integrated_adt")
  #############
  srsc <- RunPCA(srsc, features = VariableFeatures(object = srsc))
  
  pcaGG1 <- VizDimLoadings(srsc, dim = 1:6, reduction = "pca")
  ggsave(plot = pcaGG1, filename = paste(plotPf, "PCA1.png", sep = ""), dpi = res, width = 9, height = 16)
  
  pcaGG2 <- DimPlot(srsc, reduction = "pca")
  ggsave(plot = pcaGG2, filename = paste(plotPf, "PCA2.png", sep = ""), dpi = res, width = 9, height = 6)
  
  png(paste(plotPf, "PCA_HEATMAP.png", sep = ""), res = res, width = 9, height = 16, units = "in")
  DimHeatmap(srsc, dims = 1:15, cells = 500, balanced = TRUE)
  gar = dev.off()
  
  if (jsFlag) {
    jst <- Sys.time()
    srsc <- JackStraw(srsc, num.replicate = 100)
    srsc <- ScoreJackStraw(srsc, dims = 1:20)
    jsGG <- JackStrawPlot(srsc, dims = 1:20)
    ggsave(plot = jsGG, filename = paste(plotPf, "JackStraw.png", sep = ""), dpi = res, width = 9, heigh = 6)
    cat("Jack Straw test cost:")
    print(Sys.time()-jst)
  }
  elGG <- ElbowPlot(srsc, ndims = 50)
  ggsave(plot = elGG, filename = paste(plotPf, "Elbow.png", sep = ""), dpi = res, width = 9, heigh = 6)
  
  srsc <- RunUMAP(srsc, dims = 1:10)
  srsc <- FindNeighbors(srsc, dims = 1:10)
  srsc <- FindClusters(srsc, resolution = 0.1)
  
  umapGG <- DimPlot(srsc, reduction = "umap", label = TRUE)
  ggsave(plot = umapGG, filename = paste(plotPf, "UMAP.png", sep = ""), dpi = res, width = 9, heigh = 6)
  
  if (tsneFlag) {
    srsc = RunTSNE(srsc, dims = 1:ndim)
    tsneGG = DimPlot(srsc, reduction = "tsne", label = TRUE)
    ggsave(plot = tsneGG, filename = paste(plotPf, "TSNE.png", sep = ""), dpi = res, width = 9, heigh = 6)
  }
  
  outData <- FetchData(srsc, c("ident", "nGene", "UMAP1", "UMAP2"))
  write.csv(outData, file = paste(plotPf, "useful_results.csv", sep = ""))
  
  if (adtPlot) {
    allAdt <- rownames(srsc[["integrated_adt"]])
    for (iadt in allAdt) {
      cat("\t", iadt, "\n")
      tmpVln <- VlnPlot(srsc, iadt, assay = "integrated_adt")
      ggsave(plot = tmpVln, filename = paste(plotPf, iadt, "_violin.png", sep = ""), 
             dpi = res, width = 12, heigh = 6)
      tmpFeat <- FeaturePlot(srsc, features = iadt, order=TRUE) #, assay = "integrated_adt")
      ggsave(plot = tmpFeat, filename = paste(plotPf, iadt, "_feature_umap.png", sep = ""), 
             dpi = res, width = 9, heigh = 6)
      
      irna <- str_replace_all(iadt, "\\-Ab", "")
      irna <- str_replace_all(irna, "\\-Isotype", "")
      irna <- toupper(irna)
      
      if (irna %in% rownames(srsc[['RNA']])) {
        tmpFeatSct <- FeatureScatter(srsc, feature1 = iadt, feature2 = irna)+ theme_bw()
        ggsave(plot = tmpFeatSct, filename = paste(plotPf, iadt, "_feature_scatter.png", sep = ""), 
               dpi = res, width = 7.5, heigh = 6)
      }
    }
  }
  
  if (!is.null(plotFeat)) {
    cat("Start to plot feature plots and violin plots\n")
    plotFeat = paste0("rna_", plotFeat)
    vlnGG = VlnPlot(srsc, features = plotFeat, assay = "RNA", log = TRUE, ncol = 3, pt.size = 0.5)
    ggsave(plot = vlnGG, filename = paste(plotPf, "violin_plot_for_cell_annotation.png", sep = ""), 
           dpi = res, width = 32, heigh = 40)
    featGG = FeaturePlot(srsc, features = plotFeat, ncol = 3, pt.size = 0.2, order = TRUE)
    ggsave(plot = featGG, filename = paste(plotPf, "umap_feature_plot_for_cell_annotation.png", sep = ""), 
           dpi = res, width = 20, heigh = 32)
    dotGG = DotPlot(srsc, features = plotFeat, assay = "RNA", cols = c("#6600CC", "gold")) + 
      theme_bw()+
      RotatedAxis() + 
      coord_flip()
    ggsave(plot = dotGG, filename = paste(plotPf, "ident_dotplot_gene_list.png", sep = ""), 
           dpi = res, width = 9, heigh = 15)
    
  }
  
  if (!is.null(metaPlot)) {
    cat("Start to plot reduction dim plots\n")
    for (imp in metaPlot) {
      print(imp)
      dimGG = DimPlot(srsc, group.by = imp, pt.size = 0.2, order = imp)
      ggsave(plot = dimGG, filename = paste(plotPf, "umap_dim_plot_", imp, ".png", sep = ""), 
             dpi = res, width = 9, heigh = 6)
      
      splDimGG = DimPlot(srsc, group.by = imp, split.by = imp, pt.size = 0.2, order = imp)
      ggsave(plot = splDimGG, filename = paste(plotPf, "umap_split_dim_plot_", imp, ".png", sep = ""), 
             dpi = res, width = 18, heigh = 6)
      
      splDimGG = DimPlot(srsc, group.by = "ident", split.by = imp, pt.size = 0.2, order = imp)
      ggsave(plot = splDimGG, filename = paste(plotPf, "umap_double_split_dim_plot_", imp, ".png", sep = ""), 
             dpi = res, width = 18, heigh = 6)
      
      tmpCts = as.data.frame.matrix(table(Idents(srsc), srsc@meta.data[,imp]))
      write.csv(tmpCts, paste(plotPf, "cell_counts_cross_idents_", imp, ".csv", sep = ""))
    }
  }
  
  if (!is.null(metaFPlot)) {
    cat("Start to plot reduction feature plots\n")
    for (imp in metaFPlot) {
      print(imp)
      dimGG = FeaturePlot(srsc, features = imp, pt.size = 0.2, order = TRUE)
      ggsave(plot = dimGG, filename = paste(plotPf, "umap_feature_plot_", imp, ".png", sep = ""), 
             dpi = res, width = 9, heigh = 6)
      
      splDimGG = DimPlot(srsc, features = imp, split.by = 'seurat_clusters', pt.size = 0.2, order = TRUE)
      ggsave(plot = splDimGG, filename = paste(plotPf, "umap_double_split_feature_plot_", imp, ".png", sep = ""), 
             dpi = res, width = 18, heigh = 6)
    }
  }
  
  if (length(levels(outData$ident)) >= 2) {
    srsc.markers = FindAllMarkers(srsc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = logFCThr)
    topMarkers = srsc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    print(head(topMarkers))
    write.csv(srsc.markers, file = paste(plotPf, "ALL_POS_MARKERS.csv", sep = ""))
    write.csv(topMarkers, file = paste(plotPf, "top_pos_markers.csv", sep = ""))
    png(paste(plotPf, "marker_heatmap.png", sep = ""), res = res, width = 32, height = 24, units = "in")
    print(DoHeatmap(object = srsc, features = topMarkers$gene))
    gar = dev.off()
    dotGG = DotPlot(srsc, features = unique(topMarkers$gene), assay = "RNA", cols = c("#6600CC", "gold")) +
      theme_bw()+
      RotatedAxis() + 
      coord_flip()
    ggsave(plot = dotGG, filename = paste(plotPf, "top_marker_dotplot.png", sep = ""), 
           dpi = res, width = 9, heigh = 30)#原来height = 21
    
    srsc.markers = FindAllMarkers(srsc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = logFCThr, assay = "integrated_adt")
    topMarkers = srsc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
    print(head(topMarkers))
    write.csv(srsc.markers, file = paste(plotPf, "ALL_POS_ADT_MARKERS.csv", sep = ""))
    write.csv(topMarkers, file = paste(plotPf, "top_pos_adt_markers.csv", sep = ""))
    png(paste(plotPf, "adt_marker_heatmap.png", sep = ""), res = res, width = 32, height = 24, units = "in")
    print(DoHeatmap(object = srsc, features = topMarkers$gene, assay = "integrated_adt"))
    gar = dev.off()
    dotGG = DotPlot(srsc, features = unique(topMarkers$gene), assay = "integrated_adt", cols = c("#6600CC", "gold")) + 
      theme_bw()+
      RotatedAxis() + 
      coord_flip()
    ggsave(plot = dotGG, filename = paste(plotPf, "top_adt_pos_marker_dotplot.png", sep = ""), 
           dpi = res, width = 9, heigh = 16)
    srsc.markers = FindAllMarkers(srsc, only.pos = FALSE, min.pct = 0.25, logfc.threshold = logFCThr, assay = "integrated_adt")
    write.csv(srsc.markers, file = paste(plotPf, "ALL_ADT_MARKERS.csv", sep = ""))
  }
  if (rdsSave) {
    saveRDS(srsc, file = paste(plotPf, "Seurat_Objects_Clustered.RDS", sep = ""))
  }
  print(srsc)
  return(srsc)
}

suppressMessages(library(DoubletFinder))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratData))
suppressMessages(library(SeuratDisk))
suppressMessages(library(patchwork))
suppressMessages(library(forcats))
suppressMessages(library(showtext))
suppressMessages(library(SCP))

inputDir <- "6CloupeFile/results/"

sigMarker <- c("PTPRC", "EPCAM", "COL1A1", "CD3D", "CD8A", "CD4", "PDCD1", "TCF7", "NKG7", "MS4A1", "CD27", "CD19", 
               "CD68", "FCER1A", "HLA-DRA", "NCAM1", "TPSAB1", "GNLY", "FOXP3", "TBX21", "GATA3","CD33")
mycolor <- c("#f9d423","#fc913a","#6a60a9")#黄橙紫

expID <- "LCCHEMOICBPBMC_4_SINGLET_221202"
filePat <- "singlet.rds"


clusterFlag <- TRUE
rnaInteMethod <- "RPCA"
adtInteMethod <- "RPCA"

dir.create(file.path(inputDir, expID), showWarnings = FALSE)
expDir <- paste(inputDir, expID, sep = "/")
ppf <- paste(expDir, "/", expID, "_", sep = "")


inteObjs <- readRDS(paste(ppf, "adt_", adtInteMethod,  "_rna_", rnaInteMethod, "_integrated.rds", sep = ""))
inteObjs@meta.data$pre_post <- str_split_fixed(inteObjs@meta.data$orig.ident, "_", n = 3)[,1]
inteObjs@meta.data$pid <- str_split_fixed(inteObjs@meta.data$orig.ident, "_", n = 3)[,2]
print(head(rownames(inteObjs[['RNA']])))
print(inteObjs)

if (clusterFlag) {
  clstFolder <- paste(expID, "cluster", sep = "_")
  dir.create(file.path(expDir, clstFolder), showWarnings = FALSE)
  clstDir <- paste(expDir, clstFolder, sep = "/")
  tmp_pf <- paste(clstDir, "/", clstFolder, "_", sep = "")
  
  inteProObj <- srPostPro(inteObjs, plotPf = tmp_pf, rdsSave = TRUE, ndim = 30, jsFlag = TRUE, qcFlag = TRUE, adtPlot = TRUE,
                          plotFeat = sigMarker, metaPlot = c("orig.ident", "pre_post", "pid"))
}

#######supplementary figure-------------
intgene <- c("CD33","CD14","CD163","CD3D","CD8A","CD4","NCAM1","FCGR3A","NKG7","CD19","CD79A","MS4A1","PTPRC")
umapGG = FeatureDimPlot(srt = inteProObj, features = intgene,theme_use = "theme_blank",assay = "RNA",ncol = 3)
ggsave(plot = umapGG, filename = paste(tmp_pf, "UMAP_intgene_RNA.png", sep = ""), dpi = 600, width = 12, heigh = 20,units = "in")
ggsave(plot = umapGG, filename = paste(tmp_pf, "UMAP_intgene_RNA.pdf", sep = ""), width = 12, heigh = 20,units = "in")

use_df <- inteProObj@meta.data
use_df <- use_df %>%
  mutate(anno_gzk = case_when(
    seurat_clusters %in% c(0,1,2,3,4,7,11,12,13,14,16,18,20,21,22,23) ~ "Myeloid",
    seurat_clusters %in% c(5,10) ~ "T cell", 
    seurat_clusters %in% c(6,15,24) ~ "B cell",
    seurat_clusters %in% c(8,17,19) ~ "Other",
    seurat_clusters %in% c(9) ~ "NK",
    TRUE ~ "unknown"  
  ))
table(use_df$anno_gzk)
inteProObj <- AddMetaData(inteProObj,metadata = use_df)

umapGG = CellDimPlot(inteProObj, group.by = "anno_gzk",reduction = "UMAP",theme_use = "theme_scp",label = T)
ggsave(filename = paste(tmp_pf,"UMAP_2.png",sep = ""),umapGG,dpi = 600, width = 8, heigh = 6,units = "in")

###Fig.S1
umapGG = CellDimPlot(inteProObj, group.by = "anno_gzk",reduction = "UMAP",theme_use = "theme_blank",label = F)
ggsave(filename = paste(tmp_pf,"UMAP_final.png",sep = ""),umapGG,dpi = 600, width = 6, heigh = 4,units = "in")
ggsave(filename = paste(tmp_pf,"UMAP_final.pdf",sep = ""),umapGG,width = 6, heigh = 4,units = "in")
write.csv(use_df,file = paste(tmp_pf,"anno_gzk_all_cell.csv",sep = ""),row.names = T)


Myeloid_genes <- c("CD33","CD14","CD163")#FCGR3A
T_genes <- c("CD3D","CD3E","CD8A")
NK_genes <- c("NCAM1",'FCGR3A',"NKG7")
B_genes <- c("CD19","CD79A","MS4A1")#"MX1"
Other_genes <- c("PTPRC")

features <- list("Myeloid" = Myeloid_genes,
                 "T" = T_genes,
                 "NK" = NK_genes,
                 "B" = B_genes,
                 "Other" = Other_genes)

dotGG <- DotPlot(object = inteProObj, features = features, assay = "RNA", group.by = "anno_gzk",
                 scale.by = "size",scale.min = 0,scale.max = 100) +
  geom_point(aes(size = pct.exp),shape  = 21,stroke = 0.5)+
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", family = "sans"),
        axis.text.y = element_text(color = "black", family = "sans"),
        panel.grid.major = element_line(colour = "#ebebeb",
                                        linetype = "dashed"), 
        panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
        panel.background = element_rect(linetype = "dashed"),
        plot.background = element_rect(colour = NA),
        legend.key = element_rect(fill = NA),
        legend.background = element_rect(fill = NA))+
  labs(x = NULL, y = NULL, colour = NULL,size = NULL)+
  scale_color_distiller(palette = "Spectral") +
  RotatedAxis() +
  scale_alpha_continuous(name = "Per Exp.") +
  guides(color = guide_colorbar(title = c("Avg Exp.")))+
  scale_y_discrete(limits = c("Other", "B cell", "NK", "T cell", "Myeloid"))
print(dotGG)
ggsave(plot = dotGG, filename = paste0(tmp_pf,"dot_ALL_CELL_2.png"),width = 8, heigh = 3,units = "in",dpi = 600)
ggsave(plot = dotGG, filename = paste0(tmp_pf,"dot_ALL_CELL_2.pdf"),width = 8, heigh = 3,units = "in")

####save Myeloid cell data------------
Myeloid <- subset(inteProObj,subset = anno_gzk%in%c("Myeloid"))
saveRDS(Myeloid,file = "LCCHEMOICBPBMC_4_SINGLET_221202_cluster_230110_Seurat_Objects_Clustered_anno.RDS")
