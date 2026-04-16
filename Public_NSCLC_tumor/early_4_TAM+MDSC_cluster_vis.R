# only subset TAM and MDSC in CD33 data to recluster------------
rm(list = ls())

srPostPro = function(srsc, plotPf, res = 600, ndim = 15, rdsSave = FALSE, plotFeat = NULL, metaPlot = NULL, metaFPlot = NULL, adtPlot = FALSE, compCols = NULL,
                     jsFlag = FALSE, tsneFlag = FALSE, logFCThr = 0.25, qcFlag = FALSE, bm1 = NULL, bm2 = NULL) {
  cat("Start to QC and pre-processing Seurat object...\n")
  
  if (qcFlag) {
    cat("\tVisualize the basic QC parameters for integrated object...\n")
    qcgg <- VlnPlot(srsc, features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb", "percent_plat"), ncol = 3)
    ggsave(plot = qcgg, filename = paste(plotPf, "QCPlot1.png", sep = ""), dpi = res, width = 18, height = 18)
    
    srsc <- subset(srsc, subset = nFeature_RNA > 200 & percent_mito < 25 & percent_ribo > 5 & percent_hb < 5)
    qcgg <- VlnPlot(srsc, features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb", "percent_plat"), ncol = 3)
    ggsave(plot = qcgg, filename = paste(plotPf, "after_QCPlot1.png", sep = ""), dpi = res, width = 18, height = 18)
    
    
  }
  ##############
  srsc <- ScaleData(srsc, assay = "integrated",verbose = FALSE)
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
  
  srsc <- RunUMAP(srsc, dims = 1:ndim)
  srsc <- FindNeighbors(srsc, dims = 1:ndim)
  srsc <- FindClusters(srsc, resolution = 0.5)
  
  umapGG <- DimPlot(srsc, reduction = "umap", label = TRUE)
  ggsave(plot = umapGG, filename = paste(plotPf, "UMAP.png", sep = ""), dpi = res, width = 9, heigh = 6)
  
  if (tsneFlag) {
    srsc = RunTSNE(srsc, dims = 1:ndim)
    tsneGG = DimPlot(srsc, reduction = "tsne", label = TRUE)
    ggsave(plot = tsneGG, filename = paste(plotPf, "TSNE.png", sep = ""), dpi = res, width = 9, heigh = 6)
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
    dotGG = DotPlot(srsc, features = plotFeat, assay = "RNA") +
      theme_classic()+
      theme(axis.text.x = element_text(color = "black",family = "sans"),
            axis.text.y = element_text(color = "black",family = "sans"))+
      scale_color_distiller(palette = "RdBu")+
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
  
  if (rdsSave) {
    saveRDS(srsc, file = paste(plotPf, "Seurat_Objects_Clustered.RDS", sep = ""))
  }
  print(srsc)
  return(srsc)
}

suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratData))
suppressMessages(library(SeuratDisk))
suppressMessages(library(SCP))
suppressMessages(library(RColorBrewer))
suppressMessages(library(MAST))
suppressMessages(library(clustree))

set.seed(1)
inputDir <- "results/"

sigMarker <- c("PTPRC", "EPCAM", "COL1A1", "CD3D", "CD8A", "CD4", "PDCD1", "TCF7", "NKG7", "MS4A1", "CD27", "CD19", 
               "CD68", "FCER1A", "HLA-DRA", "NCAM1", "TPSAB1", "GNLY", "FOXP3", "TBX21", "GATA3","CD33","CCL5")

expID <- "LCNAIVETUMOR_49_SINGLET_241115"
filePat <- "singlet.rds"

clusterFlag <- TRUE
rnaInteMethod <- "RPCA"
adtInteMethod <- "RPCA"

dir.create(file.path(inputDir, expID), showWarnings = FALSE)
expDir <- paste(inputDir, expID, sep = "/")
ppf <- paste(expDir, "/", expID, "_", sep = "")

inteObjs <- readRDS(file = "LCNAIVETUMOR_49_SINGLET_241115_cluster_CD33_Annotation_251017.rds")
inteObjs <- subset(inteObjs,subset = cell_anno%in%c("Mph_APOE","Mph_FBP1","Mono_FCN1","MDSC-like","Mph_C1Q"))
inteObjs@meta.data <- inteObjs@meta.data[,c(1:14,23:24)]
print(inteObjs)
DefaultAssay(inteObjs) <- "integrated"

if (clusterFlag) {
  clstFolder <- paste(expID, "cluster_CD33_sub_TAM_rmMono", sep = "_")
  dir.create(file.path(expDir, clstFolder), showWarnings = FALSE)
  clstDir <- paste(expDir, clstFolder, sep = "/")
  tmp_pf <- paste(clstDir, "/", clstFolder, "_", sep = "")
  
  inteProObj <- srPostPro(inteObjs, plotPf = tmp_pf, rdsSave = TRUE, ndim = 30, jsFlag = TRUE, qcFlag = TRUE, adtPlot = TRUE,
                          plotFeat = sigMarker, metaPlot = c("disease"))
}

inteProObj <- readRDS(file = paste(tmp_pf, "Seurat_Objects_Clustered.RDS", sep = ""))

ress <- c(seq(0.1,0.2, 0.1))
for (iress in ress) {
  inteProObj <- FindClusters(inteProObj, resolution = iress)
}

cluster <- paste0("integrated_snn_res.",ress)
for (icluster in cluster) {
  umapGG = CellDimPlot(inteProObj, group.by = icluster, reduction = "UMAP", label = TRUE, label_insitu = TRUE)
  ggsave(plot = umapGG, filename = paste(tmp_pf, "UMAP_",icluster,".png", sep = ""), dpi = 300, width = 8, heigh = 6,units = "in")
}

inteProObj$cell_anno2 = "NA"
celltype = c("Mph_APOE","Mph_IL1B","Mph_SPP1","Mph_MARCO","MDSC-like","TAM_TNF","Mph_monocyte","MDSC-like")

sub_length = length(unique(inteProObj$integrated_snn_res.0.2)) - 1
for (i in 0:sub_length){ 
  inteProObj$cell_anno2[inteProObj$integrated_snn_res.0.2==i] = celltype[i+1]
}

inteProObj$cell_anno2 <- factor(inteProObj$cell_anno2,levels = c("MDSC-like","Mph_IL1B","Mph_SPP1","Mph_monocyte","Mph_APOE","Mph_MARCO","TAM_TNF"))
color_use <- c("#d04c50","#57a258","#a4c877","#e5d6e1","#F17F42","#90cbfb","#f9d423")

###Fig.3C-------------------
umapGG = CellDimPlot(inteProObj, group.by = "cell_anno2", reduction = "UMAP", label = F, label_insitu = TRUE,palcolor = color_use)
ggsave(plot = umapGG,filename = paste(tmp_pf,"TAM_anno_2.png", sep = ""), dpi = 600, width = 8, heigh = 6,units = "in")
ggsave(plot = umapGG,filename = paste(tmp_pf,"TAM_anno_2.pdf", sep = ""), width = 8, heigh = 6,units = "in")

#saveRDS(inteProObj,file = "LCNAIVETUMOR_49_SINGLET_241115_cluster_CD33_sub_TAM_rmMono_annotation_260130.RDS")
#inteProObj <- readRDS(file = "LCNAIVETUMOR_49_SINGLET_241115_cluster_CD33_sub_TAM_rmMono_annotation_260130.RDS")

use_gene <- c("CCL5","NKG7","GZMA","IL32","CD3E",#MDSC
              "CD14","FCGR3A",
              # "JCHAIN","IGHA1","CD79A","CD19","MS4A1",#B-like cluster
              # "PLD4","IL3RA","LILRA4","CLEC4C",#pDC-like cluster
              "IL1B","AREG","EREG","TIMP1",#IL1B
              "SPP1","CTSB","CTSL","FABP5",#SPP1
              "VCAN","S100A8","S100A9","LYZ",#monocyte
              "APOE","APOC1","SEPP1","LGMN",#TAM-APOE
              "MARCO","MRC1","FBP1","NUPR1",#MARCO cluster
              "PLD4","AXL","TNF","RGS1")#TAM_TNF

###Fig.3D-------------------
ht1 <- GroupHeatmap(inteProObj,features = use_gene,assay = "RNA",heatmap_palette = "viridis",group.by = "cell_anno2",
                    show_row_names = T)
ht1$plot
ggsave(filename = paste(tmp_pf, "heatmap_select_gene_USE_cell_anno2_251124.png", sep = ""), dpi = 600, width = 5, heigh = 5,units = "in")
ggsave(filename = paste(tmp_pf, "heatmap_select_gene_USE_cell_anno2_251124.pdf", sep = ""), width = 5, heigh = 5,units = "in")
