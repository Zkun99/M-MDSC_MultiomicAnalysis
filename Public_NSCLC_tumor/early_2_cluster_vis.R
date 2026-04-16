###All cell cluster and annotation

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


suppressMessages(library(DoubletFinder))
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
               "CD68", "FCER1A", "HLA-DRA", "NCAM1", "TPSAB1", "GNLY", "FOXP3", "TBX21", "GATA3","CD33")

expID <- "LCNAIVETUMOR_49_SINGLET_241115"
filePat <- "singlet.rds"

clusterFlag <- TRUE
rnaInteMethod <- "RPCA"
adtInteMethod <- "RPCA"

dir.create(file.path(inputDir, expID), showWarnings = FALSE)
expDir <- paste(inputDir, expID, sep = "/")
ppf <- paste(expDir, "/", expID, "_", sep = "")

inteObjs <- readRDS(paste(ppf, "_rna_", rnaInteMethod, "_integrated.rds", sep = ""))
inteObjs@meta.data <- inteObjs@meta.data[,c(1:10,14)]
inteObjs@meta.data$tissue <- str_split_fixed(inteObjs@meta.data$sample_name, "_", n = 3)[,2]
inteObjs@meta.data$pid <- str_split_fixed(inteObjs@meta.data$sample_name, "_", n = 3)[,1]
inteObjs@meta.data$disease <- str_split_fixed(inteObjs@meta.data$sample_name, "_", n = 4)[,3]
print(head(rownames(inteObjs[['RNA']])))
print(inteObjs)


if (clusterFlag) {
  clstFolder <- paste(expID, "cluster", sep = "_")
  dir.create(file.path(expDir, clstFolder), showWarnings = FALSE)
  clstDir <- paste(expDir, clstFolder, sep = "/")
  tmp_pf <- paste(clstDir, "/", clstFolder, "_", sep = "")
  
  inteProObj <- srPostPro(inteObjs, plotPf = tmp_pf, rdsSave = TRUE, ndim = 30, jsFlag = TRUE, qcFlag = TRUE, adtPlot = TRUE,
                          plotFeat = sigMarker, metaPlot = c("pid","disease"))
}

inteProObj <- readRDS(file = paste(tmp_pf,"Seurat_Objects_Clustered.RDS", sep = ""))

use_df <- inteProObj@meta.data
use_df <- use_df %>%
  mutate(anno_2 = case_when(
    integrated_snn_res.0.5 %in% c(15) ~ "Epithelial",
    integrated_snn_res.0.5 %in% c(1,6,9,11,14) ~ "Myeloid",
    integrated_snn_res.0.5 %in% c(4) ~ "Plasma",
    integrated_snn_res.0.5 %in% c(0,3,5,7,8,20) ~ "T cell",
    integrated_snn_res.0.5 %in% c(10,16) ~ "NK",
    integrated_snn_res.0.5 %in% c(19) ~ "Fibroblast",
    integrated_snn_res.0.5 %in% c(2,18,12) ~ "B cell",
    integrated_snn_res.0.5 %in% c(17) ~ "Endothelial",
    integrated_snn_res.0.5 %in% c(13) ~ "Cycling cell",
    integrated_snn_res.0.5 %in% c(21) ~ "pDC",
    TRUE ~ "unknown"  # 相当于else
  ))
table(use_df$anno_2)
inteProObj <- AddMetaData(inteProObj,metadata = use_df)

inteProObj$anno_2 <- factor(inteProObj$anno_2,levels = c("Myeloid","T cell","B cell","NK","pDC","Epithelial","Endothelial","Fibroblast","Plasma","Cycling cell"))

###Fig.S3D------------------
col_use = c(brewer.pal(5, 'Paired'),brewer.pal(8, 'Set2')[c(7,4,3,2,1)])#
umapGG = CellDimPlot(inteProObj, group.by = "anno_2", reduction = "UMAP", theme_use = "theme_blank",label = F,palcolor = col_use)
ggsave(umapGG,filename = paste(tmp_pf,"UMAP_anno_2_260121_theme.png", sep = ""),dpi = 600, height = 8, width = 6,units = "in")
ggsave(umapGG,filename = paste(tmp_pf,"UMAP_anno_2_260121_theme.pdf", sep = ""),height = 8, width = 6,units = "in")

intgene <- c("PTPRC","EPCAM","CDH1","KRT14","KRT16","CAPS","SNTN",#epithelial
             "CD33","CD14","CD163","CD68","S100A8","C1QA","APOE","KIT",
             "CD3D","CD8A","CD4",
             "NCAM1","FCGR3A","NKG7",
             "CD19","CD79A","MS4A1","IGHA1","JCHAIN",
             "COL1A1","COL1A2",#fibroblast
             "PECAM1","VWF",#endothelial
             "SELP","ITGA2B","PPBP",
             "CLEC4C","LILRA4",#pDC
             "SCGB1A1","ACTA2",
             "MARCO", "FOLR2", "PPARG", "CD169", "CCL18")
dotGG = DotPlot(inteProObj, features = intgene, assay = "RNA",group.by = "anno_2") + 
  theme_classic()+
  theme(axis.text.x = element_text(color = "black",family = "sans"),
        axis.text.y = element_text(color = "black",family = "sans"))+
  scale_color_distiller(palette = "RdBu")+
  RotatedAxis() +coord_flip()+scale_alpha_continuous(name = "Per Exp.")+
  guides(color = guide_colorbar(title = "Avg Exp."))
ggsave(plot = dotGG, filename = paste(tmp_pf, "dot_anno_notsure.png", sep = ""), units = "in",dpi = 600, width = 8, heigh = 12)

###subset myeloid cluster#########
intePro_CD33 <- subset(inteProObj,subset = anno_2%in%c("Myeloid"))
saveRDS(intePro_CD33,file = paste(tmp_pf,"CD33.rds",sep = ""))
