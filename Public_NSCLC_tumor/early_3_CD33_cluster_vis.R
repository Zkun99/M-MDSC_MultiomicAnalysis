# Analyze public tumor CD33 cluster

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
suppressMessages(library(ggsci))

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

inteObjs <- readRDS(file = "/Public_NSCLC_early_Tumor/LCNAIVETUMOR_49_SINGLET_241115_cluster_CD33.rds")
inteObjs@meta.data <- inteObjs@meta.data[,c(1:14)]

if (clusterFlag) {
  clstFolder <- paste(expID, "cluster_CD33", sep = "_")
  dir.create(file.path(expDir, clstFolder), showWarnings = FALSE)
  clstDir <- paste(expDir, clstFolder, sep = "/")
  tmp_pf <- paste(clstDir, "/", clstFolder, "_", sep = "")
  
  inteProObj <- srPostPro(inteObjs, plotPf = tmp_pf, rdsSave = TRUE, ndim = 30, jsFlag = TRUE, qcFlag = TRUE, adtPlot = TRUE,
                          plotFeat = sigMarker, metaPlot = c("pid","disease"))
}

inteProObj <- readRDS(file = paste(tmp_pf, "Seurat_Objects_Clustered.RDS", sep = ""))

###recluster CD33 Object and DEG----------
resScreen <- function(srsc, plotPf, plotDir, batchName = NULL, res = 300, ndim = 18, rdsSave = FALSE, plotFeat = NULL, AbFeat = NULL,jsFlag = FALSE, resolutions) {
  qcgg <- VlnPlot(srsc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = batchName)
  ggsave(plot = qcgg, filename = paste(plotPf, "QCPlot1.png", sep = ""), dpi = res, width = 9, height = 6) 
  
  for (ires in resolutions) {
    resFolder <- paste("CD33_sub_cluster", ires, sep = "_")
    dir.create(file.path(plotDir, resFolder), showWarnings = FALSE)
    resDir <- paste(plotDir, resFolder, sep = "/")
    rppf <- paste(resDir, "/", resFolder, "_", sep = "")
    
    srsc <- FindClusters(srsc, resolution = ires)
    
    srsc <- BuildClusterTree(srsc)
    png(paste(rppf, "cluster_tree.png", sep = ""), res = 300, width = 12, height = 18, unit = 'in')
    PlotClusterTree(object = srsc)
    gar <- dev.off()
    
    umapGG <- DimPlot(srsc, reduction = "umap", label = TRUE)
    ggsave(plot = umapGG, filename = paste(rppf, "UMAP.png", sep = ""), dpi = res, width = 9, heigh = 6)
    
    
    if (!is.null(plotFeat)) {
      cat("Start to plot feature plots and violin plots\n")
      vlnGG <- VlnPlot(srsc, features = plotFeat, slot = "data", log = FALSE, ncol = 3, pt.size = 0.5)
      ggsave(plot = vlnGG, filename = paste(rppf, "violin_plot_for_cell_annotation.png", sep = ""), 
             dpi = res, width = 16, height = 16) 
      
      featGG <- FeaturePlot(srsc, features = plotFeat, ncol = 3, pt.size = 0.2, sort.cell = TRUE)
      ggsave(plot = featGG, filename = paste(rppf, "umap_feature_plot_for_cell_annotation.png", sep = ""), 
             dpi = res, width = 16, height = 16)
      
      dotGG <- DotPlot(srsc, features = plotFeat) + RotatedAxis()+theme_classic()
      ggsave(plot = dotGG, filename = paste(rppf, "DotPlot_for_cell_annotation.png", sep = ""),
             dpi = res, width = 9, height = 6)
      
      rdgGG <- RidgePlot(srsc, features = plotFeat, ncol = 3)
      ggsave(paste(rppf, "ridge_plot_for_cell_annotation.png", sep = ""), rdgGG, dpi = res, width = 16, height = 18)
      
    }
    
    srsc.markers <- FindAllMarkers(srsc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    topMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    
    write.csv(srsc.markers, file = paste(rppf, "ALL_POS_MARKERS.csv", sep = "")) #ST1x
    write.csv(topMarkers, file = paste(rppf, "top10_pos_markers.csv", sep = "")) 
    
    topMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
    hm <- DoHeatmap(object = srsc, features = topMarkers$gene) + NoLegend()
    ggsave(paste(rppf, "top5_marker_heatmap.png", sep = ""), hm, dpi = res, width = 24, height = 18)
    
    topDotMarkers <- srsc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
    geneMarkers <- unique(as.character(topDotMarkers$gene))
    
    dotMarkerGG <- DotPlot(srsc, features = geneMarkers,assay = "RNA") + 
      theme_classic()+
      scale_color_distiller(palette = "RdBu")+
      RotatedAxis() + coord_flip()
    ggsave(plot = dotMarkerGG, filename = paste(rppf, "Top5Marker_DotPlot.png", sep = ""), 
           dpi = res, width = 9, height = 16)
    
  }
  return(srsc)
}


subFolder <- "CD33_sub_cluster"
dir.create(file.path(clstDir, subFolder), showWarnings = FALSE)
subDir <- paste(clstDir, subFolder, sep = "/")
subPrefix <- paste(subDir, "/", expID, "_", subFolder, "_", sep = "")

pst <- Sys.time()
ress <- c(seq(0.1,0.9, 0.1))

geneOi <- c("CD14", "CD68", "HLA-DRA", "CD274", "SIGLEC15", "PDCD1LG2", "CD1C", "ITGAX", "FCER1A", 
            "TPSAB1", "KIT", "FUT4","CD33","CCL5","GNLY","NKG7")


inteProObj <- resScreen(inteProObj, plotPf = subPrefix, plotDir = subDir, rdsSave = TRUE, jsFlag = FALSE, ndim = 10,
                        batchName = c("pid"), plotFeat = geneOi,#AbFeat = AbOi ,
                        resolutions = ress)
cat("Analysis time cost:")
print(Sys.time()-pst)

umapGG = CellDimPlot(inteProObj, group.by = "integrated_snn_res.0.3", reduction = "UMAP", label = TRUE, label_insitu = TRUE)
ggsave(plot = umapGG, filename = paste(tmp_pf, "UMAP_res0.3.png", sep = ""), dpi = 300, width = 8, heigh = 6)

intgene <- c("CD33","CD14","FCGR3A","CD68","S100A8","S100A9","SPP1","C1QA","C1QC","APOE","VCAN","FCN1",#单核巨噬
             "KIT","TPSAB1", #MAST cell
             "MS4A1", "CD27", "CD19", #plasma
             "CCL5","NKG7","IL32","CD3E","GZMA","GZMB",#CCL5 M-MDSC
             "PLD4","CD1C","ITGAX","LILRA4","HLA-DRA")#DC
umapGG = FeatureDimPlot(srt = inteProObj, features = intgene,theme_use = "theme_blank",ncol = 7)
ggsave(plot = umapGG, filename = paste(tmp_pf, "UMAP_int_gene.png", sep = ""), dpi = 300, width = 35, heigh = 16)

inteProObj$cell_anno = "NA"
celltype = c("Mph_APOE","Mph_FBP1","Mono_FCN1","Mono_CD14","cDC","MDSC-like","Plasma","Mast","Mph_C1Q",
             "Mono_CD16","DC2","DC_LAMP3")

sub_length = length(unique(inteProObj$integrated_snn_res.0.3)) - 1
for (i in 0:sub_length){ 
  inteProObj$cell_anno[inteProObj$integrated_snn_res.0.3==i] = celltype[i+1]
}

#######save CD33 data with annotation----------
saveRDS(inteProObj,file = "LCNAIVETUMOR_49_SINGLET_241115_cluster_CD33_Annotation_251017.rds")

inteProObj <- readRDS(file = "LCNAIVETUMOR_49_SINGLET_241115_cluster_CD33_Annotation_251017.rds")
use_df <- inteProObj@meta.data
use_df <- use_df %>%
  mutate(anno_2 = case_when(
    cell_anno %in% c("cDC","DC_LAMP3","DC2") ~ "DC",
    cell_anno %in% c("Mast") ~ "Mast",
    cell_anno %in% c("Mono_CD14") ~ "cMono",
    cell_anno %in% c("Mono_CD16") ~ "ncMono",
    cell_anno %in% c("Mph_APOE","Mono_FCN1","Mph_FBP1","Mph_C1Q","MDSC-like") ~ "Macrophage",
    cell_anno %in% c("Plasma") ~ "Plasma cell",
    TRUE ~ "unknown"  # 相当于else
  ))
table(use_df$anno_2)
inteProObj <- AddMetaData(inteProObj,metadata = use_df)

###Fig.S3E----------------------
col_use = c(brewer.pal(9, 'Set3')[c(3:9)])
umapGG = CellDimPlot(inteProObj, group.by = "anno_2", reduction = "UMAP", theme_use = "theme_blank",label = F,palcolor = col_use)
ggsave(umapGG,filename = paste(tmp_pf,"UMAP_anno_2_260128_theme.png", sep = ""),dpi = 600, height = 8, width = 6,units = "in")
ggsave(umapGG,filename = paste(tmp_pf,"UMAP_anno_2_260128_theme.pdf", sep = ""),height = 8, width = 6,units = "in")

Macro_genes <- c("CD68","C1QA","APOE")#,"MARCO","FABP4","PPARG"
cMono_genes <- c("FCN1","S100A8","S100A9")
ncMono_genes <- c("FCGR3A","LST1","LILRB2")
DC_genes <- c("CD1C","CST3","HLA-DQA1")
Mast_genes <- c("KIT","TPSAB1", "TPSB2")
plasma_genes <- c("JCHAIN","IGHA1","IGKC")

features <- list("Macro" = Macro_genes,
                 "cMono" = cMono_genes,
                 "ncMono" = ncMono_genes,
                 "DC" = DC_genes,
                 "Mast" = Mast_genes,
                 "plasma" = plasma_genes)

###Fig.S3F----------------------
dotGG <- DotPlot(object = inteProObj, features = features, assay = "RNA", group.by = "anno_2",
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
  scale_y_discrete(limits = c("Plasma cell","Mast","DC","ncMono","cMono","Macrophage"))
print(dotGG)
ggsave(plot = dotGG, filename = paste0(tmp_pf,"dot_ALL_CELL.png"),width = 8.5, heigh = 3,units = "in",dpi = 600)
ggsave(plot = dotGG, filename = paste0(tmp_pf,"dot_ALL_CELL.pdf"),width = 8.5, heigh = 3,units = "in")
