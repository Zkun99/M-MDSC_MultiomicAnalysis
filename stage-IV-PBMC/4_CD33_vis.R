####recluster CD33 Myeloid cell subset#######

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
  srsc <- FindClusters(srsc, resolution = 0.1)#change
  
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
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggThemeAssist))


inputDir <- "6CloupeFile/results/"
expID <- "LCCHEMOICBPBMC_12_SINGLET"

dir.create(file.path(inputDir, expID), showWarnings = FALSE)
expDir <- paste(inputDir, expID, sep = "/")
ppf <- paste(expDir, "/", expID, "_", sep = "")

clstFolder <- paste(expID, "cluster_01", sep = "_")
dir.create(file.path(expDir, clstFolder), showWarnings = FALSE)
clstDir <- paste(expDir, clstFolder, sep = "/")
tmp_pf <- paste(clstDir, "/", clstFolder, "_", sep = "")


inteObjs <- readRDS("LCCHEMOICBPBMC_4_SINGLET_221202_cluster_230110_Seurat_Objects_Clustered_anno.RDS")

if (clusterFlag) {
  clstFolder <- paste(expID, "cluster_230110", sep = "_")
  dir.create(file.path(expDir, clstFolder), showWarnings = FALSE)
  clstDir <- paste(expDir, clstFolder, sep = "/")
  tmp_pf <- paste(clstDir, "/", clstFolder, "_", sep = "")
  
  inteProObj <- srPostPro(inteObjs, plotPf = tmp_pf, rdsSave = TRUE, ndim = 30, jsFlag = TRUE, qcFlag = TRUE, adtPlot = TRUE,
                          plotFeat = sigMarker, metaPlot = c("orig.ident", "pre_post", "pid"))
}

use_df <- inteProObj@meta.data
use_df <- use_df %>%
  mutate(celltype = case_when(
    seurat_clusters %in% c(0) ~ "Classical monocyte-C1",
    seurat_clusters %in% c(1) ~ "Classical monocyte-C2",
    seurat_clusters %in% c(2) ~ "M-MDSC",
    seurat_clusters %in% c(3) ~ "Non-classical monocyte",
    seurat_clusters %in% c(4) ~ "Classical monocyte-C3", 
    seurat_clusters %in% c(5) ~ "cDC",
    seurat_clusters %in% c(6) ~ "mDC",
    seurat_clusters %in% c(7) ~ "pDC",
    seurat_clusters %in% c(8) ~ "HSC",
    seurat_clusters %in% c(9) ~ "Plasma cell",
    TRUE ~ "unknown"  
  ))
table(use_df$anno_gzk)
inteProObj <- AddMetaData(inteProObj,metadata = use_df)
inteProObj$celltype <- factor(inteProObj$celltype, levels = c("Classical monocyte-C1","Classical monocyte-C2","Classical monocyte-C3",
                                                              "Non-classical monocyte","M-MDSC","cDC","mDC","pDC","HSC","Plasma cell"))

MDSC_Obj <- subset(inteProObj,subset = celltype%in%c("M-MDSC"))
saveRDS(MDSC_Obj,file = "LCCHEMOICBPBMC_12_SINGLET_cluster_230110_MDSC.rds")

###Fig.1B########
col_ten <- c("#deebf7","#9ecae1","#4292c6","#08519c","#DC143C","#fee6ce","#FCED82","#f16913","#B383B9","#a50f15")
CellGG = CellDimPlot(inteProObj, group.by = "celltype", reduction = "UMAP",palcolor = col_ten,theme_use = "theme_blank")
ggsave(CellGG,filename = paste(tmp_pf,"celltype.png", sep = ""),width = 9, heigh = 6,units = "in",dpi = 600)
ggsave(CellGG,filename = paste(tmp_pf,"celltype.pdf", sep = ""),width = 9, heigh = 6,units = "in")


###Fig.1C-------------
inteProObj$celltype2 <- as.character(inteProObj$celltype)
inteProObj$celltype2 <- factor(inteProObj$celltype2,levels = c(A_reversed <- c("Plasma cell", "HSC", "pDC", "mDC", "cDC","Non-classical monocyte", 
                                                                               "Classical monocyte-C3","Classical monocyte-C2", 
                                                                               "Classical monocyte-C1", "M-MDSC")))

c1_genes <- c("LGALS2", "ID1","CDKN1A")
c2_genes <- c("S100A12","S100A8","S100A9")
c3_genes <- c("ISG15","IFI44L" ,"IFI6")#"MX1"
nc_genes <- c("FCGR3A", "LST1","LILRB2")
MDSC_genes <- c("CD14","CD68","ICAM1")
cDC_genes <- c("CD1C","HLA-DQA1", "FCER1A")
mDC_genes <- c("BEST1","HOOK2","S100A6")
pDC_gene <- c("PLD4","LILRA4","GZMB")
HSC_gene <- c("CD34","GATA2","CDK6")
plasma_gene <- c("JCHAIN","IGHA1","IGHA2")

features <- list("MDSC" = MDSC_genes,
                 "c1" = c1_genes,
                 "c2" = c2_genes,
                 "c3" = c3_genes,
                 "nc" = nc_genes,
                 "cDC" = cDC_genes,
                 "mDC" = mDC_genes,
                 "pDC" = pDC_gene,
                 "HSC" = HSC_gene,
                 "Plasma" = plasma_gene)

dotGG <- DotPlot(object = inteProObj, features = features, assay = "RNA", group.by = "celltype2") +
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
  guides(color = guide_colorbar(title = c("Avg Exp.")))
print(dotGG)
ggsave(plot = dotGG, filename = paste0(tmp_pf,"dot_myeloid_subset_251013.png"),width = 11, heigh = 3,units = "in",dpi = 600)
ggsave(plot = dotGG, filename = paste0(tmp_pf,"dot_myeloid_subset_251013.pdf"),width = 11, heigh = 3,units = "in")



####Fig.1F####
inteProObj_sel <- subset(inteProObj,subset = celltype%in%c("Classical monocyte-C1","Classical monocyte-C2","Classical monocyte-C3",
                                                           "Non-classical monocyte","M-MDSC"))
table(inteProObj_sel$celltype)
inteProObj_sel$type_major <- ifelse(inteProObj_sel$celltype%in%c("Classical monocyte-C1","Classical monocyte-C2",
                                                                 "Classical monocyte-C3"),"Classical monocyte",
                                    as.character(inteProObj_sel$celltype))
table(inteProObj_sel$type_major)
inteProObj_sel$type_major <- factor(inteProObj_sel$type_major,levels = c("Non-classical monocyte","Classical monocyte","M-MDSC"))

int_ab <- c("CD14-Ab","CD16-Ab","CD54-Ab","CD328-Ab")
dotGG = DotPlot(inteProObj_sel, features = int_ab, assay = "integrated_adt",group.by = "type_major",
                scale.min = 0,scale.max = 100) + 
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.text.x = element_text(color = "black",family = "sans"),
        axis.text.y = element_text(color = "black",family = "sans"),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        axis.ticks = element_line(colour = "black"),
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.1, size =16))+
  scale_color_distiller(palette = "RdBu")+
  RotatedAxis()+
  theme(plot.margin = margin(t = 0.5, r = 0.1, b = 0.1,l = 0.1,unit = "in"))
dotGG
ggsave(plot = dotGG, filename = paste0(tmp_pf,"Mono_integrated_Ab_dotplot_type_major_3.png"),dpi = 600, width = 6, heigh = 3,units = "in")
ggsave(plot = dotGG, filename = paste0(tmp_pf,"Mono_integrated_Ab_dotplot_type_major_3.pdf"),dpi = 600, width = 6, heigh = 3)

####Fig.7A-------------------
inteProObj_Major <- inteProObj
table(inteProObj_Major$celltype)
inteProObj_Major <- subset(inteProObj_Major,subset = celltype%in%c("Classical monocyte-C1","Classical monocyte-C2","Classical monocyte-C3",
                                                                   "Non-classical monocyte","M-MDSC","cDC","mDC","pDC"))
meta_df <- as.data.frame(inteProObj_Major@meta.data)
meta_df <- meta_df %>%
  mutate(celltype = case_when(
    celltype %in% c("Classical monocyte-C1","Classical monocyte-C2","Classical monocyte-C3") ~ c("Classical monocyte"),
    TRUE              ~ celltype
  ))
meta_df$type_major <- meta_df$celltype
meta_df <- meta_df %>%
  mutate(type_major = case_when(
  type_major %in% c("cDC","mDC","pDC") ~ c("DC"),
    TRUE              ~ type_major
  ))
table(meta_df$type_major)


MDSC_anno <- read.csv(file = "E:/NSCLC_PBMC_scRNA_seq/Monocle_1018/LCCHEMOICBPBMC_12_SINGLET_230110_MDSC_use_df_anno.csv",row.names = 1)
MDSC_anno <- MDSC_anno[,c("cell_id","anno")]
rownames(MDSC_anno) <- MDSC_anno$cell_id
MDSC_anno$anno <- as.character(MDSC_anno$anno)

meta_df <- merge(meta_df, MDSC_anno, by = "cell_id", all.x = TRUE)
rownames(meta_df) <- meta_df$cell_id
meta_df$anno <- ifelse(is.na(meta_df$anno), meta_df$type_major, meta_df$anno)
table(meta_df$anno)

meta_df$anno2 <- ifelse(meta_df$anno%in%c("DC"),meta_df$celltype,meta_df$anno)

inteProObj_Major <- AddMetaData(inteProObj_Major,metadata = meta_df)
# saveRDS(inteProObj_Major,file = "LCCHEMOICBPBMC_12_SINGLET_cluster_230110_Seurat_Objects_Clustered_anno_Mono+MDSC+DC_250905.RDS")
# inteProObj_Major <- readRDS(file = "LCCHEMOICBPBMC_12_SINGLET_cluster_230110_Seurat_Objects_Clustered_anno_Mono+MDSC+DC_250905.RDS")

inteProObj_Major$anno <- factor(inteProObj_Major$anno,levels = c("DC",
                                                                 "Non-classical monocyte",
                                                                 "Classical monocyte","S100A8+ MDSC",
                                                                 "GZMB+ MDSC","JCHAIN+ MDSC","CCL5+ MDSC"))
Co_inhibitory <- c("B7-H4-Ab","CD274-Ab","CD273-Ab","CD276-Ab","CD155-Ab","CD112-Ab")
Co_stimulatory <- c("CD80-Ab","CD86-Ab","CD154-Ab","CD252-Ab","CD275-Ab","CD70-Ab","CD137L-Ab")
Siglec_family <- c("CD169-Ab","CD328-Ab","Siglec-10-Ab")

features <- list("Co-inhibitory ligand" = Co_inhibitory,
                 "Co-stimulatory ligand" = Co_stimulatory,
                 "Siglec family" = Siglec_family)

dotGG = DotPlot(inteProObj_Major, features = features, assay = "integrated_adt",group.by = "anno",
                scale.min = 0,scale.max = 100) + 
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
  scale_color_distiller(palette = "RdBu")+
  RotatedAxis()+
  theme(plot.margin = margin(t = 0.5, r = 0.1, b = 0.1,l = 0.1,unit = "in"))
dotGG
ggsave(plot = dotGG, filename = paste0(tmp_pf,"Major_integrated_Ab_coLR_ALL_251124.png"),dpi = 600, width = 9, heigh = 3.5,units = "in")
ggsave(plot = dotGG, filename = paste0(tmp_pf,"Major_integrated_Ab_coLR_ALL_251124.pdf"),width = 9, heigh = 3.5)

