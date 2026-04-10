###recluster MDSC subset 

rm(list = ls())
suppressMessages(library(viridis))
suppressMessages(library(tradeSeq))
suppressMessages(library(slingshot))
suppressMessages(library(ggpubr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(DoubletFinder))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratData))
suppressMessages(library(SeuratDisk))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(scater))
suppressMessages(library(scran))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(SCP))
suppressMessages(library(ComplexHeatmap))

outpfil = "./MDSC_0526/"
filna = "LCCHEMOICBPBMC_12_SINGLET_cluster_230110_MDSC_"
tmp_pf <- paste0(outpfil,filna)
set.seed(1)

###### Load raw data ####
Lung_MDSC <- readRDS('LCCHEMOICBPBMC_12_SINGLET_cluster_230110_MDSC.rds')
Lung_MDSC@meta.data$response <- ifelse(Lung_MDSC@meta.data$pid %in% c("2","4","6"), "NR", "R")
Lung_MDSC@meta.data$condition <- str_c(Lung_MDSC@meta.data$pre_post , "+", Lung_MDSC@meta.data$response)

Lung_MDSC <- JackStraw(Lung_MDSC, num.replicate = 100)
Lung_MDSC <- ScoreJackStraw(Lung_MDSC, dims = 1:20)
jsGG <- JackStrawPlot(Lung_MDSC, dims = 1:20)+theme_bw()
ggsave(plot = jsGG, filename = paste(outpfil, filna,"JackStraw.png", sep = ""), dpi = 600, width = 9, heigh = 6)
elGG <- ElbowPlot(Lung_MDSC, ndims = 50)+theme_bw()
ggsave(plot = elGG, filename = paste(outpfil, filna, "Elbow.png", sep = ""), dpi = 600, width = 9, heigh = 6)

Lung_MDSC <- RunUMAP(Lung_MDSC, dims = 1:10)
Lung_MDSC <- FindNeighbors(Lung_MDSC, dims = 1:10)
Lung_MDSC <- FindClusters(Lung_MDSC, resolution = 0.3)#change
saveRDS(object = Lung_MDSC,file = paste(outpfil,filna,"res0.3.rds",sep = ""))

####load recluter MDSC dataset----------------
rm(list = ls())
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(Seurat))
suppressMessages(library(SCP))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ComplexHeatmap))

inputDir <- "6CloupeFile/results/"
expID <- "LCCHEMOICBPBMC_12_SINGLET_230110"

dir.create(file.path(inputDir, expID), showWarnings = FALSE)
expDir <- paste(inputDir, expID, sep = "/")
ppf <- paste(expDir, "/", expID, "_", sep = "")

clstFolder <- paste(expID, "MDSC", sep = "_")
dir.create(file.path(expDir, clstFolder), showWarnings = FALSE)
clstDir <- paste(expDir, clstFolder, sep = "/")
tmp_pf <- paste(clstDir, "/", clstFolder, "_", sep = "")

####1.load MDSC dataset#######
late_PBMC_MDSC <- readRDS("E:/NSCLC_PBMC_scRNA_seq/MDSC_0526/LCCHEMOICBPBMC_12_SINGLET_cluster_230110_MDSC_res0.3.rds")
print(late_PBMC_MDSC)

anno_MDSC <- c("CCL5+ MDSC","B-like MDSC","S100A8+ MDSC","pDC-like MDSC")
late_PBMC_MDSC$anno <- late_PBMC_MDSC$integrated_snn_res.0.3
levels(late_PBMC_MDSC$anno) <- anno_MDSC

late_PBMC_MDSC$anno_1 <- as.character(late_PBMC_MDSC$anno)
late_PBMC_MDSC$anno_1 <- factor(late_PBMC_MDSC$anno_1,levels = c("CCL5+ MDSC",
                                                                 "B-like MDSC",
                                                                 "pDC-like MDSC",
                                                                 "S100A8+ MDSC"))

table(late_PBMC_MDSC$anno,late_PBMC_MDSC$pre_post)

#####save data for RNA velocity##########
meta_RNA_velo <- as.data.frame(late_PBMC_MDSC@meta.data)
meta_RNA_velo <- meta_RNA_velo[,c("pre_post","orig.ident","response","integrated_snn_res.0.3","anno")]
meta_RNA_velo$cell_id <- rownames(meta_RNA_velo)

Embd_df <- as.data.frame(Embeddings(late_PBMC_MDSC, reduction = "umap"))
Embd_df$cell_id <- rownames(Embd_df)
meta_RNA_velo <- merge(meta_RNA_velo,Embd_df,by = "cell_id")
#meta_RNA_velo$cell_id <- substr(meta_RNA_velo$cell_id,1,16)
write.csv(meta_RNA_velo,file = paste(tmp_pf,"use_df_anno.csv",sep = ""),row.names = T)

####2..1 vis#######
###Fig.2A-----------
col_four <- c("#E41A1C","#377EB8","#66C2A5","#E78AC3")
CellGG = CellDimPlot(late_PBMC_MDSC, group.by = "anno", reduction = "UMAP",palcolor = col_four, theme_use = "theme_blank")
ggsave(CellGG,filename = paste(tmp_pf,"anno_241218.png", sep = ""),width = 5, heigh = 4,units = "in",dpi = 600)
ggsave(CellGG,filename = paste(tmp_pf,"anno_241218.pdf", sep = ""),width = 5, heigh = 4,units = "in")



####2.2 DEG###########
DefaultAssay(late_PBMC_MDSC) <- "integrated"
Idents(late_PBMC_MDSC) <- late_PBMC_MDSC$anno_1
late_PBMC_MDSC.markers = FindAllMarkers(late_PBMC_MDSC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,group.by = "anno_1")
topMarkers = late_PBMC_MDSC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(late_PBMC_MDSC.markers, file = paste(tmp_pf, "ALL_pos_markers.csv", sep = ""))
write.csv(topMarkers, file = paste(tmp_pf, "ALL_top_pos_markers.csv", sep = ""))
print(head(topMarkers))

Idents(late_PBMC_MDSC) <- late_PBMC_MDSC$anno_1
late_PBMC_MDSC.markers = FindAllMarkers(late_PBMC_MDSC,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,assay = "integrated_adt")
write.csv(late_PBMC_MDSC.markers, file = paste(tmp_pf, "ADT_ALL_POS_MARKERS_integrated.csv", sep = ""))


###Fig.2B---------------
use_gene <- c("CCL5","NKG7","GZMA","IL32","CD3E",#CCL5 cluster 
              "JCHAIN","IGHA1","CD79A","CD19","MS4A1",#B-like cluster
              "PLD4","IL3RA","LILRA4","CLEC4C",#pDC-like cluster
              "S100A8","CTSD","VCAN","CSF3R","FCN1")#S100A8 cluster
ht1 <- GroupHeatmap(late_PBMC_MDSC,features = use_gene,assay = "RNA",heatmap_palette = "viridis",group.by = "anno_1",
                    show_row_names = T)
ht1$plot
ggsave(filename = paste(tmp_pf, "heatmap_select_gene_251013.png", sep = ""), dpi = 600, width = 5, heigh = 5,units = "in")
ggsave(filename = paste(tmp_pf, "heatmap_select_gene_251013.pdf", sep = ""), width = 5, heigh = 5,units = "in")

#2.3.Maturation plot-----
#Fig.2I------------------
library(ggbreak)
intAb <- c("CD54-Ab")
late_PBMC_MDSC$anno_2 <- factor(late_PBMC_MDSC$anno,levels = c("CCL5+ MDSC","pDC-like MDSC","B-like MDSC","S100A8+ MDSC"))
featGG_CD54 = FeatureStatPlot(srt = late_PBMC_MDSC, stat.by = c("CD54-Ab"),group.by = "anno_2",plot_type = "box",
                              assay = "integrated_adt",bg_palcolor = "white",palcolor = c("#E41A1C","#E78AC3","#377EB8","#66C2A5"),
                              # comparisons = list(c("CCL5+ MDSC","pDC-like MDSC"),c("CCL5+ MDSC","B-like MDSC"),c("CCL5+ MDSC","S100A8+ MDSC")),
                              # sig_label = "p.format",pairwise_method = "wilcox.test"
                              )+
  theme_scp()+
  coord_cartesian(ylim = c(0, 2))+
  scale_y_continuous(breaks = seq(0, 2, 0.5))
featGG_CD54
ggsave(plot = featGG_CD54, filename = paste0(tmp_pf,"dot_CD54-Ab_251125.png"),width = 5, heigh = 4,units = "in",dpi = 600)
ggsave(plot = featGG_CD54, filename = paste0(tmp_pf,"dot_CD54-Ab_251125.pdf"),width = 5, heigh = 4,units = "in")


#2.4 T/B/DC marker Radar----
suppressMessages(library(scales))
suppressMessages(library(dplyr))
suppressMessages(library(pheatmap))
suppressMessages(library(gridExtra))

T_marker <- c("CD8A","CD3D","CD3E","CD3G","GZMA","TRBC1","TRBC2","TRAC","CD4","GNLY",
              "PFA","CD7","NKG7","PRF1","CST7")
B_marker <- c("JCHAIN","CD79A","CD19","MS4A1")
pDC_marker <- c("PLD4","IL3RA","LILRA4","CLEC4C","SLC15A4","CCDC50","LY9", "SELL", "GAS6")
Mono_marker <- c("CD14","S100A8","S100A9","FCN1","LYZ")

# cal_Macro_sig_score function--------------
# Cheng S, Li Z, Gao R et. al., A pan-cancer single-cell transcriptional atlas of tumor infiltrating myeloid cells, Cell
cal_Macro_sig_score <- function(seurat_obj, clusters, genes) {
  genes_filt <- genes[genes %in% rownames(seurat_obj@assays$RNA@data)]
  expression_matrix <- seurat_obj@assays$RNA@data[genes_filt, seurat_obj@meta.data$anno %in% clusters]

  scores <- as.data.frame(colMeans(as.matrix(expression_matrix)))
  colnames(scores) <- "value"
  scores$cluster <- seurat_obj@meta.data[colnames(expression_matrix), "anno"]

  average_score <- scores %>% group_by(cluster) %>% summarize(sig_score = mean(value))
  average_score$sig_score_scaled <- rescale(average_score$sig_score, to = c(0, 5))
  return(average_score)
}

T_score <- cal_Macro_sig_score(late_PBMC_MDSC, unique(late_PBMC_MDSC@meta.data$anno), T_marker)
B_score <- cal_Macro_sig_score(late_PBMC_MDSC, unique(late_PBMC_MDSC@meta.data$anno), B_marker)
pDC_score <- cal_Macro_sig_score(late_PBMC_MDSC, unique(late_PBMC_MDSC@meta.data$anno), pDC_marker)
Mono_score <- cal_Macro_sig_score(late_PBMC_MDSC, unique(late_PBMC_MDSC@meta.data$anno), Mono_marker)

res_raw_1 <- merge(T_score[, c("cluster", "sig_score")],B_score[, c("cluster", "sig_score")],by = "cluster", suffixes = c("_T", "_B"))
res_raw_2 <- merge(pDC_score[, c("cluster", "sig_score")],Mono_score[, c("cluster", "sig_score")],by = "cluster", suffixes = c("_pDC", "_Mono"))
res_raw <- merge(res_raw_1,res_raw_2,by = "cluster")

res_scaled_1 <- merge(T_score[, c("cluster", "sig_score_scaled")],B_score[, c("cluster", "sig_score_scaled")],by = "cluster", suffixes = c("_T", "_B"))
res_scaled_2 <- merge(pDC_score[, c("cluster", "sig_score_scaled")],Mono_score[, c("cluster", "sig_score_scaled")],by = "cluster", suffixes = c("_pDC", "_Mono"))
res_scaled <- merge(res_scaled_1,res_scaled_2,by = "cluster")

row.names(res_raw) <- res_raw$cluster
row.names(res_scaled) <- res_scaled$cluster
res_raw$cluster <- NULL
res_scaled$cluster <- NULL
write.csv(res_raw,file = paste0(tmp_pf,"_score_raw_250305.csv"),row.names = T)
write.csv(res_scaled,file = paste0(tmp_pf,"_score_scaled_250305.csv"),row.names = T)

##2.4.1.Radar--------
suppressMessages(library(ggradar))
library(tidyr)
library(dplyr)

res_scaled <- read.csv(file = paste0(tmp_pf,"_score_scaled_250305.csv"),row.names = 1)
Radar_df <- res_scaled%>%rownames_to_column("group")
colors <- c("#377EB8","#E41A1C","#E78AC3","#66C2A5")
titles <- c("B-like MDSC", "CCL5+ MDSC", "pDC-like MDSC","S100A8+ MDSC")

df <- Radar_df

##Fig.2C-------------
for (j in 1:length(titles)) {
  current_title <- titles[j]
  current_color <- colors[j]
  
  row_data <- df[df$group %in% current_title, ]
  
  Radarplot <- ggradar(row_data, 
                       values.radar = c("0", "2.5", "5"),
                       grid.min = 0, 
                       grid.mid = 2.5, 
                       grid.max = 5,
                       group.line.width = 2, 
                       group.point.size = 4,
                       background.circle.colour = "white",
                       gridline.mid.colour = "lightgrey",
                       group.colours = current_color
  ) + ggtitle(label = current_title)
  
  ggsave(plot = Radarplot,filename = paste0(tmp_pf, current_title, "_Rader.pdf"), height = 5,width = 8,units = "in")
  ggsave(plot = Radarplot,filename = paste0(tmp_pf, current_title, "_Rader.png"), dpi = 600,height = 5,width = 8,units = "in")
}

##SCENIC---------------
library(pheatmap)
library(RColorBrewer)
library(viridis)

auc_df <- read.csv("MDSC_0526/pyscenic/output/auc_mtx.csv", row.names = 1)
colnames(auc_df) <- gsub("\\.\\.\\.$", "", colnames(auc_df))
head(auc_df)
colnames(auc_df)
MDSC_AUC <- AddMetaData(late_PBMC_MDSC,metadata = auc_df)

matched_auc <- auc_df[colnames(late_PBMC_MDSC), ]
auc_matrix <- as.matrix(t(matched_auc))  

MDSC_AUC[["AUC"]] <- CreateAssayObject(data = auc_matrix)
DefaultAssay(MDSC_AUC) <- "AUC"
MDSC_AUC@misc$AUC_matrix <- matched_auc
MDSC_AUC@misc$AUC_assay <- auc_matrix
cell_types <- unique(MDSC_AUC@meta.data$anno_1)

top_regulons <- c("TBX21","STAT5A","SMAD4","FOXO1",
                  "NFATC2","NFIA","RFX5",
                  "LEF1","RFX1","TCF7L2","IKZF1",
                  "SPI1","CEBPA","CEBPD","CEBPG","IRF2","IRF3","IRF8")

celltype_mean_auc <- lapply(cell_types, function(ct) {
  cells <- rownames(MDSC_AUC@meta.data)[MDSC_AUC@meta.data$anno_1 == ct]
  colMeans(matched_auc[cells, top_regulons, drop = FALSE], na.rm = TRUE)
})
names(celltype_mean_auc) <- cell_types

# 
mean_auc_matrix <- do.call(cbind, celltype_mean_auc)
mean_auc_matrix <- t(mean_auc_matrix)  

# Z-score
zscore_normalize_columns <- function(mat) {
  col_means <- apply(mat, 2, mean, na.rm = TRUE)
  col_sds <- apply(mat, 2, sd, na.rm = TRUE)
  col_sds[col_sds == 0] <- 0.001  
  
  # Z-score = (x - mean) / sd
  scaled_mat <- t((t(mat) - col_means) / col_sds)
  
  scaled_mat[scaled_mat > 3] <- 3
  scaled_mat[scaled_mat < -3] <- -3
  
  return(scaled_mat)
}


scaled_matrix <- zscore_normalize_columns(mean_auc_matrix)
desired_celltype_order <- c("CCL5+ MDSC", "B-like MDSC", "pDC-like MDSC", "S100A8+ MDSC")  
desired_celltype_order <- desired_celltype_order[desired_celltype_order %in% cell_types]

scaled_matrix <- scaled_matrix[desired_celltype_order, , drop = FALSE]

cell_annotation <- data.frame(
  CellType = rownames(scaled_matrix),
  row.names = rownames(scaled_matrix)
)

celltype_colors <- setNames(c("#E78AC3","#66C2A5","#E41A1C","#377EB8"), cell_types)
ann_colors <- list(CellType = celltype_colors)

heat2 = pheatmap(scaled_matrix,
                 cluster_rows = FALSE,  
                 cluster_cols = FALSE,   
                 show_colnames = TRUE,
                 show_rownames = TRUE,
                 annotation_row = cell_annotation,
                 annotation_colors = ann_colors,
                 main = "Regulon Activity (Mean by Cell Type, Z-score normalized)",
                 fontsize_row = 10,
                 fontsize_col = 9,
                 angle_col = 45,
                 border_color = NA,
                 scale = "none",  
                 display_numbers = FALSE,
                 gaps_row = 1:(nrow(scaled_matrix)-1))  
heat2
ggsave(plot = heat2,filename = paste(tmp_pf,"select_fuction_regulon_cluster_row_FALSE_3.png",sep = ""),width = 8, heigh = 4,units = "in",dpi = 600)
ggsave(plot = heat2,filename = paste(tmp_pf,"select_fuction_regulon_cluster_row_FALSE_3.pdf",sep = ""),width = 8, heigh = 4,units = "in")

