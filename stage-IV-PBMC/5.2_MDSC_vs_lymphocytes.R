### Fig.2D+E MDSC compared with lymphocytes and monocytes-------

rm(list = ls())
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
suppressMessages(library(data.table))
suppressMessages(library(ggbreak))

inputDir <- "E:/NSCLC_PBMC_scRNA_seq/PBMC_ALL_CELL"
expID <- "LCCHEMOICBPBMC_12_SINGLET_221202"
tmp_pf <- paste(inputDir,expID,sep = "/")
set.seed(1)


CD33_All <- readRDS("LCCHEMOICBPBMC_4_SINGLET_221202_cluster_230110_Seurat_Objects_Clustered_CD33_subtype_anno.RDS")
cellty_meta <- read.csv("E:/NSCLC_PBMC_scRNA_seq/PBMC_ALL_CELL/celltypist_immune_low/predicted_labels.csv",row.names = 1)
cellty_meta$major_cell_type <- ifelse(cellty_meta$majority_voting%in%c("Tcm/Naive cytotoxic T cells",
                                                                       "Tcm/Naive helper T cells",
                                                                       "Tem/Temra cytotoxic T cells",
                                                                       "Tem/Trm cytotoxic T cells"),
                                      "T cells",ifelse(cellty_meta$majority_voting%in%c("CD16+ NK cells","NK cells"),
                                                       "NK cells",ifelse(cellty_meta$majority_voting%in%c("DC1","DC2"),
                                                              "DC",cellty_meta$majority_voting)))

table(cellty_meta$major_cell_type)
CD33_All <- AddMetaData(CD33_All,metadata = cellty_meta)

Tcell_Obj <- subset(CD33_All,subset = major_cell_type%in%c("T cells"))
Bcell_Obj <- subset(CD33_All,subset = anno_1%in%c('Naive B cells',"Memory B cells"))
Tcell_Obj@meta.data$anno <- c("T cell")
Bcell_Obj@meta.data$anno <- c("B cell")

##load M-MDSC data
late_PBMC_MDSC <- readRDS("MDSC_0526/LCCHEMOICBPBMC_12_SINGLET_cluster_230110_MDSC_res0.3.rds")
anno_MDSC <- c("CCL5+ MDSC","B-like MDSC","S100A8+ MDSC","pDC-like MDSC")
late_PBMC_MDSC$anno <- late_PBMC_MDSC$integrated_snn_res.0.3
levels(late_PBMC_MDSC$anno) <- anno_MDSC


#load monocyte and pDC data
inteProObj_sel <- readRDS(file = "LCCHEMOICBPBMC_12_SINGLET_cluster_230110_Seurat_Objects_Clustered_anno_Myeloid.RDS")
Mono_Obj <- subset(inteProObj_sel,subset = celltype%in%c("Classical monocyte-C3","Classical monocyte-C2","Classical monocyte-C1"))
pDC_Obj <- readRDS(file = "LCCHEMOICBPBMC_12_SINGLET_cluster_230110_pDC.RDS")


Mono_Obj@meta.data$anno <- c("cMono")
pDC_Obj@meta.data$anno <- c("pDC")

CCL5_Obj <- subset(late_PBMC_MDSC,subset = anno%in%"CCL5+ MDSC")
B_like_Obj <- subset(late_PBMC_MDSC,subset = anno%in%"B-like MDSC")
pDC_like_Obj <- subset(late_PBMC_MDSC,subset = anno%in%"pDC-like MDSC")
S100A8_Obj <- subset(late_PBMC_MDSC,subset = anno%in%"S100A8+ MDSC")

###CCL5+ MDSC compared with cMono and cytotoxic T cells
Sel_T_cell <- subset(Tcell_Obj,subset = majority_voting%in%c("Tem/Temra cytotoxic T cells","Tem/Trm cytotoxic T cells"))
table(Sel_T_cell$majority_voting)
Sel_T_cell$anno <- c("Tem cell")

Sel_T_vs_MDSC <- merge(CCL5_Obj,Sel_T_cell)
Sel_T_vs_MDSC <- merge(Sel_T_vs_MDSC,Mono_Obj)
table(Sel_T_vs_MDSC$anno)
Sel_T_vs_MDSC$anno <- factor(Sel_T_vs_MDSC$anno,levels = c("CCL5+ MDSC","cMono","Tem cell"))
intgene <- c("CCL5","NKG7","GZMA","CD3E","IL32")
compare_list <- list(c("CCL5+ MDSC","cMono"),c("cMono","Tem cell"),c("CCL5+ MDSC","Tem cell"))

#Fig.2D T cell gene------
featGG_T = FeatureStatPlot(srt = Sel_T_vs_MDSC, group.by = "anno",stat.by = intgene,plot_type = c("violin"),add_point = F,pt.color = "grey90",
                         pt.size = 0.01,pt.alpha = 1,jitter.width = 0.4,
                         bg_palcolor = "white",assay = "RNA",palcolor = c("#E41A1C","#76b0d6","#984EA3"),add_box = T,
                         # box_color = "white",box_width = 0.1,
                         # box_ptsize = 2,
                         ncol = 5,
                         comparisons = compare_list,sig_label = c("p.format"))&NoLegend()
ggsave(featGG_T,filename = paste(tmp_pf,"fig2_Tem_pvalue_gene_251121.png", sep = ""),dpi = 600, width = 20, heigh = 4,units = "in")
ggsave(featGG_T,filename = paste(tmp_pf,"fig2_Tem_pvalue_gene_251121.pdf", sep = ""),width = 20, heigh = 4,units = "in")

#Fig.2E T cell Antibody------
featGG_CD8 = FeatureStatPlot(srt = Sel_T_vs_MDSC, stat.by = c("CD8-Ab"),group.by = "anno",plot_type = "box",
                         assay = "integrated_adt",bg_palcolor = "white",palcolor = c("#E41A1C","#76b0d6","#984EA3"),
                         # add_point = T,pt.color = "grey",
                         # pt.size = 0.01,
                         # pt.alpha = 0.5,
                         # comparisons = compare_list,sig_label = "p.format",pairwise_method = "wilcox.test"
                         )+
  theme_classic()+
  theme(axis.text.x = element_text(color = "black",family = "sans",angle = 30,vjust = 0.85,hjust = 0.75),
        axis.text.y = element_text(color = "black",family = "sans"),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",linetype=1,color="black",size=0.25))+
  coord_cartesian(ylim = c(-0.5, 2.25))+
  scale_y_continuous(breaks = seq(-0.5, 2.25,0.5))+
  scale_y_break(c(1,2),space = 0.1,scales = 0.1)

ggsave(featGG_CD8,filename = paste(tmp_pf,"fig2_Tem_CD8-Ab_box_251125.png", sep = ""), dpi = 600,width = 4, heigh = 4,units = "in")
ggsave(featGG_CD8,filename = paste(tmp_pf,"fig2_Tem_CD8-Ab_box_251125.pdf", sep = ""), width = 4, heigh = 4,units = "in")

###B-like MDSC vs cMono and B cell
B_vs_MDSC <- merge(B_like_Obj,Bcell_Obj)
B_vs_MDSC <- merge(B_vs_MDSC,Mono_Obj)
table(B_vs_MDSC$anno)
color_3 <- c("#377EB8","#76b0d6","#74AC38")
intgene <- c("CD19","MS4A1","CD79A")

B_vs_MDSC$anno <- factor(B_vs_MDSC$anno,levels = c("B-like MDSC","cMono","B cell"))
compare_list <- list(c("B-like MDSC","cMono"),c("cMono","B cell"),c("B-like MDSC","B cell"))

#Fig.2D B cell gene------
featGG_B = FeatureStatPlot(srt = B_vs_MDSC, group.by = "anno",stat.by = c("CD19","MS4A1","CD79A"), 
                         plot_type = c("violin"),add_point = F,pt.color = "grey90",
                         pt.size = 0.01,pt.alpha = 1,jitter.width = 0.4,bg_palcolor = "white",assay = "RNA",
                         palcolor = c("#377EB8","#76b0d6","#74AC38"),add_box = T,
                         # box_color = "white",box_width = 0.1,
                         # box_ptsize = 2,
                         ncol = 3,
                         comparisons = compare_list,sig_label = c("p.format"))&NoLegend()

ggsave(featGG_B,filename = paste(tmp_pf,"fig2_Bcell_pvalue_251121.png", sep = ""), dpi = 600, width = 15, heigh = 10,units = "in")
ggsave(featGG_B,filename = paste(tmp_pf,"fig2_Bcell_pvalue_251121.pdf", sep = ""), width = 15, heigh = 10,units = "in")

#Fig.2E B cell Antibody------
intAb <- c("CD19-Ab")
featGG_CD19 = FeatureStatPlot(srt = B_vs_MDSC, stat.by = c("CD19-Ab"),group.by = "anno",plot_type = "box",
                              assay = "integrated_adt",bg_palcolor = "white",palcolor = c("#377EB8","#76b0d6","#74AC38"),
                              # comparisons = compare_list,sig_label = "p.format",pairwise_method = "wilcox.test"
                              )+
  theme_classic()+
  theme(axis.text.x = element_text(color = "black",family = "sans",angle = 30,vjust = 0.85,hjust = 0.75),
        axis.text.y = element_text(color = "black",family = "sans"),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",linetype=1,color="black",size=0.25))+
  scale_y_continuous(limits = c(-0.5, 3),breaks = seq(-0.5, 3,0.5))+
  scale_y_break(c(2,3),space = 0.1,scales = 0.1)

ggsave(featGG_CD19,filename = paste(tmp_pf,"fig2_Bcell_Ab_251125.png", sep = ""), dpi = 600, width = 4, heigh = 4,units = "in")
ggsave(featGG_CD19,filename = paste(tmp_pf,"fig2_Bcell_Ab_251125.pdf", sep = ""), width = 4, heigh = 4,units = "in")

###:B-like M-MDSC vs plasma cell
plasma_Obj <- readRDS(file = "LCCHEMOICBPBMC_12_SINGLET_cluster_230110_plasma.RDS")
plasma_Obj@meta.data$anno <- c("plasma")

Plasma_vs_MDSC <- merge(B_like_Obj,plasma_Obj)
Plasma_vs_MDSC <- merge(Plasma_vs_MDSC,Mono_Obj)
table(Plasma_vs_MDSC$anno)

Plasma_vs_MDSC$anno <- factor(Plasma_vs_MDSC$anno,levels = c("B-like MDSC","cMono","plasma"))
color_3 <- c("#377EB8","#76b0d6","#a50f15")

#Fig.2D plasma gene------
compare_list <- list(c("B-like MDSC","cMono"),c("cMono","plasma"),c("B-like MDSC","plasma"))
featGG_PLAS = FeatureStatPlot(srt = Plasma_vs_MDSC, group.by = "anno",stat.by = c("JCHAIN","IGHA1"), 
                         plot_type = c("violin"),add_point = F,pt.color = "grey90",
                         pt.size = 0.01,pt.alpha = 1,jitter.width = 0.4,bg_palcolor = "white",assay = "RNA",
                         palcolor = c("#377EB8","#76b0d6","#a50f15"),add_box = T,
                         # box_color = "white",box_width = 0.1,
                         # box_ptsize = 2,
                         ncol = 2,
                         comparisons = compare_list,sig_label = c("p.format"))&NoLegend()

ggsave(featGG_PLAS,filename = paste(tmp_pf,"fig2_Bcellvs_plasma_pvalue_251121.png", sep = ""), dpi = 600, width = 6, heigh = 4,units = "in")
ggsave(featGG_PLAS,filename = paste(tmp_pf,"fig2_Bcellvs_plasma_pvalue_251121.pdf", sep = ""), width = 6, heigh = 4,units = "in")


###PLD4+ M-MDSC vs pDC
pDC_vs_MDSC <- merge(pDC_like_Obj,pDC_Obj)
pDC_vs_MDSC <- merge(pDC_vs_MDSC,Mono_Obj)
color_3 <- c("#E78AC3","#76b0d6","#f16913")
pDC_vs_MDSC$anno <- factor(pDC_vs_MDSC$anno,levels = c("pDC-like MDSC","cMono","pDC"))

#Fig.2D pDC gene------
compare_list <- list(c("pDC-like MDSC","cMono"),c("cMono","pDC"),c("pDC-like MDSC","pDC"))
featGG_pDC = FeatureStatPlot(srt = pDC_vs_MDSC, group.by = "anno",stat.by = c("PLD4","IL3RA","LILRA4","CLEC4C"), 
                         plot_type = c("violin"),add_point = F,pt.color = "grey90",
                         pt.size = 0.01,pt.alpha = 1,jitter.width = 0.4,bg_palcolor = "white",assay = "RNA",
                         palcolor = c("#E78AC3","#76b0d6","#f16913"),
                         add_box = T,ncol = 5,
                         comparisons = compare_list,sig_label = c("p.format"))&NoLegend()
ggsave(featGG_pDC,filename = paste(tmp_pf,"fig2_pDC_pvalue_251121.png", sep = ""), dpi = 600, width = 15, heigh = 10,units = "in")
ggsave(featGG_pDC,filename = paste(tmp_pf,"fig2_pDC_pvalue_251121.pdf", sep = ""), width = 15, heigh = 10,units = "in")


#Fig.2E pDC Anitibody------
intAb <- c("CD123-Ab")
featGG_CD123 = FeatureStatPlot(srt = pDC_vs_MDSC, stat.by = c("CD123-Ab"),group.by = "anno",plot_type = "box",
                         assay = "integrated_adt",bg_palcolor = "white",palcolor = c("#E78AC3","#76b0d6","#f16913"),
                         # comparisons = compare_list,sig_label = "p.format",pairwise_method = "wilcox.test"
                         )+
  theme_classic()+
  theme(axis.text.x = element_text(color = "black",family = "sans",angle = 30,vjust = 0.85,hjust = 0.75),
        axis.text.y = element_text(color = "black",family = "sans"),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",linetype=1,color="black",size=0.25))+
  coord_cartesian(ylim = c(-0.5, 3.5))+
  scale_y_continuous(breaks = seq(-0.5,3.5,0.5))+
  scale_y_break(c(1.0,1.6),space = 0.1,scales = 0.5)
featGG_CD123
ggsave(featGG_CD123,filename = paste(tmp_pf,"fig2_pDC_Ab_251126.png", sep = ""), dpi = 600, width = 4, heigh = 4,units = "in")
ggsave(featGG_CD123,filename = paste(tmp_pf,"fig2_pDC_Ab_251126.pdf", sep = ""), width = 4, heigh = 4,units = "in")

