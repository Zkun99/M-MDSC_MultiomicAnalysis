### public late-stage MDSC and public early-stage MDSC ######
### function ref:Lodi et.al.,CRM 2025 

rm(list = ls())
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(Seurat))
suppressMessages(library(ggpubr))
suppressMessages(library(ggstatsplot))
suppressMessages(library(patchwork))
suppressMessages(library(forcats))
suppressMessages(library(showtext))
suppressMessages(library(SCP))


inputDir = "results/"
expID = "NSCLCTUMOR_Public"
set.seed(1)

dir.create(file.path(inputDir, expID), showWarnings = FALSE)
expDir <- paste(inputDir, expID, sep = "/")
clstFolder <- paste(expID, "cluster_260316", sep = "_")
dir.create(file.path(expID, clstFolder), showWarnings = FALSE)
clstDir <- paste(expDir, clstFolder, sep = "/")
tmp_pf <- paste(clstDir, "/", clstFolder, "_", sep = "")

####Load data----------------
useObj = readRDS(file = "data/Public_late+early_MergeObj_260130.RDS")
table(useObj$stage)

DefaultAssay(useObj) <- "RNA"
df <- readxl::read_xlsx("Lodi_CRM_2025_signature/mono+mac_fuction_use.xlsx")

column_list <- lapply(df, function(x) x[!is.na(x)])
print(column_list)
print(names(column_list))

ignored_gene_sets <- list()

for (iscore in names(column_list)) {
  print(paste("Adding module score for:", iscore))
  
  gene_set <- column_list[[iscore]]
  
  existing_genes <- intersect(gene_set, rownames(useObj))
  missing_genes <- setdiff(gene_set, rownames(useObj))
  
  if (length(missing_genes) > 0) {
    ignored_gene_sets[[iscore]] <- missing_genes
  }
  
  if (length(existing_genes) == 0) {
    print(paste("No genes found in object for:", iscore))
    next  
  }
  
  existing_genes <- as.data.frame(existing_genes)
  colnames(existing_genes) <- iscore
  existing_genes <- as.list(existing_genes)
  
  useObj <- AddModuleScore(object = useObj,
                           features = existing_genes,
                           ctrl = 5,
                           name = iscore,
                           assay = "RNA")
}

if (length(ignored_gene_sets) > 0) {
  print("Ignored gene sets:")
  print(ignored_gene_sets)
  save(ignored_gene_sets,file = "Ignored_gene_set.Rdata")
} else {
  print("No gene sets were ignored.")
}

colnames(useObj@meta.data)
Idents(useObj) <- useObj$stage

Function <- c("Anti.inflammatory1")

####Fig.S3G-------------------
col_use <- c("#80bf34","#ffd17a")
VlnGG = FeatureStatPlot(useObj,stat.by = Function,group.by = "stage",add_box = TRUE,bg_palcolor = "white",
                        palcolor = col_use,comparisons = list(c("Early","Late")),pairwise_method = "wilcox.test",ncol = 6)&NoLegend()
ggsave(filename = paste(tmp_pf,"addscore_violin_p.sign.png"),plot = VlnGG,width = 18,height = 12,units = "in",dpi = 600)
ggsave(filename = paste(tmp_pf,"addscore_violin_p.sign.pdf"),plot = VlnGG,width = 18,height = 12,units = "in")

