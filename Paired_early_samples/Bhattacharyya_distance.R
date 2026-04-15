# Pair inhouse early-stage PBM+Tumor Bhattacharyya Distance

rm(list = ls())

suppressMessages(library(Seurat))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(distdimscr))
suppressMessages(library(glmnet))
suppressMessages(library(ggbreak))

set.seed(1)
inputDir <- "results/"

expID <- "LC_inhouse_Paired_MDSC_250527"

clusterFlag <- TRUE
rnaInteMethod <- "RPCA"
adtInteMethod <- "RPCA"

dir.create(file.path(inputDir, expID), showWarnings = FALSE)
expDir <- paste(inputDir, expID, sep = "/")
ppf <- paste(expDir, "/", expID, "_", sep = "")

clstFolder <- paste(expID, "cluster", sep = "_")
dir.create(file.path(expDir, clstFolder), showWarnings = FALSE)
clstDir <- paste(expDir, clstFolder, sep = "/")
tmp_pf <- paste(clstDir, "/", clstFolder, "_", sep = "")

#######Function1:Statistical analysis for Bhattacharyya distance results-------------------
#' @param bhatt_results Dataframe from bhastdis function
#' @param group_name Name of the comparison group
#' @param test_type Type of statistical test ("wilcoxon" or "t_test")
#'
#' @return A list with statistical results
#' @export
#'
analyze_bhatt_distance <- function(bhatt_results, group_name, test_type = "wilcoxon") {
  
  # Extract real and random distances
  real_dist <- bhatt_results$MDSC.cells.distance[bhatt_results$comparison == "real"]
  random_dist <- bhatt_results$MDSC.cells.distance[bhatt_results$comparison == "random"]
  
  # Perform statistical test
  if (test_type == "wilcoxon") {
    # One-sided Wilcoxon rank-sum test (real > random)
    test_result <- wilcox.test(real_dist, random_dist, 
                               alternative = "greater", 
                               exact = FALSE)
  } else if (test_type == "t_test") {
    # One-sided t-test
    test_result <- t.test(real_dist, random_dist, 
                          alternative = "greater")
  }
  
  # Calculate fold change and other statistics
  stats_summary <- list(
    group = group_name,
    test_type = test_type,
    p_value = test_result$p.value,
    fold_change = median(real_dist) / median(random_dist),
    mean_fold_change = mean(real_dist) / mean(random_dist),
    real_median = median(real_dist),
    random_median = median(random_dist),
    real_mean = mean(real_dist),
    random_mean = mean(random_dist),
    real_sd = sd(real_dist),
    random_sd = sd(random_dist),
    n_iterations = length(real_dist)
  )
  
  return(stats_summary)
}

#' Generate comprehensive statistical report
#'
#' @param stats_list List of statistical results from analyze_bhatt_distance
#'
#' @return A formatted dataframe with all statistical results
#' @export
#'
generate_stats_report <- function(stats_list) {
  report_df <- data.frame(
    Group = sapply(stats_list, function(x) x$group),
    P_value = sapply(stats_list, function(x) x$p_value),
    Fold_Change = sapply(stats_list, function(x) x$fold_change),
    Real_Median = sapply(stats_list, function(x) x$real_median),
    Random_Median = sapply(stats_list, function(x) x$random_median),
    Real_Mean = sapply(stats_list, function(x) x$real_mean),
    Random_Mean = sapply(stats_list, function(x) x$random_mean),
    Real_SD = sapply(stats_list, function(x) x$real_sd),
    Random_SD = sapply(stats_list, function(x) x$random_sd),
    N_Iterations = sapply(stats_list, function(x) x$n_iterations),
    Significance = sapply(stats_list, function(x) {
      ifelse(x$p_value < 0.001, "***",
             ifelse(x$p_value < 0.01, "**",
                    ifelse(x$p_value < 0.05, "*", "ns")))
    })
  )
  
  # Format p-values for better display
  report_df$P_value_formatted <- sapply(report_df$P_value, function(p) {
    ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))
  })
  
  return(report_df)
}

#####Function2:Enhanced visualization with statistical annotations----------
#' @param bhatt_data Combined bhatt distance data
#' @param stats_report Statistical report from generate_stats_report
#'
#' @return A ggplot object with statistical annotations
#' @export
#'
plot_enhanced_bhatt <- function(bhatt_data, stats_report) {
  
  p <- ggplot(data = bhatt_data) +
    geom_boxplot(aes(x = comparison, y = MDSC.cells.distance), 
                 col = "black", outlier.shape = NA, fill = "lightblue", alpha = 0.7) +
    geom_jitter(aes(x = comparison, y = MDSC.cells.distance, color = MDSC.cells.distance), 
                size = 0.8, shape = 21, alpha = 0.6, width = 0.2) +
    scale_color_distiller(palette = "Blues", guide = "none") +
    theme_bw() +
    ggtitle("Pair early-stage MDSC vs. Macrophage") +
    ylab("Bhattacharrya distance") +
    theme(
      axis.text.x = element_text(color = "black", family = "sans", angle = 30, 
                                 vjust = 0.85, hjust = 0.75),
      axis.text.y = element_text(color = "black", family = "sans"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", linetype = 1, color = "black", size = 0.25),
      strip.text = element_text(size = 10, face = "bold")
    ) +
    facet_wrap(~ group, nrow = 1) +
    scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5)) +
    scale_x_discrete(limits = c("real", "random"))
  

  annotation_data <- data.frame(
    group = stats_report$Group,
    x = 1.5,  
    y = 2.8,  
    label = paste0("P ", stats_report$P_value_formatted, " ", stats_report$Significance,
                   "\nFC = ", round(stats_report$Fold_Change, 2))
  )
  
  p <- p + geom_text(data = annotation_data, 
                     aes(x = x, y = y, label = label),
                     size = 3, hjust = 0.5, vjust = 1,
                     inherit.aes = FALSE)
  
  return(p)
}
####Load data------------------
overallobj <- readRDS(file = paste(tmp_pf,"Seurat_Objects_Clustered.RDS", sep = ""))
overallobj$major_type <- ifelse(overallobj$celltype%in%c("MDSC"),"MDSC-PBMC",
                                ifelse(overallobj$celltype%in%c("M-MDSC"),"MDSC-Tumor","Macro-Tumor"))
table(overallobj$major_type)

overall.metadata <- as.data.frame(overallobj@meta.data)
overall.metadata <- overall.metadata[,c("orig.ident","major_type","tissue","MDSC_substype")]

overall.umap <- overallobj@reductions$umap@cell.embeddings %>%  as.data.frame() 
overall.data <- cbind(overall.umap,overall.metadata)
ggplot(overall.data,aes(x=UMAP_1,y=UMAP_2,colour=major_type)) +
  geom_point() +
  theme_classic() +
  facet_wrap(~tissue)
ggsave(paste(tmp_pf,"PBMC_T_split_0.2.png", sep = ""),dpi = 600, height = 6, width = 12)

overall.pca <- as.data.frame(overallobj@reductions$pca@cell.embeddings)

save(overallobj,overall.metadata,overall.data,overall.umap,overall.pca,file = paste(tmp_pf,"merge_MDSC_ALL_Mac_Bhatta.RData", sep = ""))


##### reload rdata used in Bhatta distance ######
load(file = paste(tmp_pf,"merge_MDSC_ALL_Mac_Bhatta.RData", sep = ""))

print(overallobj[["pca"]], dims = 1:10, nfeatures = 5)
DimPlot(overallobj, reduction = "pca")
DimHeatmap(overallobj, dims = 1:15, cells = 500, balanced = TRUE)

overall.pca <- overall.pca[,c(1:11)]

DimPlot(overallobj, reduction = "pca",group.by = "orig.ident",label = T)

##########Compare circulating MDSC vs.TAMs and circulating MDSC vs. intratumoral MDSC-------------
nCell = 200
bhastdis = function(ctrl_pca,treat_pca){
  bhatt.dist <- bhatt.dist.rand <- vector("logical",length = nCell)
  for (i in 1:nCell) {
    bhatt.dist[[i]] <- dim_dist(embed_mat_x = ctrl_pca,embed_mat_y = treat_pca,dims_use = 1:10,
                                num_cells_sample = nCell,distance_metric = "bhatt_dist",random_sample = FALSE)
    bhatt.dist.rand[[i]] <- dim_dist(embed_mat_x = ctrl_pca,embed_mat_y = treat_pca,dims_use = 1:10,
                                     num_cells_sample = nCell,distance_metric = "bhatt_dist",random_sample = TRUE)
  }
  bhatt.dist <- data.frame(MDSC.cells.distance = bhatt.dist,comparison = "real")
  bhatt.dist.rand <- data.frame(MDSC.cells.distance = bhatt.dist.rand,comparison = "random")
  bhatt.res <- rbind(bhatt.dist,bhatt.dist.rand)
}
table(overall.metadata$major_type)


PBMC_MDSC <- rownames(overall.data)[overall.data$major_type=="MDSC-PBMC"]
Tumor_Macro <- rownames(overall.data)[overall.data$major_type=="Macro-Tumor"]
PBMC_MDSC_pca <- overall.pca[PBMC_MDSC,]
Tumor_Macro_pca <- overall.pca[Tumor_Macro,]

PBMC_Macro_dis = bhastdis(ctrl_pca = PBMC_MDSC_pca,treat_pca = Tumor_Macro_pca)
PBMC_Macro_dis$group = "PBMC_MDSC_Macro"

Tumor_MDSC <- rownames(overall.data)[overall.data$major_type=="MDSC-Tumor"]
Tumor_MDSC_pca <- overall.pca[Tumor_MDSC,]
MDSC_dis = bhastdis(ctrl_pca = PBMC_MDSC_pca,treat_pca = Tumor_MDSC_pca)
MDSC_dis$group = "Between_MDSC"

bhatt.all = rbind(PBMC_Macro_dis,MDSC_dis)
bhatt.all$comparison <- factor(bhatt.all$comparison,levels = c("real","random"))
write.csv(bhatt.all,file = paste(tmp_pf,"PBMC_MDSC_distance_PC1-11_random200.csv",sep = ""),row.names = FALSE)

####Fig.3H-------------
bhatt.all <- read.csv(file = paste(tmp_pf,"PBMC_MDSC_distance_PC1-11_random200.csv",sep = ""))

bhattGG = ggplot(data = bhatt.all)+
  geom_boxplot(aes(x = comparison,y = MDSC.cells.distance,col = comparison),col = "black",outlier.shape = NA)+
  geom_jitter(aes(x = comparison,y = MDSC.cells.distance,col = MDSC.cells.distance),size=0.5,shape = 21,alpha=0.5)+
  #scale_colour_gradient(low="#6600CC",high="gold")+
  scale_color_distiller(palette = "Blues")+#trans = "reverse"
  theme_bw()+
  ggtitle("Pair early-stage MDSC vs. Macrophage") +
  ylab("Bhattacharrya distance")+
  theme(axis.text.x = element_text(color = "black",family = "sans",angle = 30,vjust = 0.85,hjust = 0.75),
        axis.text.y = element_text(color = "black",family = "sans"),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",linetype=1,color="black",size=0.25))+
  facet_wrap(~group,nrow = 1)+
  scale_y_continuous(limits = c(0, 3),breaks = seq(0,3,0.1))+
  scale_y_break(c(0.3,2),
                space = 0.3,
                scales = 1.8)+
  scale_x_discrete(limits = c("real","random"))

bhattGG
ggsave(paste(tmp_pf,"PBMC_MDSC_distance_PC1-11_random200_251103.png", sep = ""),dpi = 600, height = 5, width = 4,units = "in")
ggsave(paste(tmp_pf,"PBMC_MDSC_distance_PC1-11_random200_251103.pdf", sep = ""),height = 5, width = 4,units = "in")

########## Statistical Analysis Section  ##########
PBMC_Macro_dis <- bhatt.all[bhatt.all$group%in%c("PBMC_MDSC_Macro"),]
MDSC_dis <- bhatt.all[bhatt.all$group%in%c("Between_MDSC"),]

PBMC_Macro_stats <- analyze_bhatt_distance(PBMC_Macro_dis, "PBMC_MDSC_Macro",test_type = "wilcoxon")
MDSC_stats <- analyze_bhatt_distance(MDSC_dis, "Between_MDSC", test_type = "wilcoxon")


stats_report <- generate_stats_report(list(PBMC_Macro_stats, MDSC_stats))
print("=== Bhattacharyya Distance Statistical Analysis ===")
print(stats_report)
write.csv(stats_report,file = paste(tmp_pf, "PBMC_MDSC_distance_statistical_analysis.csv", sep = ""), row.names = FALSE)


enhanced_plot <- plot_enhanced_bhatt(bhatt.all, stats_report)
enhanced_plot
ggsave(paste(tmp_pf, "PBMC_MDSC_distance_enhanced_stats_251103.png", sep = ""), dpi = 600, height = 5, width = 6, units = "in")
ggsave(paste(tmp_pf, "PBMC_MDSC_distance_enhanced_stats_251103.pdf", sep = ""), height = 5, width = 6, units = "in")

