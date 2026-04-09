rm(list = ls())
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

inputDir <- "6CloupeFile/results/"
expID <- "LCCHEMOICBPBMC_12_SINGLET_230110"

dir.create(file.path(inputDir, expID), showWarnings = FALSE)
expDir <- paste(inputDir, expID, sep = "/")
ppf <- paste(expDir, "/", expID, "_", sep = "")

clstFolder <- paste(expID, "cluster_01", sep = "_")
dir.create(file.path(expDir, clstFolder), showWarnings = FALSE)
clstDir <- paste(expDir, clstFolder, sep = "/")
tmp_pf <- paste(clstDir, "/", clstFolder, "_", sep = "")


main_dir <- "NSCLC_PBMC_scRNA_seq/scrublet"
sub_dirs <- list.dirs(main_dir, recursive = FALSE)
all_data <- data.frame()
for (i in seq_along(sub_dirs)) {
  sub_dir <- sub_dirs[i]

  file_path <- file.path(sub_dir, "doublet.txt")

  if (file.exists(file_path)) {
    data <- read.csv(file_path)
    
    folder_name <- basename(sub_dir)
    
    data$barcode <- paste0(data$barcode, "_", i)
    
    data$Folder <- folder_name
    
    all_data <- rbind(all_data, data)
  } else {
    warning(paste("File not found:", file_path))
  }
}

head(all_data)
write.csv(all_data,file = paste(tmp_pf,"doublet_scrublet.csv",sep = ""))
all_data <- read.csv(file = paste(tmp_pf,"doublet_scrublet.csv",sep = ""))

all_data <- all_data[!duplicated(all_data$barcode),]

rownames(all_data) <- all_data$barcode

all_data <- all_data[,c("doublet_scores","predicted_doublets")]

###MDSC doublet cell proportions----------
late_PBMC_MDSC <- readRDS("MDSC_0526/LCCHEMOICBPBMC_12_SINGLET_cluster_230110_MDSC_res0.3.rds")

doublet_find <- late_PBMC_MDSC
doublet_find <- AddMetaData(doublet_find,metadata = all_data)

table(doublet_find$anno,doublet_find$predicted_doublets)

cat("Cell proportions\n")
group_ratio <- doublet_find@meta.data %>%
  group_by(anno, predicted_doublets) %>%
  summarize(n = n())
group_ratio$predicted_doublets <- factor(group_ratio$predicted_doublets,levels = c("True","False"))

plot_sample <- ggplot(group_ratio,aes(x = anno,weight = n,fill = predicted_doublets))+
  geom_bar(position = "fill")+
  coord_flip()+
  scale_fill_manual(values = c("#881400","#bbb7b7")) + ##69b2cf
  theme_scp()+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 0.5),
        axis.ticks = element_line(colour = "black"),
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.1, size =16))+
  labs(y = "Proportion")+
  scale_y_continuous(expand = c(0,0.02))
plot_sample
ggsave(filename = paste(tmp_pf,"doublet_percentage_barplot_flip_MDSC.png",sep = ""),plot = plot_sample,width = 7,height = 2,dpi = 600,units = "in") 
ggsave(filename = paste(tmp_pf,"doublet_percentage_barplot_flip_MDSC.pdf",sep = ""),plot = plot_sample,width = 7,height = 2,units = "in") 
