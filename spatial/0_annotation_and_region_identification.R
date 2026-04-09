library(Seurat)
library(plyr)
library(dplyr)

################################(1) Quality filtering #######################################
#sem is the object from AtoMx platform
#<1> quality filtering based on nCount, nFeature, area size, and qcFlag
sem <- AddMetaData(object = sem,     
                    metadata = as.numeric(sem@meta.data$nCount_RNA>=50),   
                    col.name = "nCount_larger_than_50")
sem <- AddMetaData(object = sem,    
                    metadata = as.numeric(sem@meta.data$nFeature_RNA>=10), 
                    col.name = "nFeature_larger_than_10")
#genes expressed in fewer than 10 cells
table(data.frame(rowSums(sem[["RNA"]]@counts>1))[,1]>10)
#6175, all genes expressed in more than 10 cells
sem <- AddMetaData(object = sem,     #seurat对象
                    metadata = as.numeric(sem@meta.data$Area<5*mean(sem@meta.data$Area)), 
                    col.name = "cell_area_5mean")
sem <- AddMetaData(object = sem, 
                    metadata = data.frame(colSums(sem[["RNA"]]$counts))[,1], 
                    col.name = "library_size")

#Filtering of the cells
sem1=subset(sem,cells = rownames(sem@meta.data)[data.frame(rowSums(sem@meta.data[,c("cell_area_5mean","nCount_larger_than_50","nFeature_larger_than_10")])==3)[,1]])
sem2=subset(sem1,cells=rownames(sem1@meta.data)[sem1@meta.data$qcFlagsCellCounts=="Pass" & sem1@meta.data$qcFlagsCellPropNeg=="Pass" & sem1@meta.data$qcFlagsCellComplex=="Pass" & as.numeric(sem1@meta.data$qcCellsFlagged)==0])

###############################(2) Data Normalization ######################################
#We conducted log normalization, total count normalization, and SCTransform (regression out library size and cell area). We selected log normalization as it generated an UMAP with higher discrimination between cell types.

sem2_log_norm=NormalizeData(sem2)
sem2_log_norm <- FindVariableFeatures(sem2_log_norm, selection.method = "vst", nfeatures = 2000)
sem2_log_norm <- ScaleData(sem2_log_norm, verbose = FALSE)
sem2_log_norm <- RunPCA(sem2_log_norm, npcs = 30, verbose = FALSE)
ElbowPlot(sem2_log_norm,ndims=30)
#threshold is 20
sem2_log_norm <- FindNeighbors(object = sem2_log_norm, dims = 1:20)
for(res in seq){
    sem2_log_norm <- FindClusters(sem2_log_norm, resolution = res)
}
sem2_log_norm <- RunUMAP(sem2_log_norm, dims = 1:20)

clustree(sem2_log_norm@meta.data,prefix="RNA_snn_res.")
#resolution=0.3
DimPlot(object = sem2_log_norm, reduction = "umap",group.by = "RNA_snn_res.0.3",label = TRUE)



####################################(3) Remove cells at FOV boudaries#####################################
fov_annotation=read.table("../flatFiles/215548240653247574215548/215548240653247574215548_fov_positions_file.csv",header=TRUE,sep=",")

var_fov_edge=matrix(0,nrow=dim(sem2_log_norm@meta.data),ncol=1)
rownames(var_fov_edge)=rownames(sem2_log_norm@meta.data)

for(eles in unique(sem2_log_norm@meta.data$fov)){

    for(things in rownames(sem2_log_norm@meta.data)[sem2_log_norm@meta.data$fov==eles]){
        temp_inside=round(bbox(sem2_log_norm@images$X215548.240653.247574.215548@boundaries$segmentation@polygons[[things]])*1000/0.12028)
        fov_info=fov_annotation[eles,]

        if(abs(temp_inside[1,1]-fov_info[,2])<5 || abs(temp_inside[1,2]-fov_info[,2]-4256)<5 || abs(temp_inside[2,1]-fov_info[,3]+4256)<5 || abs(temp_inside[2,2]-fov_info[,3])<5){
            var_fov_edge[things,1]=1
            print(things)
        }
        
    }
    print(eles)
}

unique(match(rownames(var_fov_edge),rownames(sem2_log_norm@meta.data))==1:dim(sem2_log_norm@meta.data)[1])

#remove cells
sem2_log_norm_rm=subset(sem2_log_norm,cells = rownames(sem2_log_norm@meta.data)[var_fov_edge==0])
sem2_total_norm_rm=subset(sem2_total_norm,cells = rownames(sem2_total_norm@meta.data)[var_fov_edge==0])
sem2_sc_norm_rm=subset(sem2_sctransform,cells = rownames(sem2_sctransform@meta.data)[var_fov_edge==0])

##############################(4) Annotation using InsituType and custom methods#############
var1=as.matrix(sem2_log_norm[["RNA"]]$counts)
var1=t(var1)
var1=as.matrix(var1)
#var2: negative probe mean control value
var2=Matrix::colMeans(sem2_log_norm@assays$negprobes$counts)
#var3: reference matrix
first_level=read.table("./0_insitutype_anno_ref/250107_Public_Late_Tumor_RAW_for_SMI_anno_with_nanostring.csv",sep=",",row.names=1,header=TRUE)
first_level=as.matrix(first_level)


#<1> define the cohorts on image based on panCK and CD45 staining results from CosMx
cohort=read.table("../Napari_Stitched_Files/_metadata_FISH.csv",header=TRUE,sep=",",row.names = 1)
#Cohort definition is left to the user, and many approaches are reasonable. For example, you might:
#Define cohort by gating cells’ immunofluorescence stains, e.g. PanCK+/- or CD45+/-
#Cluster cells on their continuous immunofluorescence values
#Use an autoencoder to extract features from cells’ images, then cluster on those features.
#Use the results of spatial clustering / “niches” to define cohorts
#Cluster on information from multiple data types together, such as the average gene expression profile of cell’s neighbors AND immunofluorescence.

#<2> Semi-supervised clustering
semisup <- insitutype(
  x = var1,
  neg = var2,
  cohort = as.vector(cohort[rownames(sem2_log_norm@meta.data),1]),
  # Enter your own per-cell background estimates here if you
  # have them; otherwise insitutype will use the negprobes to
  # estimate background for you.
  bg = NULL,
  # condensed to save time. n_clusts = 5:15 would be more optimal
  n_clusts = 0, #pure supervised
  reference_profiles = first_level,
  update_reference_profiles = TRUE,
  # choosing inadvisably low numbers to speed the vignette; using the defaults
  # in recommended.
  n_phase1 = 5000,
  n_phase2 = 20000,
  n_phase3 = 100000,
  n_starts = 20,
  max_iters = 5,
  min_anchor_llr = 0.03, 
  insufficient_anchors_thresh = 5,
  n_anchor_cells=100,
  min_anchor_cosine = 0.3
) 

#Structure of semisup:
# a vector of cells’ cluster assignments: $clust
# a vector of the posterior probability for each cell typing call: $prob
# a matrix of cells’ log-likelihoods under each cluster: $logliks
# a matrix of the mean background-subtracted expression profile of each cluster: $profiles

#<2> generate a flightpath plot:
flightpath_plot(flightpath_result = NULL, insitutype_result = semisup)

#<3> we next conducted multiple rounds of re-clustering on cell types that does not fit into clustering results.

####################################(5) Region identification on slides##################################
#<1> tumor region definition
#result is an sf object (which could be operated by sf package) and stores the boundaries of each cell, result$expanced are the boundaries expanded by 3px

#first, extract the boundaries of tumor cells, and store in ttt variable
library(igraph)
temp=intersect(result[[1]],rownames(sem2_total_norm_rm@meta.data)[sem2_total_norm_rm$Cluster_Anno_new %in% c("Epithelial (Club Cell)","Epithelial (Ciliated Cell)")])
ttt=result[!is.na(match(result[[1]],temp)),]

#find the intersections of tumor cells (after expanding the boundaries by 3px)
intersections=st_intersects(ttt$expanced)

#use graph object to find the clusters 
g <- graph_from_adj_list(intersections)
ttt$cluster <- components(g)$membership

ttt2=st_drop_geometry(ttt)
temp=data.frame(table(ttt2[,2]))
ttt2[,2]=as.character(ttt2[,2])
cluster_list=ttt2[grep(paste("^",paste(rownames(temp)[temp[,2]>10],collapse="$|^"),"$",sep=""),ttt2[,2]),1]

#temp_list: tumor cells that with more than 15 tumor cells around
#epi_data stores the neighborhood of tumor cells
temp_list=rownames(epi_data)[rowSums(epi_data!=0)>=15]

#ttt3 stores individual cells in tumor region
ttt3=result[!is.na(match(result[[1]],intersect(temp_list,cluster_list))),]

#aggregate near regions from ttt3
split_list <- ttt3 %>% group_split(expanced.sample)
sample_names <- ttt3 %>% group_keys(expanced.sample) %>% pull()

merged_list <- map2(split_list, sample_names, ~ {
  merged_geom <- st_union(st_buffer(.x, dist = 20))       # 合并相邻区域
  merged_sf <- st_as_sf(st_cast(merged_geom, "POLYGON"))  # 拆成 polygon
  merged_sf$sample <- .y                                   # 添加 sample 列
  merged_sf
})

#merged_sf_all: the final tumor region
merged_sf_all <- bind_rows(merged_list)
merged_sf_all$ID <- seq_len(nrow(merged_sf_all))

#generate the tumor-proximal regions
library(sf)
library(dplyr)
library(purrr)
library(terra)

ttt=st_buffer(merged_sf_all, dist = 416.6667*3)
ttt=st_union(ttt)
ttt2=st_union(merged_sf_all)
round_info=st_difference(ttt,ttt2)


############################endothelial cell aggregate region##############################
epi_data=neighbors_nearest_50[sem2_total_norm_rm$Cluster_Anno_new %in% c("Endothelial Cell"),sem2_total_norm_rm$Cluster_Anno_new %in% c("Endothelial Cell")]
temp_list=rownames(epi_data)[rowSums(epi_data!=0)>=5]
endo_region=result[!is.na(match(result[[1]],temp_list)),]
endo_region2=st_union(st_buffer(endo_region, dist = 1))

##################################unsupervised clustering on immune cells ###############
#Identify Immune-cell aggregated regions
#select immune cell types
immune_cells=c("DC","CD8..cytotoxic.T.cell","NK Cell","Neutrophil","MDSC","Mast Cell","Lymphocyte","Alveolar Macrophage","Treg","CD4..T.cell","Monocyte","Macrophage","B Cell","Kappa Plasma Cell","Lambda Plasma Cell")

library(dbscan)
library(dplyr)
library(ggplot2)
ttt=data.frame(sem2_log_norm@meta.data[,c("cell_ID","x_slide_mm","y_slide_mm","Cluster_Anno_new")])
immune_df_filtered <- ttt %>%
  filter(Cluster_Anno_new %in% immune_cells)

db <- dbscan(immune_df_filtered[, c("x_slide_mm", "y_slide_mm")], eps = 0.05,minPts=50)
immune_df_filtered$cluster <- db$cluster

immune_df_filtered <- immune_df_filtered %>%
  group_by(cluster) %>%
  mutate(cluster_size = n(),
         cluster = ifelse(cluster_size < 100, NA, cluster)) %>%
  ungroup() %>%
  select(-cluster_size)  

immune_df_filtered$cluster[immune_df_filtered$cluster=="0"]=NA
immune_df_filtered=na.omit(immune_df_filtered)

#the annotation of the immune regions was conducted by the 

#immune cell numbers in each immune-cell aggregated regions
composition <- immune_df_filtered %>%
  filter(cluster != 0) %>%
  group_by(cluster, Cluster_Anno_new) %>%
  tally() %>%
  tidyr::pivot_wider(names_from = Cluster_Anno_new, values_from = n, values_fill = 0)

#percentage
composition_pct <- immune_df_filtered %>%
    filter(cluster != 0) %>%
    group_by(cluster, Cluster_Anno_new) %>%
    tally() %>%
    tidyr::pivot_wider(names_from = Cluster_Anno_new, values_from = n, values_fill = 0) %>%
    mutate(total = rowSums(across(where(is.numeric)))) %>%
    mutate(across(where(is.numeric) & !matches("total"), ~ round(.x / total * 100, 2)))  
composition_pct=data.frame(composition_pct)
rownames(composition_pct)=composition_pct$cluster

#annotate the regions by unsupervised clustering
set.seed(123)
p=Heatmap(pheatmap:::scale_rows(composition_pct[,2:16]),row_km = 5,show_row_dend = FALSE,show_column_dend = FALSE,border=TRUE,col =colorRamp2(c(-2, 0, 2), hcl_palette = "Blue-Red"),show_row_names = FALSE)
immune_region_backup=p

row_clusters=row_order(p)
heatmap_cell_ids <- rownames(composition_pct)
library(tibble)
cluster_map <- bind_rows(lapply(seq_along(row_clusters), function(i) {
  ids <- heatmap_cell_ids[row_clusters[[i]]]
  tibble(cluster = ids, heatmap_cluster = paste0("Group_", i))
}))
cluster_map$cluster=as.integer(cluster_map$cluster)
cluster_map=data.frame(cluster_map)

immune_df_filtered2 <- immune_df_filtered %>%
  left_join(cluster_map, by = "cluster")

immune_df_filtered2$heatmap_cluster <- gsub("Group_1", "Treg-enriched", immune_df_filtered2$heatmap_cluster)
immune_df_filtered2$heatmap_cluster <- gsub("Group_2", "Killing-associated", immune_df_filtered2$heatmap_cluster)
immune_df_filtered2$heatmap_cluster <- gsub("Group_3", "TLS-like", immune_df_filtered2$heatmap_cluster)
immune_df_filtered2$heatmap_cluster <- gsub("Group_4", "MDSC-enriched", immune_df_filtered2$heatmap_cluster)
immune_df_filtered2$heatmap_cluster <- gsub("Group_5", "Plasma cell-enriched", immune_df_filtered2$heatmap_cluster)

immune_regions=result[!is.na(match(result[[1]],immune_df_filtered2$cell_ID)),]