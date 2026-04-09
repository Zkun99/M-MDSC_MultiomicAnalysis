####################identification of M-MDSC subpopulations############
#proportion_list：the proportion of different types of cells in the 20 cell gene neighborhood

#cluster M-MDSC based on the neighborhood
mat <- proportion_list[rownames(sem2_log_norm_rm@meta.data)[sem2_log_norm_rm$Cluster_Anno_new=="MDSC"],c(1:24)]

set.seed(123)
p=Heatmap(pheatmap:::scale_rows(mat),show_row_dend = FALSE,show_row_names = FALSE,row_km = 8,show_column_dend = FALSE,col =colorRamp2(c(-2, 0, 2), hcl_palette = "Blue-Red"),border = TRUE)
set.seed(123)
order_mdsc = row_order(p)

###################transcriptional programs of the M-MDSC subpopulations######
library(TCseq)
result_list <- list()
#calculate the average expression of each M-MDSC subpopulation
for(things in names(order_mdsc)){
    temp_list=data.frame(order_mdsc[[things]])[,1]
    rowmean=rowMeans(sem2_log_norm_rm[["RNA"]]@data[,temp_list])
    result_list[[things]]=rowmean
}
temp_val=names(result_list)
result_list=as.data.frame(do.call(cbind,result_list))
colnames(result_list)=temp_val

set.seed(666)
aaa=timeclust(as.matrix(result_list), algo = "km", k = 8, standardize = TRUE,dist="distance",dist.method="euclidean")

clust_var=data.frame(aaa@cluster)
clust_var$genes=rownames(clust_var)
clust_var=clust_var[order(clust_var[,1]),]

