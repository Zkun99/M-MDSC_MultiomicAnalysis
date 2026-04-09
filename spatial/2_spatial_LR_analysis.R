###########################Direct contact formula################################
replace_genes <- function(gvec) {
    replacements <- c("IFNA4" = "IFNA4/10/17", "IFNA10" = "IFNA4/10/17", "IFNA17" = "IFNA4/10/17",
                      "KLRC1" = "KLRC1/2", "KLRC2" = "KLRC1/2", "LILRB3" = "LILRB3/A6",
                      "TNXB" = "TNXA/B", "XCL1" = "XCL1/2", "XCL2" = "XCL1/2")
    #replace the L/R information from CellChat to genes in SMI panel
    gvec[gvec %in% names(replacements)] <- replacements[gvec[gvec %in% names(replacements)]]
    return(gvec)
}

score_jie <- function(genes, mat) {
    vals <- rowMeans(mat[genes, , drop = FALSE])
    if (length(genes) > 1) prod(vals)^(1 / length(genes)) else vals
}
#score_jie: average expression

gene_nonzero_count <- function(genes, mat) {
    nonzero <- rowSums(mat[genes, , drop = FALSE] > 0)
    if (length(genes) == 1) as.numeric(nonzero) else min(nonzero)
}
#gene_nonzero_count: number of non-zero expression cells

direct_contact_lr <- function(lr_info2,expr_mat,expr_mdsc_near,expr_endo_near,count_mdsc_near,count_endo_near) {

	#function for direct contact
    lr_scores_1= numeric(nrow(lr_info2))
	lr_scores_2= numeric(nrow(lr_info2))
    for (i in seq_len(nrow(lr_info2))) {
	    temp_ligand <- strsplit(lr_info2[i, "ligand_symbol"], "&")[[1]]
	    temp_receptor <- strsplit(lr_info2[i, "receptor_symbol"], "&")[[1]]
	    
	    temp_ligand <- intersect(rownames(expr_mat), replace_genes(temp_ligand))
	    temp_receptor <- intersect(rownames(expr_mat), replace_genes(temp_receptor))

	    if (length(temp_ligand) > 0 && length(temp_receptor) > 0) {
	        # direction 1：
	        aa <- score_jie(temp_ligand, expr_endo_near)
	        bb <- score_jie(temp_receptor, expr_mdsc_near)
	        aa_pct <- gene_nonzero_count(temp_ligand, count_endo_near) / ncol(expr_endo_near)
	        bb_pct <- gene_nonzero_count(temp_receptor, count_mdsc_near) / ncol(expr_mdsc_near)
	        
	        # direction 2：
	        cc <- score_jie(temp_ligand, expr_mdsc_near)
	        dd <- score_jie(temp_receptor, expr_endo_near)
	        cc_pct <- gene_nonzero_count(temp_ligand, count_mdsc_near) / ncol(expr_mdsc_near)
	        dd_pct <- gene_nonzero_count(temp_receptor, count_endo_near) / ncol(expr_endo_near)
	        
	        score1 <- (aa * bb * aa_pct * bb_pct)^(1/4)
	        score2 <- (cc * dd * cc_pct * dd_pct)^(1/4)
	        
	        lr_scores_1[i] <- score1
	        lr_scores_2[i] <- score2
	    }
	    print(i)
	}

	result_temp1=cbind(lr_info2,lr_scores_1,lr_scores_2)
	result_temp1=result_temp1[order(result_temp1$lr_scores_1,decreasing = TRUE),]
	result_temp1=result_temp1[result_temp1$lr_scores_1!=0 & result_temp1$lr_scores_2!=0,]
	result_temp1$rank_lr1 <- rank(-result_temp1$lr_scores_1, ties.method = "min")
	result_temp1$rank_lr2 <- rank(-result_temp1$lr_scores_2, ties.method = "min")
	return(result_temp1)
}

############################Indirect contact formula#####################################
get_top_cell_by_genes <- function(gene_list, expr_mat) {
  valid_genes <- intersect(gene_list, rownames(expr_mat))
  if (length(valid_genes) == 0) return(NA)
  
  avg_expr <- colMeans(expr_mat[valid_genes, , drop = FALSE])
  
  top_cell <- names(avg_expr)[which.max(avg_expr)]
  return(top_cell)
}

indirect_contact_lr <- function(lr_info2,expr_mat,expr_mdsc_notnear,expr_endo_notnear,count_mdsc_notnear,count_endo_notnear) {

	lr_scores_1= numeric(nrow(lr_info2))
	lr_scores_2= numeric(nrow(lr_info2))
	for (i in seq_len(nrow(lr_info2))) {
	    temp_ligand <- strsplit(lr_info2[i, "ligand_symbol"], "&")[[1]]
	    temp_receptor <- strsplit(lr_info2[i, "receptor_symbol"], "&")[[1]]
	    
	    temp_ligand <- intersect(rownames(expr_mat), replace_genes(temp_ligand))
	    temp_receptor <- intersect(rownames(expr_mat), replace_genes(temp_receptor))

	    if (length(temp_ligand) > 0 && length(temp_receptor) > 0) {
	        # direction 1
	        aa <- score_jie(temp_ligand, expr_endo_notnear)
	        bb <- score_jie(temp_receptor, expr_mdsc_notnear)
	        aa_pct <- gene_nonzero_count(temp_ligand, count_endo_notnear) / ncol(expr_endo_notnear)
	        bb_pct <- gene_nonzero_count(temp_receptor, count_mdsc_notnear) / ncol(expr_mdsc_notnear)
	        
	        # direction 2
	        cc <- score_jie(temp_ligand, expr_mdsc_notnear)
	        dd <- score_jie(temp_receptor, expr_endo_notnear)
	        cc_pct <- gene_nonzero_count(temp_ligand, count_mdsc_notnear) / ncol(expr_mdsc_notnear)
	        dd_pct <- gene_nonzero_count(temp_receptor, count_endo_notnear) / ncol(expr_endo_notnear)
	        
	        first_dist=st_distance(result[!is.na(match(result[[1]],get_top_cell_by_genes(temp_ligand,expr_endo_notnear))),],result[!is.na(match(result[[1]],get_top_cell_by_genes(temp_receptor,expr_mdsc_notnear))),])
	        second_dist=st_distance(result[!is.na(match(result[[1]],get_top_cell_by_genes(temp_ligand,expr_mdsc_notnear))),],result[!is.na(match(result[[1]],get_top_cell_by_genes(temp_receptor,expr_endo_notnear))),])
	        score1 <- ((aa * bb * aa_pct * bb_pct)^(1/4))/first_dist
	        score2 <- ((cc * dd * cc_pct * dd_pct)^(1/4))/second_dist

	        lr_scores_1[i] <- score1
	        lr_scores_2[i] <- score2
	        print(i)
	    } 
	}

	result_temp2=cbind(lr_info2,lr_scores_1,lr_scores_2)
	result_temp2=result_temp2[order(result_temp2$lr_scores_1,decreasing = TRUE),]
	result_temp2=result_temp2[result_temp2$lr_scores_1!=0 & result_temp2$lr_scores_2!=0,]
	result_temp2$rank_lr1 <- rank(-result_temp2$lr_scores_1, ties.method = "min")
	result_temp2$rank_lr2 <- rank(-result_temp2$lr_scores_2, ties.method = "min")
	return(result_temp2)
}

#Example calling of the function:
#read the LR information from CellChat
#lr_info2=read.csv("./6_LR_data/LR_info_adjusted.csv",header=TRUE)

#sem2_log_norm_rm is the main object of spatial transcriptomics
#expr_mat <- sem2_log_norm_rm[["RNA"]]@data
#count_mat <- sem2_log_norm_rm[["RNA"]]@counts

#endo_near is the cell list of endothelial cells; endo_mdsc_withnear is the cell list of MDSCs
#expr_endo_near <- expr_mat[, endo_near]
#expr_mdsc_near <- expr_mat[, endo_mdsc_withnear]
#count_endo_near <- count_mat[, endo_near]
#count_mdsc_near <- count_mat[, endo_mdsc_withnear]
#result_temp1=direct_contact_lr(lr_info2,expr_mat,expr_mdsc_near,expr_endo_near,count_mdsc_near,count_endo_near)

#Example output (or return of the function) for spatial L/R analysis:
#"pathway_name"	"ligand"	"receptor"	"annotation"	"ligand_symbol"	"receptor_symbol"	"lr_scores_1"	"lr_scores_2"	"rank_lr1"	"rank_lr2"	"pair"	"difference"
#"FN1"	"FN1"	"CD44"	"ECM-Receptor"	"FN1"	"CD44"	0.299103354247731	0.432628590320514	5	5	"FN1-CD44"	0
#"COLLAGEN"	"COL6A2"	"CD44"	"ECM-Receptor"	"COL6A2"	"CD44"	0.251348429905784	0.350540563476509	12	12	"COL6A2-CD44"	0
#"COLLAGEN"	"COL6A3"	"CD44"	"ECM-Receptor"	"COL6A3"	"CD44"	0.254650033426687	0.372003249346408	11	10	"COL6A3-CD44"	1