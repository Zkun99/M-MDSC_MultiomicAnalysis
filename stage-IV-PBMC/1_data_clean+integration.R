################################

####1:Data clean----------------

################################

rm(list = ls())
suppressMessages(library(DoubletFinder))
suppressMessages(library(stringr))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(Seurat))

tenXDir <- "6CloupeFile/4Files/"
resDir <- "6CloupeFile/results/"

expID <- "lc_chemo_icb_pbmc_221202"

read_flag <- TRUE
intePip <- "standard"

sigMarker <- c("PTPRC", "EPCAM", "COL1A1", "CD3D", "CD8A", "CD4", "PDCD1", "TCF7", "NKG7", "MS4A1", "CD27", "CD19", 
	      "CD68", "FCER1A", "HLA-DRA", "NCAM1", "TPSAB1", "GNLY", "FOXP3", "TBX21", "GATA3")

dir.create(file.path(resDir, expID), showWarnings = FALSE)
expDir <- paste(resDir, expID, sep = "/")

if (read_flag) {
	cat("Start to read the single-cell data...\n")
	inFiles <- list.files(tenXDir, pattern = "filtered_feature_bc_matrix.h5")
	saName <- str_split_fixed(inFiles, "_fil", n = 2)[,1]
	smpTbl <- as.data.frame(str_split_fixed(saName, "_", n = 3)[,1:2])
	colnames(smpTbl) <- c("Time", "PID")
	rownames(smpTbl) <- saName
	print(head(smpTbl))

	cat("All the samples:", saName, "\n")
	cleanVar <- c("raw_cell_num", "raw_feature_num", "keep_empty", "keep_mito", "keep_ribo", "keep_feature", "doublet", "singlet")
	cleanRecord <- as.data.frame(matrix(ncol = length(cleanVar), nrow = length(inFiles)))
	cleanRecord$file <- inFiles
	cleanRecord$sample <- saName

	for (isa in 1:length(saName)) {
		st <- Sys.time()
		cat(saName[isa], ":" , inFiles[isa], "\n")
		tmpDir <- paste(tenXDir, inFiles[isa], sep = "/")
		tmpData <- Read10X_h5(tmpDir)
		print(names(tmpData))

		tmpObj <- CreateSeuratObject(counts = tmpData$`Gene Expression`, project = saName[isa], min.cells = 3, min.features = 5)
		cleanRecord[isa, "raw_cell_num"] <- ncol(tmpObj)
		cleanRecord[isa, "raw_feature_num"] <- nrow(tmpObj)
#		print(dim(tmpData$`Antibody Capture`))
#		print(head(colnames(tmpObj)))

		tmp_folder <- paste(saName[isa], expID, sep = "_")
		dir.create(file.path(expDir, tmp_folder), showWarnings = FALSE)
		tmp_pf <- paste(expDir, "/", tmp_folder, "/", saName[isa], "_", sep = "")

		cat("Mitochondrial/Riposomal QC\n")
		tmpObj <- PercentageFeatureSet(tmpObj, "^MT-", col.name = "percent_mito")
		tmpObj <- PercentageFeatureSet(tmpObj, "^RP[SL]", col.name = "percent_ribo")
		tmpObj <- PercentageFeatureSet(tmpObj, "^HB[^(P)]", col.name = "percent_hb")
		tmpObj <- PercentageFeatureSet(tmpObj, "PECAM1|PF4", col.name = "percent_plat")

		feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb", "percent_plat")
		qc_vln <- VlnPlot(tmpObj, features = feats, pt.size = 0.1, ncol = 3) + NoLegend()
		ggsave(paste(tmp_pf, "qc_vln.png", sep = ""), qc_vln, dpi = 300, height = 9, width = 9)

		print(tmpObj)

		cat("Filtering...\n")
		cleanObj <- tmpObj
		selected_c <- WhichCells(tmpObj, expression = nFeature_RNA > 200)
		cleanObj <- subset(cleanObj, cells = selected_c)
		cleanRecord[isa, "keep_empty"] <- length(selected_c)

		selected_mito <- WhichCells(cleanObj, expression = percent_mito < 25)
		cleanObj <- subset(cleanObj, cells = selected_mito)
		cleanRecord[isa, "keep_mito"] <- length(selected_mito)

		selected_ribo <- WhichCells(cleanObj, expression = percent_ribo > 5)
		cleanObj <- subset(cleanObj, cells = selected_ribo)
		cleanRecord[isa, "keep_ribo"] <- length(selected_ribo)

		C <- cleanObj@assays$RNA@counts
		C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
		most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
#		print(head(most_expressed))
		me_c <- as.data.frame(as.matrix(C)[most_expressed,])
		cells <- colnames(me_c)
		me_c$gene <- rownames(me_c)
		me_c <- gather(me_c, "cell", "c", cells)
#		print(head(me_c))
		top_gg <- ggplot(me_c, aes(x = gene, y = c)) +
			geom_boxplot(aes(fill = gene)) +
			coord_flip() +
			theme_classic()
		ggsave(paste(tmp_pf, "enriched_gene.png", sep = ""), dpi = 300, width = 9, height = 6)

		cleanObj <- cleanObj[!grepl("MALAT1", rownames(cleanObj)),]
		cleanObj <- cleanObj[!grepl("^MT-", rownames(cleanObj)),]
		cleanObj <- cleanObj[!grepl("^RP[SL]", rownames(cleanObj)),]
		cleanRecord[isa, "keep_feature"] <- nrow(cleanObj)

		if (nrow(cleanObj@meta.data) > 50) {
			cleanObj <- NormalizeData(cleanObj)
			cleanObj <- CellCycleScoring(object = cleanObj, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
	#		print(head(cleanObj@meta.data))
			cc_vln <- VlnPlot(cleanObj, features = c("S.Score", "G2M.Score"), pt.size = 0.1, ncol = 2) + NoLegend()
			ggsave(paste(tmp_pf, "cell_cycle_vln.png", sep = ""), cc_vln, dpi = 300, height = 3, width = 3)
	#		print(cleanObj)
			cleanObj <- FindVariableFeatures(cleanObj, verbose=T, selection.method = "vst", nfeatures = 2000)
			cleanObj <- ScaleData(cleanObj)
			cleanObj <- RunPCA(cleanObj)
			cleanObj <- RunUMAP(cleanObj, dims = 1:10)

			sweep_res_list <- paramSweep_v3(cleanObj, PCs = 1:10, sct = F)
			sweep_stats <- summarizeSweep(sweep_res_list, GT = F)
			bcmvn <- find.pK(sweep_stats)
	#		print(bcmvn) ## TODO: give the optimal pK!
			dblPct <- 0.0008*nrow(tmpObj@meta.data)+0.0527
			nExp_poi <- round(dblPct/100*nrow(cleanObj@meta.data))

			print(cleanObj)
			cleanObj <- doubletFinder_v3(cleanObj, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = F, sct = F)
			saveRDS(cleanObj, paste(tmp_pf, "cleaned_seurat_object_w_df.rds", sep = ""))
			print(cleanObj)

			adt_assay <- CreateAssayObject(counts = tmpData$`Antibody Capture`[,colnames(cleanObj)])
			cleanObj[['ADT']] <- adt_assay

			dfCol <- colnames(cleanObj@meta.data)[grepl("DF.classifications", colnames(cleanObj@meta.data))]
			dbl_c <- rownames(cleanObj@meta.data)[cleanObj@meta.data[,dfCol] == "Doublet"]
			sgl_c <- rownames(cleanObj@meta.data)[cleanObj@meta.data[,dfCol] == "Singlet"]

			dblObj <- subset(cleanObj, cells = dbl_c)
			saveRDS(dblObj, paste(tmp_pf, "cleaned_seurat_object_df_doublet.rds", sep = ""))
			cleanRecord[isa, "doublet"] <- ncol(dblObj)
			print(dblObj)

			sglObj <- subset(cleanObj, cells = sgl_c)
			saveRDS(sglObj, paste(tmp_pf, "cleaned_seurat_object_df_singlet.rds", sep = ""))
			cleanRecord[isa, "singlet"] <- ncol(sglObj)
			print(sglObj) # NOTE: doubletFinder will remove other assays!!!!
		} else {
			saveRDS(cleanObj, paste(tmp_pf, "cleaned_seurat_object_too_few_cells.rds", sep = ""))
		}
	
	}
	write.csv(cleanRecord, paste(expDir, "/", expID, "_clean_record.csv", sep = ""))
}


################################

####2:Data integration----------

################################

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

inputDir <- "6CloupeFile/results/"
inputID <- c("lc_chemo_icb_pbmc_221202")

sigMarker <- c("PTPRC", "EPCAM", "COL1A1", "CD3D", "CD8A", "CD4", "PDCD1", "TCF7", "NKG7", "MS4A1", "CD27", "CD19", 
               "CD68", "FCER1A", "HLA-DRA", "NCAM1", "TPSAB1", "GNLY", "FOXP3", "TBX21", "GATA3")

expID <- "SINGLET"
filePat <- "singlet.rds"

inteFlag <- TRUE
rnaInteMethod <- "RPCA"
adtInteMethod <- "RPCA"

dir.create(file.path(inputDir, expID), showWarnings = FALSE)
expDir <- paste(inputDir, expID, sep = "")
ppf <- paste(expDir, "/", expID, "_", sep = "")

if (inteFlag) {
  inFiles <- c()
  for (iii in inputID) {
    cat(iii, "\n")
    tmpDirs <- list.dirs(paste(inputDir, iii, sep = "/"), recursive = F, full.names = F)
    useDirs <- tmpDirs
    useFiles <- c()
    if (length(useDirs) > 0) {
      for (iud in useDirs) {
        tmpFiles <- list.files(paste(inputDir, iii, iud, sep = "/"), pattern = filePat)
        useFiles <- c(useFiles, paste(iud, tmpFiles, sep = "/"))
      }
    }
    useFiles <- useFiles[str_detect(useFiles, ".rds")]
    inFiles <- c(inFiles, str_c(iii, "/", useFiles))
  }
  sampleDf <- as.data.frame(str_split_fixed(inFiles, "\\/", n = 3))
  colnames(sampleDf) <- c("dataset_id", "sample_name", "input_file")
  sampleDf$dir <- inFiles
  sampleDf <- sampleDf[sampleDf$input_file != "",]
  print(sampleDf)
  write.csv(sampleDf, paste(ppf, "sample_table.csv", sep = ""))
  
  saName <- sampleDf$sample_name
  
  srCleanRnaObjs <- vector(mode = "list", length = length(saName))
  names(srCleanRnaObjs) <- saName
  srCleanAdtObjs <- vector(mode = "list", length = length(saName))
  names(srCleanAdtObjs) <- saName
  
  for (isa in 1:length(saName)) {
    cat(isa, "\n")
    tmp_rds <- paste(inputDir, sampleDf[isa, "dir"], sep = "")##sep= "/"改
    tmpObj <- readRDS(tmp_rds)
    print(tmpObj)
    tmpObj@meta.data$dataset_id <- sampleDf[isa, "dataset_id"]
    tmpObj@meta.data$sample_name <- sampleDf[isa, "sample_name"]
    tmpObj@meta.data$DF_classification <- tmpObj@meta.data[,str_detect(colnames(tmpObj@meta.data), "DF.classification")]
    tmpObj@meta.data$DF_pANN_score <- tmpObj@meta.data[,str_detect(colnames(tmpObj@meta.data), "pANN_")]
    
    if (ncol(tmpObj) > 50) {
      clearRnaObj <- CreateSeuratObject(counts = tmpObj@assays$RNA@counts, project = saName[isa], meta.data = tmpObj@meta.data)
      DefaultAssay(clearRnaObj) <- "RNA"
      if (rnaInteMethod == "standard" | rnaInteMethod == "RPCA") {
        clearRnaObj <- NormalizeData(clearRnaObj)
        clearRnaObj <- FindVariableFeatures(clearRnaObj)
      }
      srCleanRnaObjs[saName[isa]] <- clearRnaObj
      
      clearAdtObj <- CreateSeuratObject(counts = tmpObj@assays$ADT@counts, project = saName[isa], meta.data = tmpObj@meta.data)
      DefaultAssay(clearAdtObj) <- "RNA"
      if (adtInteMethod == "standard" | adtInteMethod == "RPCA") {
        clearAdtObj <- NormalizeData(clearAdtObj, normalization.method = "CLR", margin = 2)
        clearAdtObj <- FindVariableFeatures(clearAdtObj)
      }
      srCleanAdtObjs[saName[isa]] <- clearAdtObj
      
    } else {
      srCleanRnaObjs[saName[isa]] <- NULL
      srCleanAdtObjs[saName[isa]] <- NULL
    }
  }
  print(srCleanRnaObjs)
  print(srCleanAdtObjs)
  
  if (rnaInteMethod == "RPCA") {
    cat("Reciprocal PCA integration\n")
    features <- SelectIntegrationFeatures(object.list = srCleanRnaObjs)
    srCleanRnaObjs <- lapply(X = srCleanRnaObjs, FUN = function(x) {
      x <- ScaleData(x, features = features, verbose = TRUE)
      x <- RunPCA(x, features = features, verbose = TRUE)
    })
    integrateAnchors <- FindIntegrationAnchors(object.list = srCleanRnaObjs, reference = c(1,2), reduction = "rpca", dims = 1:50)
    inteRnaObjs <- IntegrateData(anchorset = integrateAnchors, k.weight = 50, dims = 1:50)
    DefaultAssay(inteRnaObjs) = "integrated"
  }
  saveRDS(inteRnaObjs, paste(ppf, rnaInteMethod, "_rna_integrated.rds", sep = ""))
  
  if (adtInteMethod == "RPCA") {
    cat("Reciprocal PCA integration\n")
    features <- SelectIntegrationFeatures(object.list = srCleanAdtObjs)
    srCleanAdtObjs <- lapply(X = srCleanAdtObjs, FUN = function(x) {
      x <- ScaleData(x, features = features, verbose = TRUE)
      x <- RunPCA(x, features = features, verbose = TRUE)
    })
    integrateAnchors <- FindIntegrationAnchors(object.list = srCleanAdtObjs, reference = c(1,2), reduction = "rpca")
    inteAdtObjs <- IntegrateData(anchorset = integrateAnchors, k.weight = 50)
    DefaultAssay(inteAdtObjs) = "integrated"
  }
  
  saveRDS(inteAdtObjs, paste(ppf, adtInteMethod, "_adt_integrated.rds", sep = ""))
  
  inteObjs <- inteRnaObjs
  inteObjs[["integrated_adt"]] <- inteAdtObjs[['integrated']]
  saveRDS(inteObjs, paste(ppf, "adt_", adtInteMethod,  "_rna_", rnaInteMethod, "_integrated.rds", sep = ""))
  writeLines(capture.output(sessionInfo()), paste(ppf, "integration_sessionInfo.txt", sep = ""))
}

