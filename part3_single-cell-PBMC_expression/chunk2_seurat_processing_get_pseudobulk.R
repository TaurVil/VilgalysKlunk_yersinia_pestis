# cd ./pestis_single_cell/
# module load R/4.0.4; module load hdf5; module load hdf5_hl; R

##### Setup dataset #########
library(Seurat); library(SeuratDisk)
library(ggplot2); library(patchwork)
library(data.table); library(dplyr)

##-----------------------------------------------
## get data for each sample
##-----------------------------------------------
for (i in c(names)) {
  # the original processing also included using the souporcell output to de-multiplex cells. genotypes were clustered using souporcell and the clusters were assigned to individuals by comparing genotype data with that previuosly generated (Randolph et al. 2021, science). We only retained cells determined to be singlets. 
  
  # start with null condition
  # read in data and transform into seurat object
  NI <- Read10X(data.dir =paste("./cellranger/",i,"_NI_OUT/outs/filtered_feature_bc_matrix/", sep=""))
  sNI <-CreateSeuratObject(NI, raw.data = NI, project = paste("NI_",i,sep=""))
  # assign variable for stimulation and percent mt reads
  sNI$stim <- "UNSTIM"; sNI[["percent.mt"]] <- PercentageFeatureSet(sNI, pattern = "^MT-")
  # filter out low quality cells
  sNI <- subset(sNI, subset = nFeature_RNA > 300 & percent.mt < 20)
  
  # repeat for stimulated samples
  YP <- Read10X(data.dir =paste("./cellranger/",i,"_YP_OUT/outs/filtered_feature_bc_matrix/", sep=""))
  sYP <-CreateSeuratObject(YP, raw.data = YP, project = paste("YP_",i,sep=""))
  sYP$stim <- "STIM"; sYP[["percent.mt"]] <- PercentageFeatureSet(sYP, pattern = "^MT-")
  sYP <- subset(sYP, subset = nFeature_RNA > 300 & percent.mt < 20)
  
  # merge together data for the individual, then merge with any other completed individuals
  merged <- merge(sNI, y = c(sYP), add.cell.ids = c("unstim", "stim"), project = "merged")
  if (i > 1) {all_data <- merge(all_data, y=merged)} else {merged -> all_data}
  rm(merged, sYP, YP_soup, YP, NI, sNI, NI_soup)
}; rm(i)
my_data <- all_data; rm(all_data)

## Make sure dataset contains the 10 individuals we expect, and assign ERAP2 genotypes to them (based on previous sequencing)
unique(paste(my_data@meta.data$batch, my_data@meta.data$sample,my_data@meta.data$stim, sep="_"))
{
  my_data$genotype <- NA
  my_data$genotype[my_data$sample_name == "HMN171216"] <- 0
  my_data$genotype[my_data$sample_name == "HMN171217"] <- 2
  my_data$genotype[my_data$sample_name == "HMN171221"] <- 2
  my_data$genotype[my_data$sample_name == "HMN171222"] <- 0
  my_data$genotype[my_data$sample_name == "HMN171225"] <- 0
  my_data$genotype[my_data$sample_name == "HMN171226"] <- 2
  my_data$genotype[my_data$sample_name == "HMN171227"] <- 2
  my_data$genotype[my_data$sample_name == "HMN171228"] <- 2
  my_data$genotype[my_data$sample_name == "HMN52544"] <- 2
  my_data$genotype[my_data$sample_name == "HMN52551"] <- 2
  my_data$genotype[my_data$sample_name == "HMN83553"] <- 2
  my_data$genotype[my_data$sample_name == "HMN83565"] <- 0
  my_data$genotype[my_data$sample_name == "HMN83576"] <- 0
  my_data$genotype[my_data$sample_name == "HMN83577"] <- 0
  
}
my_data <- subset(x=my_data, subset = sample_name %in% c("HMN171216", "HMN171217", "HMN171221", "HMN171222", 
                                                         "HMN171225", "HMN171226", "HMN171227", "HMN171228", 
                                                         "HMN52544", "HMN52551", "HMN83553", "HMN83565", "HMN83576", "HMN83577"))	

## save intermediate output file that we can return too later. This is especially important if unsure about the total memory available and you're worried about crashing your session if the R objects become too large. 
saveRDS(my_data, file = "./pbmc_n10_EU_hg38pestis.rds")
##-----------------------------------------------

##-----------------------------------------------
## Predict cell types based on the reference PBMCs (following: https://satijalab.org/seurat/articles/multimodal_reference_mapping.html)
##-----------------------------------------------
readRDS("./pbmc_n10_EU_hg38pestis.rds") -> my_data
reference <- LoadH5Seurat("./pbmc_multimodal.h5seurat")
{
  my_data <- SCTransform(my_data, verbose = FALSE)
  DefaultAssay(my_data) <- 'SCT'
  anchors <- FindTransferAnchors(
    reference = reference,
    query = my_data,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50
  )
  
  ## Error message here saying it's using RNA is only for the UMAP
  my_data <- MapQuery(
    anchorset = anchors,
    query = my_data,
    reference = reference, 
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    ),
    reference.reduction = "spca", 
    reduction.model = "wnn.umap"
  )
  
}; rm(reference, anchors)
##-----------------------------------------------

##-----------------------------------------------
## Integrate infected to null samples
##-----------------------------------------------
{
  DefaultAssay(my_data) <- 'SCT'
  
  #split into separate seurat objects 
  pbmc_list <- SplitObject(my_data, split.by = "stim")
  
  ## Follow https://satijalab.org/seurat/archive/v3.0/integration.html
  # SCT transform for each part of the dataset (null/infected)
  for (i in 1:length(pbmc_list)) {
    pbmc_list[[i]] <- SCTransform(pbmc_list[[i]], verbose = FALSE)
  }; rm(i)
  # select features for downstream integration, and run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated.
  pbmc.features <- SelectIntegrationFeatures(object.list = pbmc_list, nfeatures = 5000)
  options(future.globals.maxSize= 2*891289600)
  pbmc_list <- PrepSCTIntegration(object.list = pbmc_list, anchor.features = pbmc.features, 
                                      verbose = FALSE)
  #  identify anchors and integrate the datasets
  pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc_list, normalization.method = "SCT", 
                                         anchor.features = pbmc.features, verbose = FALSE)
  pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, normalization.method = "SCT", 
                                       verbose = FALSE)
  
  # # standard preprocessing (log normalization), identify variable features independently 
  # for (i in 1:length(pbmc_list)) {
  #   pbmc_list[[i]] <- NormalizeData(pbmc_list[[i]], verbose = FALSE)
  #   pbmc_list[[i]] <- FindVariableFeatures(pbmc_list[[i]], selection.method = "vst", 
  #                                          nfeatures = 2000, verbose = FALSE)
  #   print(i)
  # }

  # switch to integrated dataset (original still stored in SCT/RNA slot)
  DefaultAssay(pbmc.integrated) <- "integrated"
  
  # cluster and visualize
  pbmc_integrated <- RunPCA(pbmc.integrated, npcs = 30, verbose = FALSE)
  
  # tSNE and clustering
  pbmc_integrated <- RunUMAP(pbmc_integrated, reduction = "pca", dims = 1:30)
  
  pbmc_integrated <- FindNeighbors(pbmc_integrated, reduction = "pca", dims = 1:20)
  pbmc_integrated <- FindClusters(pbmc_integrated, resolution = 0.5)
  
}
saveRDS(pbmc_integrated, file = "./pbmc_integrated_n10_EU_hg38pestis.rds")
##-----------------------------------------------

##-----------------------------------------------
## Reintegrate by individual 
##-----------------------------------------------
{
  DefaultAssay(my_data) <- 'SCT'
  my_data@meta.data$gt_stim <- paste(my_data@meta.data$genotype, my_data@meta.data$stim, sep="_")
  
  #split into separate seurat objects 
  pbmc_list <- SplitObject(my_data, split.by = "gt_stim"); rm(my_data); gc()
  
  ## Follow https://satijalab.org/seurat/articles/integration_rpca.html for RPCA
  # SCT transform for each part of the dataset (null/infected)
  pbmc_list <- lapply(X = pbmc_list, FUN = SCTransform) # , method = "glmGamPoi", faster
  pbmc.features <- SelectIntegrationFeatures(object.list = pbmc_list, nfeatures = 3000)
  pbmc_list <- PrepSCTIntegration(object.list = pbmc_list, anchor.features = pbmc.features)
  pbmc_list <- lapply(X = pbmc_list, FUN = RunPCA, features = pbmc.features)
  
  pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc_list, normalization.method = "SCT",
                                         anchor.features = pbmc.features, dims = 1:30, 
                                         reduction = "rpca", k.anchor = 20)
  rm(pbmc_list); gc()
  pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, normalization.method = "SCT", dims = 1:30)
  
  # switch to integrated dataset (original still stored in SCT/RNA slot)
  DefaultAssay(pbmc.integrated) <- "integrated"
  
  # cluster and visualize
  pbmc_integrated <- RunPCA(pbmc.integrated, npcs = 30, verbose = FALSE)
  
  # tSNE and clustering
  pbmc_integrated <- RunUMAP(pbmc_integrated, reduction = "pca", dims = 1:30)
  
  pbmc_integrated <- FindNeighbors(pbmc_integrated, reduction = "pca", dims = 1:20)
  pbmc_integrated <- FindClusters(pbmc_integrated, resolution = 0.5)
  
}
saveRDS(pbmc_integrated, file = "./pbmc_integrated_by_genotype-condition_n10_EU_hg38pestis.rds")
{
  DefaultAssay(my_data) <- 'SCT'
  #split into separate seurat objects 
  pbmc_list <- SplitObject(my_data, split.by = "id_batch_sample"); rm(my_data); gc()
  
  ## Follow https://satijalab.org/seurat/articles/integration_rpca.html for RPCA
  # SCT transform for each part of the dataset (null/infected)
  pbmc_list <- lapply(X = pbmc_list, FUN = SCTransform) # , method = "glmGamPoi", faster
  pbmc.features <- SelectIntegrationFeatures(object.list = pbmc_list, nfeatures = 3000)
  pbmc_list <- PrepSCTIntegration(object.list = pbmc_list, anchor.features = pbmc.features)
  pbmc_list <- lapply(X = pbmc_list, FUN = RunPCA, features = pbmc.features)
  
  pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc_list, normalization.method = "SCT",
                                         anchor.features = pbmc.features, dims = 1:30, 
                                         reduction = "rpca", k.anchor = 20)
  rm(pbmc_list); gc()
  pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, normalization.method = "SCT", dims = 1:30)
  
  # switch to integrated dataset (original still stored in SCT/RNA slot)
  DefaultAssay(pbmc.integrated) <- "integrated"
  
  # cluster and visualize
  pbmc_integrated <- RunPCA(pbmc.integrated, npcs = 30, verbose = FALSE)
  
  # tSNE and clustering
  pbmc_integrated <- RunUMAP(pbmc_integrated, reduction = "pca", dims = 1:30)
  
  pbmc_integrated <- FindNeighbors(pbmc_integrated, reduction = "pca", dims = 1:20)
  pbmc_integrated <- FindClusters(pbmc_integrated, resolution = 0.5)
  
}
saveRDS(pbmc_integrated, file = "./pbmc_integrated_by_genotype-condition_n10_EU_hg38pestis.rds")
##-----------------------------------------------

##### 
my_data <- readRDS("./pbmc_integrated_n10_EU_hg38pestis.rds")
cell_types <- c("NK", "CD8 T", "CD4 T", "B", "Mono")

##-----------------------------------------------
## Calculate all PBMC pseudobulk means
##-----------------------------------------------
pseudobulk_means <- NULL
for (x in 1:length(cell_types)) {
  celltype_i <- cell_types[x]
  ### get cell type of interest ###
  seurat_obj <- subset(my_data, predicted.celltype.l1 == cell_types[x])
  ### convert to SCE class ###
  sce <- as.SingleCellExperiment(seurat_obj, assay="RNA")
  ## QC metrics ##
  sce <- addPerCellQC(sce)
  sce <- addPerFeatureQC(sce)
  
  
  ##############################################
  ## compute logcounts values (log2(CPM + 1)) ##
  ## filter lowly-expressed genes ##
  # keep genes with median logCPM > 1
  logcounts(sce) <- log2(calculateCPM(sce, size_factors = NULL) + 1)
  keep_gene <- rowMeans(logcounts(sce)) > logCPM_filter
  sce_filt <- sce[keep_gene,]
  
  ## record number of genes
  if(x == 1){
    genes <- table(keep_gene)[2]
    number_genes <- c(genes)
    cell_type <- c(celltype_i)
  }else{
    genes <- table(keep_gene)[2]
    number_genes <- c(number_genes, genes)
    cell_type <- c(cell_type, celltype_i)
  }
  num_genes_out <- as.data.frame(cbind(cell_type, number_genes))
  rownames(num_genes_out) <- num_genes_out$cell_type
  
  ######################################
  ## normalization with scran factors ##
  # compute sum factors
  sce_filt <- computeSumFactors(sce_filt)
  sce_filt$size_factor <- sizeFactors(sce_filt)
  
  colData(sce_filt)$batch <- as.factor(colData(sce_filt)$batch)
  
  # pdf(paste0(filtered_log2CPM_dir,celltype_i,"_distribution_sizeFactors_by_batch.pdf"))
  # colData(sce_filt) %>% as.data.frame %>%
  #   ggplot(aes(x = log2(size_factor), y = batch, fill = stim)) +
  #   geom_density_ridges(alpha = 0.5) + geom_vline(xintercept = 0, linetype = 2) +
  #   theme_ridges() + ggtitle(paste0(celltype_i, ", distribution of size factors")) +
  #   scale_fill_tableau()
  # dev.off()
  
  sce_filt <- sce_filt[, (sizeFactors(sce_filt) < 8 & sizeFactors(sce_filt) > 0.125)]
  sce_filt <- logNormCounts(sce_filt)
  
  ######################################
  ## compute averages by individual and condition ##
  normalized_data <- exprs(sce_filt)
  meta_data <- as.data.frame(colData(sce_filt))
  
  sample_colname <- "id_batch_sample"
  IDs <- as.data.frame(meta_data)[, sample_colname]
  unique_ID_list <- as.list(unique(IDs))
  pseudobulk <- as.data.frame(pblapply(unique_ID_list, FUN = function(x){rowMeans(normalized_data[,IDs == x, drop = FALSE])}))
  cellcount <- as.data.frame(pblapply(unique_ID_list, FUN = function(x){ncol(normalized_data[,IDs == x, drop = FALSE])}))
  colnames(pseudobulk) <- names(cellcount) <- unique_ID_list
  rownames(pseudobulk) <- rownames(exprs(sce_filt))
  
  ######################################
  ## Output files ##
  pseudobulk_means[[celltype_i]] <- pseudobulk
  if (x == 1) {
    as.data.frame(t(cellcount)) -> all_counts; colnames(all_counts) <- celltype_i
  } else {
    cbind(all_counts, t(cellcount)) -> all_counts; colnames(all_counts)[x] <- celltype_i
  }
  
  print(celltype_i)
  rm(celltype_i, pseudobulk, genes, keep_gene, cellcount, normalized_data, meta_data, unique_ID_list, sce_filt)
}; rm(x)
rm(IDs, cell_type, dec, number_genes, seurat_obj, sce)
rm(my_data)
save.image("./pbmc_n10_pseudobulk_means.RData")
