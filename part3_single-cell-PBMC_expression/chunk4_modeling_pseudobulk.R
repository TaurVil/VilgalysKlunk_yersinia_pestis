library(ggplot2); library(patchwork)
library(data.table); library(dplyr)
library(limma); library(edgeR); library(plyr)

##-----------------------------------------------
## Start with PBMC pseudobulk means
##-----------------------------------------------
load("./DATA/pbmc_n10_pseudobulk_means.RData")
##-----------------------------------------------

##-----------------------------------------------
## get metadata for each sample
##-----------------------------------------------
## make a matrix of sample names
names <- strsplit(colnames(pseudobulk_means$B), ":"); ldply(names) -> names; 
{
  names$sample_name[names$V1 == "1_2_UNSTIM"] <- "HMN171216"
  names$sample_name[names$V1 == "1_2_STIM"] <- "HMN171216"
  names$sample_name[names$V1 == "1_0_UNSTIM"] <- "HMN171217"
  names$sample_name[names$V1 == "1_4_STIM"] <- "HMN171217"
  names$sample_name[names$V1 == "3_2_UNSTIM"] <- "HMN171221"
  names$sample_name[names$V1 == "3_5_STIM"] <- "HMN171221"
  names$sample_name[names$V1 == "3_1_UNSTIM"] <- "HMN171222"
  names$sample_name[names$V1 == "3_1_STIM"] <- "HMN171222"
  names$sample_name[names$V1 == "4_1_UNSTIM"] <- "HMN171225"
  names$sample_name[names$V1 == "4_3_STIM"] <- "HMN171225"
  names$sample_name[names$V1 == "4_2_UNSTIM"] <- "HMN171226"
  names$sample_name[names$V1 == "4_4_STIM"] <- "HMN171226"
  names$sample_name[names$V1 == "5_5_UNSTIM"] <- "HMN171227"
  names$sample_name[names$V1 == "5_1_STIM"] <- "HMN171227"
  names$sample_name[names$V1 == "5_0_UNSTIM"] <- "HMN171228"
  names$sample_name[names$V1 == "5_2_STIM"] <- "HMN171228"
  
  names$sample_name[names$V1 == "13_0_UNSTIM"] <- "HMN83565"
  names$sample_name[names$V1 == "13_5_STIM"] <- "HMN83565"
  names$sample_name[names$V1 == "14_3_UNSTIM"] <- "HMN83576"
  names$sample_name[names$V1 == "14_2_STIM"] <- "HMN83576"
  
}
## add whether a sample is stimulated
names$stim <- c(rep("0_unstim", 10), rep("1_stim", 10))
## add in target genotype
names$geno <- NA; {
  names$geno[names$sample_name == "HMN171216"] <- 0
  names$geno[names$sample_name == "HMN171217"] <- 2
  names$geno[names$sample_name == "HMN171221"] <- 2
  names$geno[names$sample_name == "HMN171222"] <- 0
  names$geno[names$sample_name == "HMN171225"] <- 0
  names$geno[names$sample_name == "HMN171226"] <- 2
  names$geno[names$sample_name == "HMN171227"] <- 2
  names$geno[names$sample_name == "HMN171228"] <- 2
  names$geno[names$sample_name == "HMN83565"] <- 0
  names$geno[names$sample_name == "HMN83576"] <- 0
}
info <- names; rm(names)
##-----------------------------------------------

##-----------------------------------------------
# plot expression per sample (no outputs)
# genes with 0 counts (after normalization) form peak at -3
##-----------------------------------------------
## Plot B cell distribution of expression per sample
par(mfrow=c(2,5)); for (i in 1:ncol(pseudobulk_means$B)) {
  tmp <- pseudobulk_means$B[,i]; tmp[tmp == 0] <- 1e-3
  hist(log10(tmp), breaks=50, xlab = "log10(pseudomeans)",main=paste(info$sample_name[i], info$stim[i], "\ncells:", all_counts[i,4], sep=" "))
}

## Plot CD4 T cell distribution of expression per sample
par(mfrow=c(2,5)); for (i in 1:ncol(pseudobulk_means$`CD4 T`)) {
  tmp <- pseudobulk_means$`CD4 T`[,i]; tmp[tmp == 0] <- 1e-3
  hist(log10(tmp), breaks=50, xlab = "log10(pseudomeans)", main=paste(info$sample_name[i], info$stim[i], "\ncells:", all_counts[i,3], sep=" "))
}

## Plot Monocyte distribution of expression per sample
par(mfrow=c(2,5)); for (i in 1:ncol(pseudobulk_means$`Mono`)) {
  tmp <- pseudobulk_means$`Mono`[,i]; tmp[tmp == 0] <- 1e-3
  hist(log10(tmp), breaks=50, xlab = "log10(pseudomeans)", main=paste(info$sample_name[i], info$stim[i], "\ncells:", all_counts[i,5], sep=" "))
}

## Plot CD8 T cell distribution of expression per sample
par(mfrow=c(2,5)); for (i in 1:ncol(pseudobulk_means$`CD8 T`)) {
  tmp <- pseudobulk_means$`CD8 T`[,i]; tmp[tmp == 0] <- 1e-3
  hist(log10(tmp), breaks=50, xlab = "log10(pseudomeans)", main=paste(info$sample_name[i], info$stim[i], "\ncells:", all_counts[i,2], sep=" "))
}
rm(i, tmp)
##-----------------------------------------------

##-----------------------------------------------
## model expression for each cell type
## voom normalization, stimulation + individual + cells
##-----------------------------------------------
voom_data <- NULL
limma_stim <- NULL ; for (i in 1:length(cell_types)) { 
  data <- pseudobulk_means[[i]]
  names <- info[order(match(info$V1, colnames(data))),]
  sum(names$V1 == colnames(data))/ncol(data)
  
  # calculate cpm and filter lowly expressed genes 
  lib.size <- colSums(data)/1e6
  cpm <- t(t(data)/lib.size)
  rm(lib.size)
  data <- data[rowMeans(cpm) >= min(c(70, mean(cpm[row.names(cpm) == "ERAP2",]), mean(cpm[row.names(cpm) == "ERAP1",]))),]
  mean_cpm <- rowMeans(cpm[rowMeans(cpm) >= min(c(70, mean(cpm[row.names(cpm) == "ERAP2",]), mean(cpm[row.names(cpm) == "ERAP1",]))),])
  
  
  # add number of cells to our metadata 
  names$cells <- NA
  for (j in 1:nrow(names)) {
    names$cells[j] <- all_counts[row.names(all_counts) == names$V1[j], i]
  }; rm(j)
  
  # create design matrix
  stim <- names$stim; cells <- names$cells; samples <- as.factor(names$sample_name)
  design <- model.matrix(~ 1 + stim + cells + samples)  
  colnames(design)[1:2] <- c("intercept", "stim")
  rm(stim, cells, samples)
  # + pop:stimulated , data = exp
  
  # voom normalize, then model using lmfit 
  v <- voom(data, design, plot=T)
  fit = lmFit(v)
  voom_data[[cell_types[i]]] <- v$E
  
  contr.matrix <- makeContrasts(
    Stimulation = stim,
    levels=colnames(design)
  )
  fit2 <- contrasts.fit(fit = fit, contrasts = contr.matrix)
  fit2 = eBayes(fit2)
  
  # extract summary table + cell type + gene names
  res_stim <- topTreat(fit2, coef=1, n=Inf)
  res_stim$test <- "stimulation"
  res_stim$celltype <- cell_types[i]
  res_stim$gene <- rownames(res_stim)
  
  SE <- sqrt(fit2$s2.post) * fit2$stdev.unscaled; colnames(SE) <- paste("SE_", colnames(SE), sep="")
  B22 <- fit2$coefficients; colnames(B22) <- paste("B_", colnames(B22), sep="")
  #fit2$p.value[row.names(data) == "ERAP2",]
  #fit2$p.value[row.names(data) == "ERAP1",]
  output <- as.data.frame(cbind(fit2$p.value, B22, SE)); rm(SE, B22)
  output$gene <- rownames(output)
  
  output <- merge(output, res_stim, by="gene")
  
  output$cpm <- mean_cpm
  limma_stim <- rbind(limma_stim, output); rm(output)
  rm(fit2, fit, v, data, cpm, mean_cpm, res_stim, contr.matrix, names, design)
  print(cell_types[i])
  
}; rm(i)

ggplot(data=limma_stim, aes(y=B_Stimulation , x=celltype)) + geom_boxplot(fill='gray') + geom_hline(yintercept = 0, col='black') + theme_classic()
ggplot(data=limma_stim, aes(y=-log10(P.Value), x=celltype)) + geom_violin(fill='gray') + geom_hline(yintercept = 0, col='black') + theme_classic()
##-----------------------------------------------

##-----------------------------------------------
## model expression for each cell type
## voom normalization, stimulation + individual + cells
##-----------------------------------------------

limma_geno <- NULL ; for (i in 1:length(cell_types)) { 
  data <- pseudobulk_means[[i]]
  names <- info[order(match(info$V1, colnames(data))),]
  sum(names$V1 == colnames(data))/ncol(data)
  
  # calculate cpm and filter lowly expressed genes 
  lib.size <- colSums(data)/1e6
  cpm <- t(t(data)/lib.size)
  rm(lib.size)
  data <- data[rowMeans(cpm) >= min(c(70, mean(cpm[row.names(cpm) == "ERAP2",]), mean(cpm[row.names(cpm) == "ERAP1",]))),]
  mean_cpm <- rowMeans(cpm[rowMeans(cpm) >= min(c(70, mean(cpm[row.names(cpm) == "ERAP2",]), mean(cpm[row.names(cpm) == "ERAP1",]))),])
  
  
  # add number of cells to our metadata 
  names$cells <- NA
  for (j in 1:nrow(names)) {
    names$cells[j] <- all_counts[row.names(all_counts) == names$V1[j], i]
  }; rm(j)
  
  # create design matrix
  stim <- names$stim; cells <- names$cells; samples <- as.factor(names$sample_name)
  gt <- names$geno
  design <- model.matrix(~ 1 + stim + gt + cells)  
  colnames(design)[1:3] <- c("intercept", "stim", "genotype")
  rm(stim, cells, samples, gt)
  # + pop:stimulated , data = exp
  
  # voom normalize, then model using lmfit 
  v <- voom(data, design, plot=T)
  fit = lmFit(v)
  
  contr.matrix <- makeContrasts(
    genotype = genotype,
    levels=colnames(design)
  )
  fit2 <- contrasts.fit(fit = fit, contrasts = contr.matrix)
  fit2 = eBayes(fit2)
  
  # extract summary table + cell type + gene names
  res_stim <- topTreat(fit2, coef=1, n=Inf)
  res_stim$test <- "genotype"
  res_stim$celltype <- cell_types[i]
  res_stim$gene <- rownames(res_stim)
  
  SE <- sqrt(fit2$s2.post) * fit2$stdev.unscaled; colnames(SE) <- paste("SE_", colnames(SE), sep="")
  B22 <- fit2$coefficients; colnames(B22) <- paste("B_", colnames(B22), sep="")
  #fit2$p.value[row.names(data) == "ERAP2",]
  #fit2$p.value[row.names(data) == "ERAP1",]
  output <- as.data.frame(cbind(fit2$p.value, B22, SE)); rm(SE, B22)
  output$gene <- rownames(output)
  
  output <- merge(output, res_stim, by="gene")
  
  output$cpm <- mean_cpm
  limma_geno <- rbind(limma_geno, output); rm(output)
  rm(fit2, fit, v, data, cpm, mean_cpm, res_stim, contr.matrix, names, design)
  print(cell_types[i])
  
}; rm(i)

ggplot(data=limma_geno, aes(y=B_genotype , x=celltype)) + geom_boxplot(fill='gray') + geom_hline(yintercept = 0, col='black') + theme_classic()
ggplot(data=limma_geno, aes(y=-log10(P.Value), x=celltype)) + geom_violin(fill='gray') + geom_hline(yintercept = 0, col='black') + theme_classic()
##-----------------------------------------------

##-----------------------------------------------
# Plot effect by expression quantile
##-----------------------------------------------
i=5
tmp_data <- limma_geno[limma_geno$celltype == unique(limma_geno$celltype)[i],]
new_tmp <- tmp_data %>%
  mutate(quantile = ntile(cpm,10))
ggplot(new_tmp, aes(x=quantile, y=B_genotype, group=quantile)) + 
  geom_violin() + geom_boxplot(width=0.2) + 
  geom_hline(yintercept=0, color="grey") + 
  ggtitle(paste("quantiles ranked by cpm\ncell type:",unique(limma_geno$cell_type)[i],sep="")) + theme_classic()
print(cor.test(new_tmp$`B_geno`, new_tmp$cpm, method='spearman'))
rm(i, new_tmp, tmp_data)
##-----------------------------------------------

##-----------------------------------------------
# Plot target genes
##-----------------------------------------------

sig <- c(
  "CTLA4", "ICOS", # chr2
  "ERAP1", "ERAP2", "LNPEP", # chr5
  "TICAM2", "TMED7", "TMED7-TICAM2", # chr5
  "NFATC1" # chr18    
)

for (tmp_gene in sig) {
  print(tmp_gene)
  print(limma_geno[limma_stim$gene == tmp_gene,])
  par(mfrow=c(2,3))
  for (i in c(1:5)) {
    nom <- cell_types[i]
    tmp_data <- voom_data[[nom]]
    tmp_data <- tmp_data[row.names(tmp_data) == tmp_gene,]
    if (length(tmp_data) > 0) {boxplot(tmp_data ~ as.factor(info$stim), xlab="stimulation",
                                       ylab="voom-normalized expression", frame=F, 
                                       main=paste(nom, tmp_gene, sep=": "))
      points(tmp_data ~ jitter(as.numeric(as.factor(info$stim))), pch=16, col='purple4')
    }}
}
rm(i, tmp_gene, nom, tmp_data, sample_colname)
##-----------------------------------------------

##-----------------------------------------------
# forest plots
##-----------------------------------------------

tmp_data <- limma_stim[limma_stim$gene == "ERAP2",]

pl1 <- ggplot(data = tmp_data, aes(x=celltype, y=B_Stimulation , color = Stimulation  < 0.05)) + 
  geom_point(size=4) + 
  geom_errorbar(aes(ymin=B_Stimulation -SE_Stimulation, ymax=B_Stimulation +SE_Stimulation)) + 
  scale_colour_manual(name = 'PC1 > 0', values = setNames(c('black','gray'),c(T, F))) +
  coord_flip() + geom_hline(yintercept=0, linetype="dashed", 
                            color = "gray", size=1.2) + 
  ggtitle("ERAP2") + 
  scale_x_discrete(labels= rev(c("NK cells", "Monocytes", "CD8 T cells", "CD4 T Cells", "B Cells"))) + 
  theme_classic() + xlab(label = "cell type") + ylab(label="standardized effect size") + theme(legend.position = "none") 


tmp_data <- limma_stim[limma_stim$gene == "ERAP1",]
pl2 <- ggplot(data = tmp_data, aes(x=celltype, y=B_Stimulation , color = Stimulation  < 0.05)) + 
  geom_point(size=4) + 
  geom_errorbar(aes(ymin=B_Stimulation -SE_Stimulation, ymax=B_Stimulation +SE_Stimulation)) + 
  scale_colour_manual(name = 'PC1 > 0', values = setNames(c('black','gray'),c(T, F))) +
  coord_flip() + geom_hline(yintercept=0, linetype="dashed", 
                            color = "gray", size=1.2) + 
  ggtitle("ERAP1") + 
  scale_x_discrete(labels= rev(c("NK cells", "Monocytes", "CD8 T cells", "CD4 T Cells", "B Cells"))) + 
  theme_classic() + xlab(label = "cell type") + ylab(label="standardized effect size") + theme(legend.position = "none") 

tmp_data <- limma_stim[limma_stim$gene == "LNPEP",]
pl3 <- ggplot(data = tmp_data, aes(x=celltype, y=B_Stimulation , color = Stimulation  < 0.05)) + 
  geom_point(size=4) + 
  geom_errorbar(aes(ymin=B_Stimulation -SE_Stimulation, ymax=B_Stimulation +SE_Stimulation)) + 
  scale_colour_manual(name = 'PC1 > 0', values = setNames(c('black','gray'),c(T, F))) +
  coord_flip() + geom_hline(yintercept=0, linetype="dashed", 
                            color = "gray", size=1.2) + 
  ggtitle("LNPEP") + 
  scale_x_discrete(labels= rev(c("NK cells", "Monocytes", "CD8 T cells", "CD4 T Cells", "B Cells"))) + 
  theme_classic() + xlab(label = "cell type") + ylab(label="standardized effect size") + theme(legend.position = "none") 

pl1 + pl2 + pl3

tmp_data <- limma_stim[limma_stim$gene %in% c("ICOS","TMED7", "NFATC1"),]
ggplot(data = tmp_data, aes(x=paste(celltype, gene, sep="_"), y=B_Stimulation , color = Stimulation  < 0.05)) + 
  geom_point(size=4) + 
  geom_errorbar(aes(ymin=B_Stimulation -SE_Stimulation, ymax=B_Stimulation +SE_Stimulation)) + 
  scale_colour_manual(name = 'PC1 > 0', values = setNames(c('black','gray'),c(T, F))) +
  coord_flip() + geom_hline(yintercept=0, linetype="dashed", 
                            color = "gray", size=1.2) + 
  # scale_x_discrete(labels= rev(c("NK cells", "Monocytes", "CD8 T cells", "CD4 T Cells", "B Cells"))) + 
  theme_classic() + xlab(label = "cell type") + ylab(label="standardized effect size") + theme(legend.position = "none")  + 
  theme(text = element_text(size=14))


##-----------------------------------------------
# forest plots -- genotype
##-----------------------------------------------

tmp_data <- limma_geno[limma_geno$gene == "ERAP2",]

pl1 <- ggplot(data = tmp_data, aes(x=celltype, y=B_genotype , color = genotype  < 0.05)) + 
  geom_point(size=4) + 
  geom_errorbar(aes(ymin=B_genotype -SE_genotype, ymax=B_genotype +SE_genotype)) + 
  scale_colour_manual(name = 'PC1 > 0', values = setNames(c('black','gray'),c(T, F))) +
  coord_flip() + geom_hline(yintercept=0, linetype="dashed", 
                            color = "gray", size=1.2) + 
  ggtitle("ERAP2") + 
  scale_x_discrete(labels= rev(c("NK cells", "Monocytes", "CD8 T cells", "CD4 T Cells", "B Cells"))) + 
  theme_classic() + xlab(label = "cell type") + ylab(label="standardized effect size") + theme(legend.position = "none") 


tmp_data <- limma_geno[limma_geno$gene == "ERAP1",]
pl2 <- ggplot(data = tmp_data, aes(x=celltype, y=B_genotype , color = genotype  < 0.05)) + 
  geom_point(size=4) + 
  geom_errorbar(aes(ymin=B_genotype -SE_genotype, ymax=B_genotype +SE_genotype)) + 
  scale_colour_manual(name = 'PC1 > 0', values = setNames(c('black','gray'),c(T, F))) +
  coord_flip() + geom_hline(yintercept=0, linetype="dashed", 
                            color = "gray", size=1.2) + 
  ggtitle("ERAP1") + 
  scale_x_discrete(labels= rev(c("NK cells", "Monocytes", "CD8 T cells", "CD4 T Cells", "B Cells"))) + 
  theme_classic() + xlab(label = "cell type") + ylab(label="standardized effect size") + theme(legend.position = "none") 

tmp_data <- limma_geno[limma_geno$gene == "LNPEP",]
pl3 <- ggplot(data = tmp_data, aes(x=celltype, y=B_genotype , color = genotype  < 0.05)) + 
  geom_point(size=4) + 
  geom_errorbar(aes(ymin=B_genotype -SE_genotype, ymax=B_genotype +SE_genotype)) + 
  scale_colour_manual(name = 'PC1 > 0', values = setNames(c('black','gray'),c(T, F))) +
  coord_flip() + geom_hline(yintercept=0, linetype="dashed", 
                            color = "gray", size=1.2) + 
  ggtitle("LNPEP") + 
  scale_x_discrete(labels= rev(c("NK cells", "Monocytes", "CD8 T cells", "CD4 T Cells", "B Cells"))) + 
  theme_classic() + xlab(label = "cell type") + ylab(label="standardized effect size") + theme(legend.position = "none") 

pl1 + pl2 + pl3

# fig3e <- pl1 + pl3 + plot_layout(guides="collect") & theme(legend.position = 'bottom')
# save(fig3e, file="~/Barreiro/pestis_aDNA/fig3E.rdata")

##-----------------------------------------------
# 
##-----------------------------------------------

