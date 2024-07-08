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

set.seed(42)
plimma_stim <- NULL ; for (rep in 1:25) {for (i in 1:length(cell_types)) { 
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
  
  # permute stimulation within a sample 
  perms <- names
  for (nomnom in unique(perms$sample_name)) {
    perms$stim[perms$sample_name == nomnom] <- sample(perms$stim[perms$sample_name == nomnom])
  }
  
  # create design matrix
  stim <- perms$stim; cells <- names$cells; samples <- as.factor(names$sample_name)
  design <- model.matrix(~ 1 + stim + cells + samples)  
  colnames(design)[1:2] <- c("intercept", "stim")
  rm(stim, cells, samples, perms)
  # + pop:stimulated , data = exp
  
  # voom normalize, then model using lmfit 
  v <- voom(data, design, plot=F)
  fit = lmFit(v)
  
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
  plimma_stim <- rbind(plimma_stim, output); rm(output)
  rm(fit2, fit, v, data, cpm, mean_cpm, res_stim, contr.matrix, names, design)
  print(cell_types[i])
  
}; rm(i)}; rm(rep)
dim(plimma_stim)

library(cobs); perm.fdr=function(input_df,perm_df,Pvals_col_name,name){
  
  pvals_index=which(colnames(input_df)==Pvals_col_name)
  ro<-input_df[order(input_df[,pvals_index]),]
  p_obs <- data.frame(pvalue=ro[,pvals_index])
  p_vector<-matrix(as.matrix(perm_df),ncol=1)
  p_vector=data.frame(p_vector[order(p_vector)])
  
  F<-p_obs[,1]
  F_o<-p_obs[,1]
  pi_hat<-p_obs[,1]
  
  j=1
  observed<-length(p_obs[,1])
  randoms<-length(p_vector[,1])
  
  for(i in 1:observed)
  {
    repeat
    {
      if((p_vector[j,1]<p_obs[i,1])&j<randoms){j<-j+1}else{break}
    }
    F[i]=i/observed
    F_o[i]=(j-1)/randoms
    if(F_o[i]<1){pi_hat[i]=(1-F[i])/(1-F_o[i])}else{pi_hat[i]=1}
  }
  tabla <-data.frame(pi_hat,pval=p_obs[,1])
  
  tabla[1,]=c(1,0)
  last_percentile_average=mean(tabla$pi_hat[as.integer(min((length(tabla[,1])*0.99),(nrow(tabla)-1)):length(tabla[,1]))])
  tabla[nrow(tabla),]=c(last_percentile_average,1)
  constraint_matrix=as.matrix(data.frame(c(0,2),c(0,1),c(1,0)))
  f_hat<-suppressWarnings(cobs(tabla$pval,tabla$pi_hat,constraint="convex",pointwise=constraint_matrix,maxiter=1000,print.warn=FALSE,print.mesg=FALSE))
  
  f_hat_serie=f_hat$fitted
  pi_o=f_hat_serie[length(f_hat_serie)]
  pi_o=min(pi_o,1)
  pi_o=max(pi_o,0)
  
  Fdr_ST_perm=pi_o*F_o/F
  
  for(i in 1:length(p_obs[,1]))
  {
    Fdr_ST_perm[i]=pi_o*F_o[i]/F[i]
    if(i>1)
    {
      for(j in 1:(i-1))
      {
        if(Fdr_ST_perm[i-j]>Fdr_ST_perm[i]){Fdr_ST_perm[i-j]=Fdr_ST_perm[i]}else{break}
      }
    }
    if(Fdr_ST_perm[i]>1)  Fdr_ST_perm[i]=1
  }
  
  fdrs_df <-data.frame(ro,q_ST_perm=Fdr_ST_perm)
  rownames(fdrs_df)=rownames(ro)
  colnames(fdrs_df)[ncol(fdrs_df)]=paste0("fdr_",name)
  
  return(fdrs_df)
}

qvals <- NULL
for (i in cell_types) {
  tmp <- perm.fdr(limma_stim[limma_stim$celltype == i,], Pvals_col_name = "P.Value", plimma_stim$P.Value[plimma_stim$celltype == i], "permstim")
  rbind(qvals, tmp) -> qvals; rm(tmp)
}; rm(i)

res <- as.data.frame(matrix(ncol=1, nrow=length(cell_types))); colnames(res) <- 'celltype'; res$celltype <- cell_types
for (nomnom in cell_types) {
  print(nomnom)
  for (fdr in c(0.01, 0.05, 0.1, 0.2)) {
    res[[paste0("FDR ",fdr*100,"%")]][res$celltype == nomnom] <- sum(qvals$celltype == nomnom & qvals$fdr_permstim <= fdr)
  }
}
# write.table(res, "./single_cell_DEG_by_celltype.txt", row.names=F, col.names = T, sep="\t", quote=F)

write.table(qvals[order(qvals$gene, qvals$celltype),c(1,5,6,7,8,12,14)], "./single_cell_results.txt", row.names=F, col.names = T, sep="\t", quote=F)

save.image("./single_cell_with_permutations.RData")

ggplot(data=limma_stim, aes(y=B_Stimulation , x=celltype)) + geom_boxplot(fill='gray') + geom_hline(yintercept = 0, col='black') + theme_classic()
ggplot(data=limma_stim, aes(y=-log10(P.Value), x=celltype)) + geom_violin(fill='gray') + geom_hline(yintercept = 0, col='black') + theme_classic()
##-----------------------------------------------

##-----------------------------------------------
## PCA
## Does ERAP2 genotype map onto major PCs for any cell type or condition? 
##-----------------------------------------------
for (i in 1:length(cell_types)) {
  print(cell_types[i])
  #chose the data to produce the PCA
  target.pca.data = t(voom_data[[cell_types[i]]])
  # code to produce PCA factors
  pca.data = as.data.frame(prcomp(target.pca.data, scale. = T)$x)
  percentVar  = stats:::summary.prcomp(prcomp(target.pca.data))$importance
  # estimate the variance of the first two components
  PC1.variance = round(100 * percentVar[2,1])
  PC2.variance = round(100 * percentVar[2,2])
  PC3.variance = round(100 * percentVar[2,3])
  PC4.variance = round(100 * percentVar[2,4])
  
  # produce the data frame for plotting
  pca.data.4plot= data.frame(PC1 = pca.data$PC1, PC2 = pca.data$PC2, PC3 = pca.data$PC3, PC4 = pca.data$PC4, conditions = info$stim, genotype = as.factor(info$geno))
  ggplot(data = pca.data.4plot, aes(x = PC1, y = PC2, color = conditions, pch=genotype)) + 
    geom_point(size =5) + 
    xlab(paste0("PC1: ", PC1.variance, "% variance")) + 
    ylab(paste0("PC2: ", PC2.variance, "% variance")) + 
    theme_bw() + ggtitle(cell_types[i]) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=12), axis.title=element_text(size=14, face='bold'), legend.text = element_text(size=14), legend.title = element_text(size=14, face='bold'), legend.position='bottom') + scale_color_discrete(name='condition', labels=c("null", "YP"))
  
  ggplot(data = pca.data.4plot, aes(x = PC3, y = PC4, color = conditions, pch=genotype)) + 
    geom_point(size =5) + 
    xlab(paste0("PC3: ", PC3.variance, "% variance")) + 
    ylab(paste0("PC4: ", PC4.variance, "% variance")) + 
    theme_bw() + ggtitle(cell_types[i]) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=12), axis.title=element_text(size=14, face='bold'), legend.text = element_text(size=14), legend.title = element_text(size=14, face='bold'), legend.position='bottom') + scale_color_discrete(name='condition', labels=c("null", "YP"))
  
  
  print("1"); print(summary(lm(pca.data$PC1 ~ conditions*genotype, data=pca.data.4plot))$coefficients)
  print("2"); print(summary(lm(pca.data$PC2 ~ conditions*genotype, data=pca.data.4plot))$coefficients)
  print("3"); print(summary(lm(pca.data$PC3 ~ conditions*genotype, data=pca.data.4plot))$coefficients)
  print("4"); print(summary(lm(pca.data$PC4 ~ conditions*genotype, data=pca.data.4plot))$coefficients)
  print("5"); print(summary(lm(pca.data$PC5 ~ conditions*genotype, data=pca.data.4plot))$coefficients)
  print("6"); print(summary(lm(pca.data$PC6 ~ conditions*genotype, data=pca.data.4plot))$coefficients)
  print("7"); print(summary(lm(pca.data$PC7 ~ conditions*genotype, data=pca.data.4plot))$coefficients)
  print("8"); print(summary(lm(pca.data$PC8 ~ conditions*genotype, data=pca.data.4plot))$coefficients)
  print("9"); print(summary(lm(pca.data$PC9 ~ conditions*genotype, data=pca.data.4plot))$coefficients)
  print("10"); print(summary(lm(pca.data$PC10 ~ conditions*genotype, data=pca.data.4plot))$coefficients)
  
  rm(pca.data.4plot, PC1.variance, PC2.variance, pca.data, percentVar)
}

## per condition: not stimulated
for (i in 1:length(cell_types)) {
  print(cell_types[i])
  #chose the data to produce the PCA
  target.pca.data = t(voom_data[[cell_types[i]]][,info$stim == "0_unstim"])
  # code to produce PCA factors
  pca.data = as.data.frame(prcomp(target.pca.data, scale. = T)$x)
  percentVar  = stats:::summary.prcomp(prcomp(target.pca.data))$importance
  # estimate the variance of the first two components
  PC1.variance = round(100 * percentVar[2,1])
  PC2.variance = round(100 * percentVar[2,2])
  PC3.variance = round(100 * percentVar[2,3])
  PC4.variance = round(100 * percentVar[2,4])
  
  # produce the data frame for plotting
  pca.data.4plot= data.frame(PC1 = pca.data$PC1, PC2 = pca.data$PC2, PC3 = pca.data$PC3, PC4 = pca.data$PC4, genotype = as.factor(info$geno[info$stim == "0_unstim"]))
  ggplot(data = pca.data.4plot, aes(x = PC1, y = PC2, color = genotype, pch=genotype)) + 
    geom_point(size =5) + 
    xlab(paste0("PC1: ", PC1.variance, "% variance")) + 
    ylab(paste0("PC2: ", PC2.variance, "% variance")) + 
    theme_bw() + ggtitle(cell_types[i]) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=12), axis.title=element_text(size=14, face='bold'), legend.text = element_text(size=14), legend.title = element_text(size=14, face='bold'), legend.position='bottom') 
  
  ggplot(data = pca.data.4plot, aes(x = PC3, y = PC4, color = genotype, pch=genotype)) + 
    geom_point(size =5) + 
    xlab(paste0("PC3: ", PC3.variance, "% variance")) + 
    ylab(paste0("PC4: ", PC4.variance, "% variance")) + 
    theme_bw() + ggtitle(cell_types[i]) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=12), axis.title=element_text(size=14, face='bold'), legend.text = element_text(size=14), legend.title = element_text(size=14, face='bold'), legend.position='bottom') 
  
  
  print("1"); print(summary(lm(pca.data$PC1 ~ genotype, data=pca.data.4plot))$coefficients)
  print("2"); print(summary(lm(pca.data$PC2 ~ genotype, data=pca.data.4plot))$coefficients)
  print("3"); print(summary(lm(pca.data$PC3 ~ genotype, data=pca.data.4plot))$coefficients)
  print("4"); print(summary(lm(pca.data$PC4 ~ genotype, data=pca.data.4plot))$coefficients)
  print("5"); print(summary(lm(pca.data$PC5 ~ genotype, data=pca.data.4plot))$coefficients)
  print("6"); print(summary(lm(pca.data$PC6 ~ genotype, data=pca.data.4plot))$coefficients)
  print("7"); print(summary(lm(pca.data$PC7 ~ genotype, data=pca.data.4plot))$coefficients)
  print("8"); print(summary(lm(pca.data$PC8 ~ genotype, data=pca.data.4plot))$coefficients)
  print("9"); print(summary(lm(pca.data$PC9 ~ genotype, data=pca.data.4plot))$coefficients)
  print("10"); print(summary(lm(pca.data$PC10 ~ genotype, data=pca.data.4plot))$coefficients)
  
  rm(pca.data.4plot, PC1.variance, PC2.variance, pca.data, percentVar)
}

## per condition: stimulated
for (i in 1:length(cell_types)) {
  print(cell_types[i])
  #chose the data to produce the PCA
  target.pca.data = t(voom_data[[cell_types[i]]][,info$stim == "1_stim"])
  # code to produce PCA factors
  pca.data = as.data.frame(prcomp(target.pca.data, scale. = T)$x)
  percentVar  = stats:::summary.prcomp(prcomp(target.pca.data))$importance
  # estimate the variance of the first two components
  PC1.variance = round(100 * percentVar[2,1])
  PC2.variance = round(100 * percentVar[2,2])
  PC3.variance = round(100 * percentVar[2,3])
  PC4.variance = round(100 * percentVar[2,4])
  
  # produce the data frame for plotting
  pca.data.4plot= data.frame(PC1 = pca.data$PC1, PC2 = pca.data$PC2, PC3 = pca.data$PC3, PC4 = pca.data$PC4, genotype = as.factor(info$geno[info$stim == "0_unstim"]))
  ggplot(data = pca.data.4plot, aes(x = PC1, y = PC2, color = genotype, pch=genotype)) + 
    geom_point(size =5) + 
    xlab(paste0("PC1: ", PC1.variance, "% variance")) + 
    ylab(paste0("PC2: ", PC2.variance, "% variance")) + 
    theme_bw() + ggtitle(cell_types[i]) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=12), axis.title=element_text(size=14, face='bold'), legend.text = element_text(size=14), legend.title = element_text(size=14, face='bold'), legend.position='bottom') 
  
  ggplot(data = pca.data.4plot, aes(x = PC3, y = PC4, color = genotype, pch=genotype)) + 
    geom_point(size =5) + 
    xlab(paste0("PC3: ", PC3.variance, "% variance")) + 
    ylab(paste0("PC4: ", PC4.variance, "% variance")) + 
    theme_bw() + ggtitle(cell_types[i]) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=12), axis.title=element_text(size=14, face='bold'), legend.text = element_text(size=14), legend.title = element_text(size=14, face='bold'), legend.position='bottom') 
  
  
  print("1"); print(summary(lm(pca.data$PC1 ~ genotype, data=pca.data.4plot))$coefficients)
  print("2"); print(summary(lm(pca.data$PC2 ~ genotype, data=pca.data.4plot))$coefficients)
  print("3"); print(summary(lm(pca.data$PC3 ~ genotype, data=pca.data.4plot))$coefficients)
  print("4"); print(summary(lm(pca.data$PC4 ~ genotype, data=pca.data.4plot))$coefficients)
  print("5"); print(summary(lm(pca.data$PC5 ~ genotype, data=pca.data.4plot))$coefficients)
  print("6"); print(summary(lm(pca.data$PC6 ~ genotype, data=pca.data.4plot))$coefficients)
  print("7"); print(summary(lm(pca.data$PC7 ~ genotype, data=pca.data.4plot))$coefficients)
  print("8"); print(summary(lm(pca.data$PC8 ~ genotype, data=pca.data.4plot))$coefficients)
  print("9"); print(summary(lm(pca.data$PC9 ~ genotype, data=pca.data.4plot))$coefficients)
  print("10"); print(summary(lm(pca.data$PC10 ~ genotype, data=pca.data.4plot))$coefficients)
  
  rm(pca.data.4plot, PC1.variance, PC2.variance, pca.data, percentVar)
}

##-----------------------------------------------


##-----------------------------------------------
## model genotype-by-expression for each cell type
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
  data <- data[rowMeans(cpm) >= min(c(15, mean(cpm[row.names(cpm) == "ERAP2",]), mean(cpm[row.names(cpm) == "ERAP1",]))),]
  mean_cpm <- rowMeans(cpm[rowMeans(cpm) >= min(c(15, mean(cpm[row.names(cpm) == "ERAP2",]), mean(cpm[row.names(cpm) == "ERAP1",]))),])
  
  
  # add number of cells to our metadata 
  names$cells <- NA
  for (j in 1:nrow(names)) {
    names$cells[j] <- all_counts[row.names(all_counts) == names$V1[j], i]
  }; rm(j)
  
  # create design matrix
  stim <- names$stim; cells <- names$cells; samples <- as.factor(names$sample_name)
  gt <- names$geno
  design <- model.matrix(~ 1 + stim:gt + cells)  
  colnames(design)[1:4] <- c("intercept", "cells", "gt_unstim", "gt_stim")
  rm(stim, cells, samples, gt)
  # + pop:stimulated , data = exp
  
  # voom normalize, then model using lmfit 
  v <- voom(data, design, plot=F)
  fit = lmFit(v)
  
  contr.matrix <- makeContrasts(
    genotype = gt_stim,
    levels=colnames(design)
  )
  fit2 <- contrasts.fit(fit = fit, contrasts = contr.matrix)
  fit2 = eBayes(fit2)
  # extract summary table + cell type + gene names
  res_stim <- topTreat(fit2, coef=1, n=Inf)
  res_stim$test <- "gt_stim"
  res_stim$celltype <- cell_types[i]
  res_stim$gene <- rownames(res_stim)
  SE <- sqrt(fit2$s2.post) * fit2$stdev.unscaled; colnames(SE) <- paste("SE_", colnames(SE), sep="")
  B22 <- fit2$coefficients; colnames(B22) <- paste("B_", colnames(B22), sep="")
  output <- as.data.frame(cbind(fit2$p.value, B22, SE)); rm(SE, B22)
  output$gene <- rownames(output)
  output <- merge(output, res_stim, by="gene")
  output$cpm <- mean_cpm
  limma_geno <- rbind(limma_geno, output); rm(output)
  
  contr.matrix <- makeContrasts(
    genotype = gt_unstim,
    levels=colnames(design)
  )
  fit2 <- contrasts.fit(fit = fit, contrasts = contr.matrix)
  fit2 = eBayes(fit2)
  # extract summary table + cell type + gene names
  res_stim <- topTreat(fit2, coef=1, n=Inf)
  res_stim$test <- "gt_null"
  res_stim$celltype <- cell_types[i]
  res_stim$gene <- rownames(res_stim)
  SE <- sqrt(fit2$s2.post) * fit2$stdev.unscaled; colnames(SE) <- paste("SE_", colnames(SE), sep="")
  B22 <- fit2$coefficients; colnames(B22) <- paste("B_", colnames(B22), sep="")
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
tmp_data <- limma_stim[limma_stim$celltype == unique(limma_stim$celltype)[i],]
new_tmp <- tmp_data %>%
  mutate(quantile = ntile(cpm,10))
ggplot(new_tmp, aes(x=quantile, y=B_Stimulation , group=quantile)) + 
  geom_violin() + geom_boxplot(width=0.2) + 
  geom_hline(yintercept=0, color="grey") + 
  ggtitle(paste("quantiles ranked by cpm\ncell type:",unique(limma_geno$cell_type)[i],sep="")) + theme_classic()
print(cor.test(new_tmp$`B_Stimulation`, new_tmp$cpm, method='spearman'))
rm(i, new_tmp, tmp_data)
##-----------------------------------------------

##-----------------------------------------------
# Plot target genes
##-----------------------------------------------

sig <- c(
  "CTLA4", "ICOS", # chr2
  "ERAP1", "ERAP2", "LNPEP", # chr5
  "TICAM2", "TMED7", # chr5
  "NFATC1" # chr18    
)

fcs <- NULL; for (tmp_gene in sig[4]) {
  print(tmp_gene)
  print(limma_stim[limma_stim$gene == tmp_gene,c(1,5,8,10:13)])
  tmp_fcs <- NULL
  for (i in c(1:5)) {
    nom <- cell_types[i]
    tmp_data <- voom_data[[nom]]
    tmp_data <- tmp_data[row.names(tmp_data) == tmp_gene,]
    if (length(tmp_data) > 0) {
      for (tmp_sample in unique(info$sample_name)) {
        this_fc <- as.data.frame(matrix(nrow=1,ncol=4))
        colnames(this_fc) <- c("name", "celltype", "stim", "unstim")
        this_fc$name <- tmp_sample
        this_fc$celltype <- nom
        this_fc$stim <- tmp_data[info$sample_name == tmp_sample & info$stim == "1_stim"]
        this_fc$unstim <- tmp_data[info$sample_name == tmp_sample & info$stim == "0_unstim"]
        rbind(this_fc, tmp_fcs) -> tmp_fcs; rm(this_fc)
      }; rm(tmp_sample)
    }}
  rbind(tmp_fcs,fcs) -> fcs; rm(tmp_fcs)
}

ggplot(data=fcs, aes(x=celltype, y=log2(stim/unstim))) + 
  geom_violin() + 
  geom_jitter(width=0.2) + 
  geom_hline(yintercept = 0, col='gray') + 
  ggtitle("ERAP2") +
  stat_summary(fun.data='mean_cl_boot', fun.args=list(conf.int=.95), geom='pointrange', col='red', size=1) +
  theme_classic()

for (tmp_gene in sig) {
  print(tmp_gene)
  print(limma_stim[limma_stim$gene == tmp_gene,])
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

limma_stim <- limma_stim[limma_stim$gene %in% c("ERAP1", "ERAP2", "LNPEP", "TICAM2", "TMED7", "CTLA4", "NFATC1", "ICOS"),]
limma_stim$gene2[limma_stim$gene %in% c("ERAP1", "ERAP2", "LNPEP")] <- paste("1", limma_stim$gene[limma_stim$gene %in% c("ERAP1", "ERAP2", "LNPEP")], sep="_")
limma_stim$gene2[limma_stim$gene %in% c("TICAM2", "TMED7")] <- paste("2", limma_stim$gene[limma_stim$gene %in% c("TICAM2", "TMED7")], sep="_")
limma_stim$gene2[limma_stim$gene %in% c("ICOS","CTLA4")] <- paste("3", limma_stim$gene[limma_stim$gene %in% c("ICOS","CTLA4")], sep="_")
limma_stim$gene2[limma_stim$gene %in% c("NFATC1")] <- paste("4", limma_stim$gene[limma_stim$gene %in% c("NFATC1")], sep="_")


panel_stimulation_by_celltype <- ggplot(limma_stim, aes(y=factor(gene2), x=factor(celltype))) +
  geom_point(aes(color=logFC, size= -log10(P.Value)), shape=16) +
  scale_color_gradient2(low='purple4', mid='snow2', high='tomato2') + 
  scale_size(range=c(5,20)) + ylab("gene") + xlab("celltype") + 
  geom_hline(yintercept = c(0.5,1.5,2.5,3.5), size = .2) +
  scale_y_discrete(limits=rev) + scale_x_discrete(position = "top")  +
  geom_text(aes(label=signif(logFC, digits = 2)), nudge_y = 0.1) +
  geom_text(aes(label=paste0("p = ",signif(P.Value, digits = 2))), nudge_y = -0.1) +
  theme_classic() + theme(legend.position = "bottom", 
                          axis.text=element_text(color='black'), text=element_text(size=16)
  )

limma_geno$test <- gsub("genotype", "rs2549794 genotype", limma_geno$test)
tmp <- limma_geno[limma_geno$gene %in% c("ERAP2", "ERAP1", "LNPEP"),]
tmp$logFC <- -tmp$logFC
panel_genotype_by_celltype <- ggplot(tmp, aes(y=factor(gene), x=factor(celltype))) +
  geom_point(aes(color=logFC, size= -log10(P.Value)), shape=16) +
  scale_color_gradient2(low='purple4', mid='snow2', high='tomato2') + 
  scale_size(range=c(5,20)) + ylab("gene") + xlab("celltype") + 
  geom_hline(yintercept = c(0.5), size = .2) +
  scale_y_discrete(limits=rev) + scale_x_discrete(position = "top")  +
  geom_text(aes(label=signif(logFC, digits = 2)), nudge_y = 0.05) +
  geom_text(aes(label=paste0("p = ",signif(P.Value, digits = 2))), nudge_y = -0.05) +
  theme_classic() + theme(legend.position = "bottom", 
                          axis.text=element_text(color='black'), text=element_text(size=16)
  )

rm(i, logCPM_filter, voom_data, tmp, target.pca.data, pseudobulk_means, num_genes_out, info, all_counts, nom, PC3.variance, PC4.variance, sample_colname, sig, tmp_data, tmp_gene, working_dir)

save.image("~/single_cell_images.RData")


#################
## eQTL per cell type, boxplot
#################
noms <- c(
  '0_unstim' = 'Not-Infected', 
  '1_stim' = 'Infected'
)
eQTL <- NULL; for (i in 1:length(cell_types)) {
  tmp_type <- cell_types[i]
  tmp_data <- voom_data[[tmp_type]]
  tmp_data <- tmp_data[row.names(tmp_data) == "ERAP2",]
  to_plot <- cbind(info, tmp_data)
  
  plot <- ggplot(data=to_plot[to_plot$stim == "0_unstim",], aes(y=tmp_data, x=as.factor(geno), fill=stim)) + 
    # facet_grid(. ~ stim, labeller = as_labeller(noms)) + 
    geom_violin() + 
    geom_jitter(width=0.2, size=2, col='gray42') + 
    stat_summary(fun.data='mean_cl_boot', geom='pointrange', col='red', size=1, show.legend = F) +
    scale_fill_manual(values=c("lightgoldenrodyellow", "lightgoldenrod2"),
                      name="", breaks=c("0_unstim", "1_stim"),
                      labels=c("Not-Infected", "Infected")) + 
    theme_classic() +
    ylab(label = paste0("normalized expression")) + #ggtitle(paste0(tmp_type))+ 
    xlab(label=paste("genotype: alternate alleles",sep="")) +
    theme(legend.position = "bottom", axis.title.x = element_blank()) + 
    theme(strip.background = element_rect(color='NA', fill='NA', size=1.5)); plot
  
  assign(paste("eQTL.",i,sep=""), plot)
}; rm(i, tmp_type)

library(patchwork); eQTL_panel <- eQTL.4 / eQTL.3 / eQTL.2 / eQTL.5 / eQTL.1 + plot_layout(guides="collect") & theme(legend.position = 'none')

save(eQTL_panel, file="./eQTL_panel.RData")


load("~/Barreiro/KlunkVilgalys_2021/pestis_aDNA/fig3D.rdata")
load("~/Barreiro/KlunkVilgalys_2021/pestis_aDNA/fig3E.rdata")

# 16" wide by 6.25" tall
fig3d + fig3e + eQTL_panel + plot_layout(widths=c(3,6,1))
