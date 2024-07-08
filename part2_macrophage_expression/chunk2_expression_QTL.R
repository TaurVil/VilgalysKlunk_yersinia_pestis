library(data.table)
library(DescTools)
library(edgeR)
library(stringr)
library(ggplot2)
library(ggridges)
library(ggrepel)
library(superheat)

##--------------------------------------
# assemble genotypes and filter for target variants
##--------------------------------------

# Genotype data for these samples was previously published in Nedelec et al. We've included a condensed version here consisting of the genotypes at our candidate variants from the ancient DNA. 

library(data.table); library(DESeq2); library(ggplot2); library(limma); library(patchwork)
load(("./DATA/eQTL_candidates.RData"))
load("./DATA/n33_pestis_infection.RData")
sig <- c(
  # "CTLA4", "ICOS", # chr2
  # "TICAM2", "TMED7", "TMED7-TICAM2" # chr5
  # "NFATC1", # chr18    
  "ERAP1", "ERAP2", "LNPEP" # chr5
)
sig_genes <- subset(res_lmfit, res_lmfit$Symbol %in% sig)

## Redefine expression based on the test we want to use
method <- "CPM"
if (method == "voom") { exp <- as.matrix(v$E) }
if (method == "CPM") { exp <- as.matrix(log2(cpm+0.01)) }
if (nrow(exp) == nrow(genes)) {genes$Gene -> row.names(exp)}

for (i in 1:length(sig)) {
  # if (i %in% seq(1,2)) {tmp_rsid <- "rs11571319"}
  # if (i %in% seq(6,8)) {tmp_rsid <- "rs17473484"}
  # if (i %in% 9) {tmp_rsid <- "rs73973415"}
  # 
  tmp_gene_name <- sig[i]
  if (tmp_gene_name %in% c("ERAP1", "ERAP2", "LNPEP")) {tmp_rsid <- "rs2549794"}
  
  print(paste("starting ", tmp_gene_name, "      \n", sep=""))
  tmp_gene_name2 <- genes$Gene[genes$Symbol == tmp_gene_name]
  
  tmp_gt <- gt[snps$V1 == tmp_rsid,]
  
  tmp_data <- cbind(names, exp[rownames(exp) == tmp_gene_name2,])
  tmp_data$snp_name <- paste(tmp_data$individual, tmp_data$individual, sep="_")
  colnames(tmp_data)[6] <- 'expression'
  tmp_data$gt <- NA; for (j in 1:nrow(tmp_data)) {
    if (tmp_data$snp_name[j] %in% snps_names$V1) {tmp_data$gt[j] <- tmp_gt[snps_names$V1 == tmp_data$snp_name[j]]}
  }; rm(j)
  
  to_plot <- tmp_data[!is.na(tmp_data$gt),]
  
  to_plot$gt <- as.numeric(to_plot$gt)
  
  print(summary(lm(data=to_plot, expression ~ gt*stimulated))$coefficients)
  print(summary(lm(data=to_plot, expression ~ gt+stimulated))$coefficients)
  print(summary(lm(data=to_plot, expression ~ stimulated:gt))$coefficients)
  
  # This was to adjust reference/alternate for ERAP2, but let's just do it as 012 across the board
  # to_plot$gt <- as.character(as.factor(to_plot$gt))
  # to_plot$gt[to_plot$gt == 0] <- '0_derived'
  # to_plot$gt[to_plot$gt == 1] <- '1_heterozygous'
  # to_plot$gt[to_plot$gt == 2] <- '2_ancestral'
  to_plot$gt <- as.factor(to_plot$gt)
  
  plot <- ggplot(data=to_plot, aes(y=expression, x=gt, fill=stimulated)) + 
    #geom_violin() + 
    geom_boxplot(width=0.25) + 
    scale_fill_manual(values=c("#999999", "#E69F00"),
                      name="",                                                                                                                               breaks=c("FALSE", "TRUE"),
                      labels=c("Not-Infected", "Infected")) + 
    theme_classic() +
    #scale_fill_manual(values=c("#999999", "#E69F00"), name="Population",breaks=c("FALSE", "TRUE"),labels=c("African", "European")) + 
    ylab(label = "normalized CPM") + ggtitle(tmp_gene_name) + xlab(label=paste(tmp_rsid, ": alternate alleles",sep="")) +
    theme(legend.position = "bottom")
  assign(paste("QTL.",i,sep=""), plot)
}

QTL.1 + QTL.2 + QTL.3 + # QTL.4 +  QTL.5 + QTL.6 + QTL.7 + QTL.8 + 
 plot_layout(ncol=4) + 
  plot_layout(guides="collect") & theme(legend.position = 'bottom')

rownames(exp) <- genes$Gene

# par(mfrow=c(1,3))
# exp1 <- as.matrix(v$E)[v$genes$Gene == "ENSG00000164307",]
# exp2 <- as.matrix(exp)[rownames(exp) == "ENSG00000164307",]
# exp3 <- as.matrix(cpm)[rownames(cpm) == "ENSG00000164307",]
# plot(exp1 ~ exp2, xlab="exp", ylab="voom", col=as.factor(names$stimulated), pch=16, main='ERAP1'); abline(a=0,b=1,col='blue')
# 
# exp1 <- as.matrix(v$E)[v$genes$Gene == "ENSG00000243414",]
# exp2 <- as.matrix(exp)[rownames(exp) == "ENSG00000243414",]
# exp3 <- as.matrix(cpm)[rownames(cpm) == "ENSG00000243414",]
# plot(exp1 ~ exp2, xlab="exp", ylab="voom", col=as.factor(names$stimulated), pch=16, main="TICAM2"); abline(a=0,b=1,col='blue')
# 
# exp1 <- as.matrix(v$E)[v$genes$Gene == "ENSG00000164308",]
# exp2 <- as.matrix(exp)[rownames(exp) == "ENSG00000164308",]
# exp3 <- as.matrix(cpm)[rownames(cpm) == "ENSG00000164308",]
# plot(exp1 ~ exp2, xlab="exp", ylab="voom", col=as.factor(names$stimulated), pch=16, main="ERAP2"); abline(a=0,b=1,col='blue')

noms <- c(
  'FALSE' = 'Not-Infected', 
  'TRUE' = 'Infected'
)

## split each into a NI and YP panel, for the 2 we want
eQTL <- NULL; eQTL_data <- NULL; for (i in c(1:3)) {
  # if (i %in% seq(1,2)) {tmp_rsid <- "rs11571319"}
  # if (i %in% seq(6,8)) {tmp_rsid <- "rs17473484"}
  # if (i %in% 9) {tmp_rsid <- "rs73973415"}
  
  tmp_gene_name <- sig[i]
  if (tmp_gene_name %in% c("ERAP1", "ERAP2", "LNPEP")) {tmp_rsid <- "rs2549794"}
  
  print(paste("starting ", tmp_gene_name, sep=""))
  tmp_gene_name2 <- genes$Gene[genes$Symbol == tmp_gene_name]
  
  tmp_gt <- gt[snps$V1 == tmp_rsid,]
  
  tmp_data <- cbind(names, exp[rownames(exp) == tmp_gene_name2,])
  tmp_data$snp_name <- paste(tmp_data$individual, tmp_data$individual, sep="_")
  colnames(tmp_data)[6] <- 'expression'
  tmp_data$gt <- NA; for (j in 1:nrow(tmp_data)) {
    if (tmp_data$snp_name[j] %in% snps_names$V1) {tmp_data$gt[j] <- tmp_gt[snps_names$V1 == tmp_data$snp_name[j]]}
  }; rm(j)
  
  to_plot <- tmp_data[!is.na(tmp_data$gt),]
  to_plot$gene <- tmp_gene_name2
  
  eQTL_data <- rbind(eQTL_data, to_plot[,c(1,4,3,6,8,9)])
  
  to_plot$gt <- as.factor(to_plot$gt)
  to_plot$gene <- tmp_gene_name
  
  plot <- ggplot(data=to_plot, aes(y=expression, x=gt, fill=stimulated)) + 
    facet_grid(. ~ stimulated, labeller = as_labeller(noms)) + 
    geom_violin() + 
    #geom_boxplot(width=0.25) + 
    geom_jitter(width=0.2, size=2, col='gray42') + 
    stat_summary(fun.data='mean_cl_boot', geom='pointrange', col='red', size=1) +
    scale_fill_manual(values=c("lightgoldenrodyellow", "lightgoldenrod2"),
                      name="", breaks=c("FALSE", "TRUE"),
                      labels=c("Not-Infected", "Infected")) + 
    theme_classic() +
    #scale_fill_manual(values=c("#999999", "#E69F00"), name="Population",breaks=c("FALSE", "TRUE"),labels=c("African", "European")) + 
    ylab(label = "normalized CPM") + ggtitle(tmp_gene_name) + xlab(label=paste(tmp_rsid, ": alternate alleles",sep="")) +
    theme(legend.position = "bottom") + 
    theme(strip.background = element_rect(color='NA', fill='NA', linewidth=1.5)); plot
  
  mod <- summary(lm(to_plot$expression ~ to_plot$stimulated:as.numeric(to_plot$gt)))
  res <- as.data.frame(matrix(nrow=2,ncol=5)); 
  colnames(res) <- c("gene", "condition", "snp", "P.Value", "Beta")
  res$gene <- tmp_gene_name
  res$condition <- c("NI", "pestis")
  res$snp <- tmp_rsid
  res$P.Value <- mod$coefficients[2:3,4]
  res$Beta <- mod$coefficients[2:3,1]
  rbind(eQTL, res) -> eQTL
  
  assign(paste("eQTL.",i,sep=""), plot)
}; rm(i)

fig3c <- eQTL.2 + plot_layout(guides="collect") & theme(legend.position = 'none'); fig3c
save(fig3c, file="./fig3C.rdata")
# save(eQTL, file="~/qtl_pestis.rdata")
# save(eQTL_data, file="~/pestis_eqtl_s8.rdata")
