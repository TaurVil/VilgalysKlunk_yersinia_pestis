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
  "CTLA4", "ICOS", # chr2
  "ERAP1", "ERAP2", "LNPEP", # chr5
  "TICAM2", "TMED7", "TMED7-TICAM2"#, # chr5
  # "NFATC1" # chr18    
)
sig_genes <- subset(res_lmfit, res_lmfit$Symbol %in% sig)

## Redefine expression based on the test we want to use
method <- "voom"
if (method == "voom") { exp <- as.matrix(v$E) }
if (method == "CPM") { exp <- as.matrix(log2(cpm+0.01)) }
if (nrow(exp) == nrow(genes)) {genes$Gene -> row.names(exp)}

for (i in 1:length(sig)) {
  if (i %in% seq(1,2)) {tmp_rsid <- "rs11571319"}
  if (i %in% seq(3,5)) {tmp_rsid <- "rs2549794"}
  if (i %in% seq(6,8)) {tmp_rsid <- "rs17473484"}
  if (i %in% 9) {tmp_rsid <- "rs73973415"}
  
  tmp_gene_name <- sig[i]
  print(paste("\nstarting ", tmp_gene_name, "\n", sep=""))
  tmp_gene_name2 <- genes$Gene[genes$Symbol == tmp_gene_name]
  
  tmp_gt <- gt[snps$V1 == tmp_rsid,]
  
  tmp_data <- cbind(names, exp[rownames(exp) == tmp_gene_name2,])
  tmp_data$snp_name <- paste(tmp_data$individual, tmp_data$individual, sep="_")
  colnames(tmp_data)[6] <- 'expression'
  tmp_data$gt <- NA; for (j in 1:nrow(tmp_data)) {
    if (tmp_data$snp_name[j] %in% snps_names$V1) {tmp_data$gt[j] <- tmp_gt[snps_names$V1 == tmp_data$snp_name[j]]}
  }; rm(j)
  
  to_plot <- tmp_data[!is.na(tmp_data$gt),]
  
  print(summary(lm(to_plot$expression ~ as.numeric(to_plot$gt)*to_plot$stimulated))$coefficients)
  print(summary(lm(to_plot$expression ~ as.numeric(to_plot$gt)+to_plot$stimulated))$coefficients)
  
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

QTL.1 + QTL.2 + QTL.3 + QTL.4 + 
  QTL.5 + QTL.6 + QTL.7 + QTL.8 + plot_layout(ncol=4) + 
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
for (i in c(4,6)) {
  if (i %in% seq(1,2)) {tmp_rsid <- "rs11571319"}
  if (i %in% seq(3,5)) {tmp_rsid <- "rs2549794"}
  if (i %in% seq(6,8)) {tmp_rsid <- "rs17473484"}
  if (i %in% 9) {tmp_rsid <- "rs73973415"}
  
  tmp_gene_name <- sig[i]
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
    theme(strip.background = element_rect(color='NA', fill='NA', size=1.5)); plot
  
  assign(paste("eQTL.",i,sep=""), plot)
}; rm(i)

fig3c <- eQTL.4 / eQTL.6 + plot_layout(guides="collect") & theme(legend.position = 'none'); fig3c
save(fig3c, file="./fig3C.rdata")
