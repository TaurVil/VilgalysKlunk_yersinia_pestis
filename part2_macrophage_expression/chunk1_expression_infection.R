# n=33 individuals, macrophages
##--------------------------------------
# read in expression data
##--------------------------------------
library(data.table); library(plyr); library(dplyr); library(data.table);library(DescTools)
library(edgeR); library(stringr); library(ggplot2); library(ggridges); library(ggrepel); library(superheat)
exp <- fread("./DATA/rawCountMatrix.csv")
exp[,1:2] -> genes; exp[,-(1:2)] -> exp
##--------------------------------------

##--------------------------------------
# get metadata
##--------------------------------------
names <- colnames(exp); names <- as.data.frame(names, ncol=1); rownames(exp) <- genes$Gene
library(stringr); names$pop <- str_detect(names$names, "EU"); names$stimulated <- str_detect(names$names, "YP")
strsplit(as.character(names$names),"_") -> V1; ldply(V1) -> V1; V1[,1] -> names$individual; rm(V1)
names$reads <- colSums(exp)
names$pop[names$pop == F] <- "AFR"; names$pop[names$pop == T] <- "EUR"
pop <- as.factor(names$pop); stimulated <- as.factor(names$stimulated); individual <- factor(names$individual)
##--------------------------------------

##--------------------------------------
# filter based on cpm
# an average of at least 1 cpm in either stimulated or null data
##--------------------------------------
cpm <- exp; for (i in 1:66) {cpm[,i] <- (exp[,..i]*1e6)/names$reads[i]}; rm(i)
as.data.frame(cpm) -> cpm
genes$cpm_null <- rowMeans(cpm[,c(names$stimulated == F)])
genes$cpm_stim <- rowMeans(cpm[,c(names$stimulated == T)])
rownames(cpm) <- rownames(exp)
# hist(log10(genes$cpm_null), breaks=50)
# plot(genes$cpm_null ~ genes$cpm_stim); abline(a=0,b=1,col='red')
genes$no.reads <- rowMeans(exp == 0)
genes <- subset(genes, (genes$cpm_null >= 1 | genes$cpm_stim >= 1) & genes$no.reads < 0.5)
exp<-subset(exp, rownames(exp) %in% genes$Gene); 
cpm <- subset(cpm, rownames(cpm) %in% genes$Gene)

# genes <- genes[!apply(cpm, 1, max) > 5*apply(cpm, 1, mean),]
# exp <- exp[!apply(cpm, 1, max) > 5*apply(cpm, 1, mean),]

genes$log2 <- log2(genes$cpm_stim/genes$cpm_null)
as.data.frame(exp) -> exp
##--------------------------------------

##--------------------------------------
# model data using:
# model as a function of stimulation while controlling for individual (modeling population as well is inappropriate since we're controlling for individual)
# fit with lmFit
# extract DE genes using eBayes contrast
##--------------------------------------
rownames(exp) <- genes$Gene

indiv_nams <- names$individual
design <- model.matrix(~ 1 + stimulated + indiv_nams)
colnames(design)[1:2] <- c("intercept", "stimulated")

v = voom(exp, design, plot = T, save.plot = T)
##--------------------------------------


##--------------------------------------
# PCA
# samples in row 32 and 33 seem to be oddly placed; let's exclude them
# they don't just appear odd due to low library size
##--------------------------------------
#chose the data to produce the PCA
target.pca.data = t(v$E)
# code to produce PCA factors
pca.data = as.data.frame(prcomp(target.pca.data, scale. = T)$x)
percentVar  = stats:::summary.prcomp(prcomp(target.pca.data))$importance

# estimate the variance of the first two components
PC1.variance = round(100 * percentVar[2,1])
PC2.variance = round(100 * percentVar[2,2])

# produce the data frame for plotting
pca.data.4plot= data.frame(PC1 = pca.data$PC1, PC2 = pca.data$PC2, Species = pop, conditions = stimulated)
ggplot(data = pca.data.4plot, aes(x = PC1, y = PC2, color = conditions)) + 
  geom_point(size =5, pch=20) + 
  xlab(paste0("PC1: ", PC1.variance, "% variance")) + 
  ylab(paste0("PC2: ", PC2.variance, "% variance")) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=12), axis.title=element_text(size=14, face='bold'), legend.text = element_text(size=14), legend.title = element_text(size=14, face='bold'), legend.position='bottom') + scale_color_discrete(name='condition', labels=c("null", "Y. pestis"))
  
rm(pca.data.4plot, PC1.variance, PC2.variance, pca.data, percentVar)
##--------------------------------------

# ##--------------------------------------
# # remodel data excluding the individual who doesn't group well in PCA -- did this to check whether the results change, they don't
# ##--------------------------------------
# exp <- exp[,-c(15,32:34)]
# names <- names[-c(15,32:34),]
# pop <- pop[-c(15,32:34)]
# stimulated <- stimulated[-c(15,32:34)]
# individual <- individual[-c(15,32:34)]
# 
# y <- DGEList(counts = exp, lib.size = colSums(exp), samples = names[,2:3], genes = genes)
# y <- calcNormFactors(y, method="TMM") #accounts for differences in library size
# # model with population and stimulation
# design <- model.matrix(~ 1 + stimulated)
# colnames(design) <- c("intercept", "stimulated")
# 
# # voom normalization then using duplicatedcorrelation function to take account of paired nature of the experiment then fitting in a fit object
# v = voom(y, design, plot = T)
# cor = duplicateCorrelation(v, design, block = names$individual)
# v = voom(y, design, plot=T, correlation = cor$consensus)
# cor = duplicateCorrelation(v, design, block = names$individual)
##--------------------------------------

# model with lmfit 
##--------------------------------------
fit = lmFit(v, design)

contr.matrix <- makeContrasts(
  Stimulation = stimulated,
  levels=colnames(design)
)
fit2 <- contrasts.fit(fit = fit, contrasts = contr.matrix)
fit2 = eBayes(fit2, robust=T)
res_lmfit <- topTreat(fit2, coef=1, n=Inf)
boxplot(res_lmfit$logFC); abline(h=0,col='black')
rm(target.pca.data, fit, fit2, contr.matrix, design, indiv_nams)
##--------------------------------------

#####
# rm(res_stim, dds)
# save.image("./DATA/n33_pestis_infection.RData")

library(data.table); library(DESeq2); library(ggplot2); library(limma); library(tidyverse)
setwd("~//GitHub/VilgalysKlunk_yersinia_pestis/part2_macrophage_expression")
load("./DATA/n33_pestis_infection.RData")
load("~/Barreiro/2022_KlunkVilgalys/time_series_pestis.RData")

## correlation between time series and original data
t1 <- time_series_results[time_series_results$test == "x4h",]
t2 <- res_lmfit
t2$gene <- row.names(t2)
t2 <- merge(t2, genes, by.x="gene", by.y="Gene")
colnames(t2)[2] <- "n33_original"
colnames(t1)[1] <- "n8_live"
tmp <- merge(t1, t2, by.x="gene", by.y="Symbol")

summary(lm(tmp$n8_live ~ tmp$n33_original))


l_vs_hk_1 <- ggplot(tmp, aes(x=n33_original, y=n8_live)) +
  geom_point(alpha=0.3, color='black') +
  # geom_density_2d(geom=alpha=0.4, contour_var="ndensity", bins=20) +
  geom_density_2d_filled(alpha=0.7,bins = 9) +
  scale_fill_brewer() + xlab("effect size: heat-killed") + ylab("effect size: live bacteria") +
  geom_smooth(method='lm') +
  # geom_point(data=tmp[tmp$gene %in% sig,], col='orange3', size=4) +
  theme_classic() + theme(legend.position ='none'); l_vs_hk_1

mean(abs(tmp$n8_live))
mean(abs(tmp$n33_original))

#####
# exp_raw <- as.matrix(log2(cpm+0.01))
exp_raw <- as.matrix(v$E)
# exp_raw <- exp; for (i in 1:66) {exp_raw[,i] <- (exp[,i]*1e6)/colSums(exp)[i]}; rm(i); exp_raw <- as.matrix(log2(cpm+0.01))
# rownames(exp_raw) <- genes$Symbol
# 
sum(res_lmfit$adj.P.Val < 0.01)
res_lmfit$gene <- rownames(res_lmfit)
  
res <- merge(res_lmfit, genes, by.x = 'gene', by.y="Gene")
sig <- c(
  # "CTLA4", "ICOS", # chr2
  # "TICAM2", "TMED7", "TMED7-TICAM2", # chr5
  # "NFATC1", # chr18  
  "ERAP1", "ERAP2", "LNPEP" # chr5
  )
tmp <- subset(res, res$Symbol %in% sig)
tmp[order(tmp$logFC),]

## expression of the 8 target genes we analyze
fcs <- NULL
for (tmp_num in 1:length(tmp$Symbol)) {
  tmp_gene_name <- tmp$Symbol[tmp_num]
  print(paste("starting:", tmp_gene_name, sep=" "))
  tmp_data <- cbind(names, exp_raw[genes$Symbol == tmp_gene_name,], t(cpm[genes$Symbol == tmp_gene_name,]))
  colnames(tmp_data)[6] <- 'expression'
  colnames(tmp_data)[7] <- 'cpm'
  
  tmp_data$pop[tmp_data$pop == F] <- 'African'; tmp_data$pop[tmp_data$pop == T] <- 'European'
  tmp_data$stimulated[tmp_data$stimulated == F] <- 'control'; tmp_data$stimulated[tmp_data$stimulated == T] <- 'infected'
  
  t1 <- tmp_data[tmp_data$stimulated == "control",]; colnames(t1)[7] <- 'NIcpm'
  t2 <- tmp_data[tmp_data$stimulated == "infected",]; colnames(t2)[7] <- 'YPcpm'
  t3 <- merge(t1, t2, by = "individual")
  t3$fc <- log2(((t3$YPcpm+0.05)/(t3$NIcpm+0.05)))
  # boxplot(t3$fc); abline(h=0,col='black')
  t3$gene <- tmp_gene_name
  rbind(fcs, t3) -> fcs; rm(t1, t2, t3)
  
  print(summary(lm(tmp_data$expression ~ tmp_data$stimulated + tmp_data$individual))$coefficients[1:3,])
  model <- summary(lm(tmp_data$expression~tmp_data$stimulated + tmp_data$pop));# print(model$coefficients)
  tmp_data$no_pop <- tmp_data$expression - (as.numeric(as.factor(tmp_data$pop))-1) * model$coefficients[3,1]; rm(model)
  
  plot <- ggplot(data=tmp_data, aes(y=expression, x=stimulated)) + 
    geom_violin(fill='lightgray') + 
    geom_boxplot(width=0.2, fill='gold3') + coord_cartesian(ylim=c(min(tmp_data$expression),max(tmp_data$expression))) +
    xlab(label = "infected") + ylab(label = "normalized CPM") + ggtitle(tmp_gene_name) + theme_classic()
  assign(x = paste("plot.",tmp_gene_name, sep=""), value = plot)
  
  rm(tmp_gene_name, tmp_data, plot)
}
rm(tmp_num)

exp_raw <- as.data.frame(time_series_cpm)

fcs_time <- NULL
for (tmp_num in 1:length(sig[sig %in% time_series_genes$genes])) {
  tmp_gene_name <- sig[sig %in% time_series_genes$genes][tmp_num]
  print("starting:"); print(tmp_gene_name)
  tmp_data <- cbind(time_series_metadata, t(exp_raw[row.names(exp_raw) == tmp_gene_name,]))
  colnames(tmp_data)[7] <- 'expression'
  
  tmp_data <- subset(tmp_data, tmp_data$Time %in% c(0,4))
  tmp_data$stimulated <- NA
  tmp_data$stimulated[tmp_data$Time == 0] <- 'control'; tmp_data$stimulated[tmp_data$Time == 4] <- 'infected'
  
  t1 <- tmp_data[tmp_data$stimulated == "control",]; colnames(t1)[7] <- 'NIcpm'
  t2 <- tmp_data[tmp_data$stimulated == "infected",]; colnames(t2)[7] <- 'YPcpm'
  t3 <- merge(t1, t2, by = "indiv_number")
  t3$fc <- log2(((t3$YPcpm+0.05)/(t3$NIcpm+0.05)))
  # boxplot(t3$fc); abline(h=0,col='black')
  t3$gene <- tmp_gene_name
  rbind(fcs_time, t3) -> fcs_time; rm(t1, t2, t3)

  rm(tmp_gene_name, tmp_data)
}
rm(tmp_num)

fcs <- fcs[,c(1:2,7,13,14,15)]
fcs_time <- fcs_time[,c(1:2,7,14,16,17)]

colnames(fcs_time) <- colnames(fcs)
fcs_time$analysis <- 'live'
fcs$analysis <- 'heat_killed'
fcs <- rbind(fcs, fcs_time)

library(patchwork)
#8.5 by 4.5 inches
(plot.LNPEP + plot.ERAP1 + plot.ERAP2) + plot_layout(ncol=4)

# (plot.TICAM2 + plot.TMED7 + `plot.TMED7-TICAM2` + plot.CTLA4 +
    # plot.ICOS + plot.LNPEP + plot.ERAP1 + plot.ERAP2) + plot_layout(ncol=4)

# fcs$fc[fcs$YPcpm == 0] <- -6
# fcs$fc[fcs$NIcpm == 0] <- 5
# fcs$fc[fcs$NIcpm == 0 & fcs$YPcpm == 0] <- 0

fcs$gene2[fcs$gene %in% c("ERAP1", "ERAP2", "LNPEP")] <- paste("1", fcs$gene[fcs$gene %in% c("ERAP1", "ERAP2", "LNPEP")], sep="_")
# fcs$gene2[fcs$gene %in% c("TICAM2", "TMED7")] <- paste("2", fcs$gene[fcs$gene %in% c("TICAM2", "TMED7")], sep="_")
# fcs$gene2[fcs$gene %in% c("ICOS")] <- paste("3", fcs$gene[fcs$gene %in% c("ICOS")], sep="_")
# fcs$gene2[fcs$gene %in% c("CTLA4")] <- paste("4", fcs$gene[fcs$gene %in% c("CTLA4")], sep="_")

fcs <- fcs[!fcs$gene == "TMED7-TICAM2",]

tmp <- fcs[fcs$gene %in% c("ERAP1", "ERAP2", "LNPEP"),]
part1 <- ggplot(tmp[tmp$analysis == 'heat_killed',], aes(x=gene2, y=fc)) + 
  # geom_violin(width=1,fill='lightgray') +
  ylab('log fold change') +
  #geom_boxplot( width=0.2, fill='gold3') + 
  # geom_jitter(width=0.2, size=2, col='gray42') + 
  stat_summary(fun.data='mean_cl_boot', geom='pointrange', col='red', size=1) +
  # geom_jitter(data=tmp[tmp$analysis == 'live',],width=0.1, size=3, col='darkorchid2', aes(x=as.numeric(as.factor(gene2))+0.2)) + 
  stat_summary(data=tmp[tmp$analysis == 'live',],fun.data='mean_cl_boot', geom='pointrange', col='darkorchid3', size=1, aes(x=as.numeric(as.factor(gene2))+0.2)) +
  
  xlab('') + ggtitle("rs2549794") + 
  theme_classic() + geom_hline(yintercept = 0, col='black') + theme(text = element_text(size=14))
part1
tmp <- fcs[fcs$gene %in% c("TMED7"),] #, "TICAM2"
part2 <- ggplot(tmp, aes(x=gene2, y=fc)) + 
  # geom_violin(width=1,fill='lightgray') +
  ylab('log fold change') +
  #geom_boxplot( width=0.2, fill='gold3') + 
  # geom_jitter(width=0.2, size=2, col='gray42') + 
  stat_summary(fun.data='mean_cl_boot', geom='pointrange', col='red', size=1) +
  # geom_jitter(data=tmp[tmp$analysis == 'live',],width=0.1, size=3, col='darkorchid2', aes(x=as.numeric(as.factor(gene2))+0.2)) + 
  stat_summary(data=tmp[tmp$analysis == 'live',],fun.data='mean_cl_boot', geom='pointrange', col='darkorchid3', size=1, aes(x=as.numeric(as.factor(gene2))+0.2)) +
  xlab('') + ggtitle("rs17473484") + 
  theme_classic() + geom_hline(yintercept = 0, col='black') + theme(text = element_text(size=14))
part2
tmp <- fcs[fcs$gene %in% c("ICOS"),] #, "CTLA4"
part3 <- ggplot(tmp, aes(x=gene2, y=fc)) + 
  # geom_violin(width=1,fill='lightgray') +
  ylab('log fold change') +
  #geom_boxplot( width=0.2, fill='gold3') + 
  # geom_jitter(width=0.2, size=2, col='gray42') + 
  stat_summary(fun.data='mean_cl_boot', geom='pointrange', col='red', size=1) +
  # geom_jitter(data=tmp[tmp$analysis == 'live',],width=0.1, size=3, col='darkorchid2', aes(x=as.numeric(as.factor(gene2))+0.2)) + 
  stat_summary(data=tmp[tmp$analysis == 'live',],fun.data='mean_cl_boot', geom='pointrange', col='darkorchid3', size=1, aes(x=as.numeric(as.factor(gene2))+0.2)) +
  xlab('') + ggtitle("rs11571319") + 
  theme_classic() + geom_hline(yintercept = 0, col='black') + theme(text = element_text(size=14))
part3

library(patchwork)
fig3a_r2r <- part1 + part2 + part3 + plot_layout(ncol = 3, widths = c(3,2,2, 1)); fig3a_r2r

sifig <- l_vs_hk_1 + fig3a_r2r + plot_layout(ncol = 2, widths = c(1,2)) + plot_annotation(tag_levels = 'A') + theme(axis.text = element_text(color='black')); sifig

plot_pestis <- part1 + part2 + part3 + plot_spacer() + plot_layout(ncol = 4, widths = c(3,2,2, 1)); plot_pestis

fig3a <- ggplot(fcs, aes(x=gene2, y=fc)) + 
  geom_violin(width=1,fill='lightgray') +
  ylab('log2( cpm_infected / cpm_null)') +
  geom_violin(data=fcs[fcs$gene %in% c("CTLA4", "ICOS"),], width=0.7,fill='lightgray') +
  geom_violin(data=fcs[fcs$gene %in% c("TICAM2"),], width=0.7,fill='lightgray') +
  #geom_boxplot( width=0.2, fill='gold3') + 
  geom_jitter(width=0.2, size=2, col='gray42') + 
  geom_vline(xintercept = 3.6, linetype='dashed', col='gray75') +
  geom_vline(xintercept = 5.5, linetype='dashed', col='gray75') +
  geom_vline(xintercept = 6.5, linetype='dashed', col='gray75') +
  # stat_summary(fun='mean', geom='crossbar', width=0.5, size=1, col='red') +
  stat_summary(fun.data='mean_cl_boot', geom='pointrange', col='red', size=1) +
  xlab('') + ggtitle("rs") + 
  theme_classic() + geom_hline(yintercept = 0, col='black') + theme(text = element_text(size=14))
fig3a

save(fig3a, file="./fig3a.rdata") # fig3B.rdata
#save(fig3a_r2r, file="~/fig3a_r2r.rdata")


rm(list=ls(pattern="plot."), part1, part2, part3, fig3c, genes, names, v, individual, pop, sig, stimulated, res, res_lmfit, fig3a, exp_raw, cpm, tmp, exp)
fcs -> fcs_pestis; rm(fcs)

##########################
## Nedelec: 
##########################
load("~/nedelec_results.RData")

fcs_nedelec$gene[fcs_nedelec$gene == "ENSG00000164307"] <- "ERAP1"
fcs_nedelec$gene[fcs_nedelec$gene == "ENSG00000164308"] <- "ERAP2"
fcs_nedelec$gene[fcs_nedelec$gene == "ENSG00000113441"] <- "LNPEP"
fcs_nedelec$gene[fcs_nedelec$gene == "ENSG00000243414"] <- "TICAM2"
fcs_nedelec$gene[fcs_nedelec$gene == "ENSG00000134970"] <- "TMED7"
fcs_nedelec$gene[fcs_nedelec$gene == "ENSG00000131196"] <- "NFATC1"


fcs_nedelec$gene[fcs_nedelec$gene == "ENSG00000172215"] <- "CXCR6"
fcs_nedelec$gene[fcs_nedelec$gene == "ENSG00000163820"] <- "FYCO1"
fcs_nedelec$gene[fcs_nedelec$gene == "ENSG00000163818"] <- "LZTFL1"
fcs_nedelec$gene[fcs_nedelec$gene == "ENSG00000173578"] <- "XCR1"


tmp <- fcs_nedelec[fcs_nedelec$gene %in% c("ERAP1", "ERAP2", "LNPEP"),]
part1 <- ggplot(tmp, aes(x=gene, y=fc, color=stimulation, fill=stimulation)) + 
  geom_violin(fill='lightgray',position=position_dodge(0.9), size=1.25) +
  ylab('log fold change') +
  scale_color_brewer(palette = "Paired") +
  #geom_boxplot( width=0.2, fill='gold3') + 
#  geom_jitter(size=2, col='gray78', position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9)) + 
  stat_summary(fun.data='mean_cl_boot', geom='pointrange', size=1,position=position_dodge(0.9)) +
  xlab('') + ggtitle("rs2549794") + 
  theme_classic() + geom_hline(yintercept = 0, col='black') + theme(text = element_text(size=14))
  #, legend.position='none'
part1

tmp <- fcs_nedelec[fcs_nedelec$gene %in% c("TICAM2", "TMED7"),]
part2 <- ggplot(tmp, aes(x=gene, y=fc, color=stimulation, fill=stimulation)) + 
  geom_violin(fill='lightgray',position=position_dodge(0.9), size=1.25) +
  ylab('log fold change') +
  scale_color_brewer(palette = "Paired") +
  #geom_boxplot( width=0.2, fill='gold3') + 
#  geom_jitter(size=2, col='gray78', position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9)) + 
  stat_summary(fun.data='mean_cl_boot', geom='pointrange', size=1,position=position_dodge(0.9)) +
  xlab('') + ggtitle("rs2549794") + 
  theme_classic() + geom_hline(yintercept = 0, col='black') + theme(text = element_text(size=14))
#, legend.position='none
part2
tmp <- fcs_nedelec[fcs_nedelec$gene %in% c("NFATC1"),]
part3 <- ggplot(tmp, aes(x=gene, y=fc, color=stimulation, fill=stimulation)) + 
  geom_violin(fill='lightgray',position=position_dodge(0.9), size=1.25) +
  ylab('log fold change') +
  scale_color_brewer(palette = "Paired") +
  #geom_boxplot( width=0.2, fill='gold3') + 
#  geom_jitter(size=2, col='gray78', position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9)) + 
  stat_summary(fun.data='mean_cl_boot', geom='pointrange', size=1,position=position_dodge(0.9)) +
  xlab('') + ggtitle("rs1052025") + 
  theme_classic() + geom_hline(yintercept = 0, col='black') + theme(text = element_text(size=14))
#, legend.position='none
part3
 
plot_nedelec <- part1 + part2 + plot_spacer() +  part3 + plot_layout(guides = "collect", ncol = 4, widths = c(3,2,1, 1))
##########################

##########################
## Quach: 
##########################
load("~/quach_results.RData")

fcs_quach$gene[fcs_quach$gene == "ENSG00000164307"] <- "ERAP1"
fcs_quach$gene[fcs_quach$gene == "ENSG00000164308"] <- "ERAP2"
fcs_quach$gene[fcs_quach$gene == "ENSG00000113441"] <- "LNPEP"
fcs_quach$gene[fcs_quach$gene == "ENSG00000243414"] <- "TICAM2"
fcs_quach$gene[fcs_quach$gene == "ENSG00000134970"] <- "TMED7"
fcs_quach$gene[fcs_quach$gene == "ENSG00000131196"] <- "NFATC1"
fcs_quach$gene[fcs_quach$gene == "ENSG00000163600"] <- "ICOS"

tmp <- fcs_quach[fcs_quach$gene %in% c("ERAP1", "ERAP2", "LNPEP"),]
part1 <- ggplot(tmp, aes(x=gene, y=fc, color=stimulation, fill=stimulation)) + 
  geom_violin(fill='lightgray',position=position_dodge(0.9), size=1.25) +
  ylab('log fold change') +
  scale_color_brewer(palette="Dark2") +
  #geom_boxplot( width=0.2, fill='gold3') + 
#  geom_jitter(size=2, col='gray78', position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9)) + 
  stat_summary(fun.data='mean_cl_boot', geom='pointrange', size=1,position=position_dodge(0.9)) +
  xlab('') + ggtitle("rs2549794") + 
  theme_classic() + geom_hline(yintercept = 0, col='black') + theme(text = element_text(size=14))
#, legend.position='none'
part1
tmp <- fcs_quach[fcs_quach$gene %in% c("TICAM2", "TMED7"),]
part2 <- ggplot(tmp, aes(x=gene, y=fc, color=stimulation, fill=stimulation)) + 
  geom_violin(fill='lightgray',position=position_dodge(0.9), size=1.25) +
  ylab('log fold change') +
  scale_color_brewer(palette = "Dark2") +
  #geom_boxplot( width=0.2, fill='gold3') + 
#  geom_jitter(size=2, col='gray78', position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9)) + 
  stat_summary(fun.data='mean_cl_boot', geom='pointrange', size=1,position=position_dodge(0.9)) +
  xlab('') + ggtitle("rs2549794") + 
  theme_classic() + geom_hline(yintercept = 0, col='black') + theme(text = element_text(size=14))
#, legend.position='none
part2
tmp <- fcs_quach[fcs_quach$gene %in% c("ICOS"),]
part3 <- ggplot(tmp, aes(x=gene, y=fc, color=stimulation, fill=stimulation)) + 
  geom_violin(fill='lightgray',position=position_dodge(0.9), size=1.25) +
  ylab('log fold change') +
  scale_color_brewer(palette = "Dark2") +
  #geom_boxplot( width=0.2, fill='gold3') + 
#  geom_jitter(size=2, col='gray78', position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9)) + 
  stat_summary(fun.data='mean_cl_boot', geom='pointrange', size=1,position=position_dodge(0.9)) +
  xlab('') + ggtitle("rs1052025") + 
  theme_classic() + geom_hline(yintercept = 0, col='black') + theme(text = element_text(size=14))
#, legend.position='none
part3
tmp <- fcs_quach[fcs_quach$gene %in% c("NFATC1"),]
part4 <- ggplot(tmp, aes(x=gene, y=fc, color=stimulation, fill=stimulation)) + 
  geom_violin(fill='lightgray',position=position_dodge(0.9), size=1.25) +
  ylab('log fold change') +
  scale_color_brewer(palette = "Dark2") +
  #geom_boxplot( width=0.2, fill='gold3') + 
#  geom_jitter(size=2, col='gray78', position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9)) + 
  stat_summary(fun.data='mean_cl_boot', geom='pointrange', size=1,position=position_dodge(0.9)) +
  xlab('') + ggtitle("rs1052025") + 
  theme_classic() + geom_hline(yintercept = 0, col='black') + theme(text = element_text(size=14))
#, legend.position='none
part4

plot_quach <- part1 + part2 +  part3 + part4 + plot_layout(guides = "collect", ncol = 4, widths = c(3,2,1, 1))
##########################

fig3a_r2r / plot_nedelec / plot_quach + plot_layout(guides='collect')

plot_nedelec / plot_quach + plot_layout(guides='collect')
fig3a_r2r



############################

load("./DATA/n33_pestis_infection.RData")
res_lmfit -> res_pestis; rm(res_lmfit, v, individual, pop, stimulated, cpm, exp, names)
genes -> genes_pestis; rm(genes)

res_pestis$gene <- row.names(res_pestis); res_pestis$type <- '00_pestis'
res_LPS$gene <- row.names(res_LPS); res_LPS$type <- "03_LPS"
res_IAV$gene <- row.names(res_IAV); res_IAV$type <- "05_IAV"
res_PAM3CSK4$gene <- row.names(res_PAM3CSK4); res_PAM3CSK4$type <- '04_PAM3CSK4'
res_R848$gene <- row.names(res_R848); res_R848$type <- "06_R848"
res_salmonella$type <- '01_salmonella'
res_listeria$type <- '02_listeria'

all_res <- NULL; for (f in ls(pattern="res_")) {all_res <- rbind(all_res, get(f))}

all_res$gene[all_res$gene == "ENSG00000164307"] <- "ERAP1"
all_res$gene[all_res$gene == "ENSG00000164308"] <- "ERAP2"
all_res$gene[all_res$gene == "ENSG00000113441"] <- "LNPEP"
all_res$gene[all_res$gene == "ENSG00000243414"] <- "TICAM2"
all_res$gene[all_res$gene == "ENSG00000134970"] <- "TMED7"
all_res$gene[all_res$gene == "ENSG00000131196"] <- "NFATC1"
all_res$gene[all_res$gene == "ENSG00000163600"] <- "ICOS"
all_res$gene[all_res$gene == "ENSG00000163599"] <- "CTLA4"

all_res <- all_res[all_res$gene %in% c("ERAP1", "ERAP2", "LNPEP", "TICAM2", "TMED7", "CTLA4", "NFATC1", "ICOS"),]

all_res$gene2[all_res$gene %in% c("ERAP1", "ERAP2", "LNPEP")] <- paste("1", all_res$gene[all_res$gene %in% c("ERAP1", "ERAP2", "LNPEP")], sep="_")
all_res$gene2[all_res$gene %in% c("TICAM2", "TMED7")] <- paste("2", all_res$gene[all_res$gene %in% c("TICAM2", "TMED7")], sep="_")
all_res$gene2[all_res$gene %in% c("ICOS","CTLA4")] <- paste("3", all_res$gene[all_res$gene %in% c("ICOS","CTLA4")], sep="_")
all_res$gene2[all_res$gene %in% c("NFATC1")] <- paste("4", all_res$gene[all_res$gene %in% c("NFATC1")], sep="_")

panel_stim_by_study <- ggplot(all_res, aes(y=factor(gene2), x=factor(type))) +
  geom_point(aes(color=logFC, size= (P.Value)), shape=16) +
  scale_color_gradient2(low='purple4', mid='snow2', high='tomato2') + 
  scale_size(range=c(20,5)) + ylab("gene") + xlab("stimulation") + 
  geom_hline(yintercept = c(0.5,1.5,3.5,5.5), size = .2) +
  scale_y_discrete(limits=rev) + scale_x_discrete(position = "top")  +
  geom_text(aes(label=signif(logFC, digits = 2)), nudge_y = 0.1) +
  geom_text(aes(label=paste0("p = ",signif(P.Value, digits = 2))), nudge_y = -0.1) +
  theme_classic() + theme(legend.position = "bottom", 
                          axis.text=element_text(color='black'), text=element_text(size=16)
                          )

save(panel_stim_by_study, file = "~/stim_by_study.RData")

load("~/quach_results2.RData")
load("~/nedelec_results.RData")
load("~/qtl_pestis.rdata")
eQTL_quach <- eQTL_quach[,c(4:8)]
eQTL_nedelec <- eQTL_nedelec[,c(6,5,7,4,1)]
colnames(eQTL_nedelec) <- colnames(eQTL_quach) <- c("gene", "condition", "snp", "P.Value", "Beta")
eQTL_quach$Beta[eQTL_quach$snp == "rs2549794"] <- -eQTL_quach$Beta[eQTL_quach$snp == "rs2549794"]
eQTL_nedelec$Beta[eQTL_nedelec$snp == "rs17473484"] <- -eQTL_nedelec$Beta[eQTL_nedelec$snp == "rs17473484"]

eQTL <- rbind(eQTL_quach[!eQTL_quach$condition %in% c("NS", "NI"),], eQTL_nedelec[!eQTL_nedelec$condition %in% c("NI","all"),], eQTL)

eQTL <- eQTL[eQTL$snp %in% c("rs2549794", "rs17473484"),]

eQTL$type[eQTL$condition == 'NI'] <- '00_non-infected'
eQTL$type[eQTL$condition == 'pestis'] <- '00_pestis'
eQTL$type[eQTL$condition == 'S'] <- '01_salmonella'
eQTL$type[eQTL$condition == 'L'] <- '02_listeria'
eQTL$type[eQTL$condition == 'LPS'] <- '03_LPS'
eQTL$type[eQTL$condition == 'PAM3CSK4'] <- '04_PAM3CSK4'
eQTL$type[eQTL$condition == 'IAV'] <- '05_IAV'
eQTL$type[eQTL$condition == 'R848'] <- '06_R848'

eQTL$gene[eQTL$gene == "ENSG00000164307"] <- "ERAP1"
eQTL$gene[eQTL$gene == "ENSG00000164308"] <- "ERAP2"
eQTL$gene[eQTL$gene == "ENSG00000113441"] <- "LNPEP"
eQTL$gene[eQTL$gene == "ENSG00000243414"] <- "TICAM2"
eQTL$gene[eQTL$gene == "ENSG00000134970"] <- "TMED7"
eQTL$gene[eQTL$gene == "ENSG00000131196"] <- "NFATC1"
eQTL$gene[eQTL$gene == "ENSG00000163600"] <- "ICOS"
eQTL$gene[eQTL$gene == "ENSG00000163599"] <- "CTLA4"

fcs_nedelec$gene[fcs_nedelec$gene == "ENSG00000172215"] <- "CXCR6"
fcs_nedelec$gene[fcs_nedelec$gene == "ENSG00000163820"] <- "FYCO1"
fcs_nedelec$gene[fcs_nedelec$gene == "ENSG00000163818"] <- "LZTFL1"
fcs_nedelec$gene[fcs_nedelec$gene == "ENSG00000173578"] <- "XCR1"

eQTL <- eQTL[eQTL$gene %in% c("ERAP2", "TICAM2"),]

panel_eqtl_by_study <- ggplot(eQTL, aes(y=paste(snp,factor(gene),sep=":"), x=factor(type))) +
  geom_point(aes(color=Beta, size= (P.Value)), shape=16) +
  scale_color_gradient2(low='purple4', mid='snow2', high='tomato2') + 
  scale_size(range=c(20,5)) + ylab("gene") + xlab("stimulation") + 
  geom_hline(yintercept = c(0.5,1.5,3.5,5.5), size = .2) +
  scale_y_discrete(limits=rev) + scale_x_discrete(position = "top")  +
  geom_text(aes(label=signif(Beta, digits = 2)), nudge_y = 0.1) +
  geom_text(aes(label=paste0("p = ",signif(P.Value, digits = 2))), nudge_y = -0.1) +
  theme_classic() + theme(legend.position = "bottom", 
                          axis.text=element_text(color='black'), text=element_text(size=16)
  ); panel_eqtl_by_study

save(panel_eqtl_by_study, file = "~/eqtl_by_study.RData")


##### 
## QTL BOXPLOTS

load("~/pestis_eqtl_s8.rdata")
eQTL_data$stimulated[eQTL_data$stimulated == F] <- "0_control"
eQTL_data$stimulated[eQTL_data$stimulated == T] <- "0_pestisl"

load("~/quach_results2.RData")
load("~/nedelec_results.RData"); rm(target_genes_nedelec)
colnames(nedelec_eqtl_s8) <- colnames(quach_eqtl_s8) <- colnames(eQTL_data)

quach_eqtl_s8$gt[quach_eqtl_s8$gene == "ERAP2"] <- 2- quach_eqtl_s8$gt[quach_eqtl_s8$gene == "ERAP2"]

eQTL_data <- rbind(eQTL_data, nedelec_eqtl_s8, quach_eqtl_s8)

eQTL_data$stimulated[eQTL_data$stimulated == 'S'] <- '01_salmonella'
eQTL_data$stimulated[eQTL_data$stimulated == 'L'] <- '02_listeria'
eQTL_data$stimulated[eQTL_data$stimulated == 'LPS'] <- '03_LPS'
eQTL_data$stimulated[eQTL_data$stimulated == 'PAM3CSK4'] <- '04_PAM3CSK4'
eQTL_data$stimulated[eQTL_data$stimulated == 'IAV'] <- '05_IAV'
eQTL_data$stimulated[eQTL_data$stimulated == 'R848'] <- '06_R848'

eQTL_data$gene[eQTL_data$gene == "ENSG00000164308"] <- "ERAP2"
eQTL_data$gene[eQTL_data$gene == "ENSG00000243414"] <- "TICAM2"

# noms <- c(
#   'FALSE' = 'Not-Infected', 
#   'TRUE' = 'Infected'
# )

eQTL_data <- eQTL_data[!is.na(eQTL_data$gt),]

plot <- ggplot(data=eQTL_data, aes(y=expression, x=as.factor(gt), fill=stimulated)) + 
  facet_grid(gene ~ stimulated) + #, labeller = as_labeller(noms)
  geom_boxplot() + 
  #geom_boxplot(width=0.25) + 
  # geom_jitter(width=0.2, size=2, col='gray42') + 
  # stat_summary(fun.data='mean_cl_boot', geom='pointrange', col='red', size=1) +
  # scale_fill_manual(values=c("lightgoldenrodyellow", "lightgoldenrod2"),
  #                   name="", breaks=c("FALSE", "TRUE"),
  #                   labels=c("Not-Infected", "Infected")) + 
  theme_classic() +
  #scale_fill_manual(values=c("#999999", "#E69F00"), name="Population",breaks=c("FALSE", "TRUE"),labels=c("African", "European")) + 
  ylab(label = "normalized CPM") + xlab(label=paste("alternate alleles",sep="")) +
  theme(legend.position = "none") + 
  theme(strip.background = element_rect(color='NA', fill='NA', size=1.5)); plot


