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
save.image("./DATA/n33_pestis_infection.RData")

library(data.table); library(DESeq2); library(ggplot2); library(limma)
load("./DATA/n33_pestis_infection.RData")
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
  "CTLA4", "ICOS", # chr2
  "ERAP1", "ERAP2", "LNPEP", # chr5
  "TICAM2", "TMED7", "TMED7-TICAM2", # chr5
  "NFATC1" # chr18    
  )
tmp <- subset(res, res$Symbol %in% sig)
tmp[order(tmp$logFC),]

## expression of the 8 target genes we analyze
fcs <- NULL
for (tmp_num in 1:length(tmp$Symbol)) {
  tmp_gene_name <- tmp$Symbol[tmp_num]
  print("starting:"); print(tmp_gene_name)
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

library(patchwork)
#8.5 by 4.5 inches
(plot.TICAM2 + plot.TMED7 + `plot.TMED7-TICAM2` + plot.CTLA4 +
    plot.ICOS + plot.LNPEP + plot.ERAP1 + plot.ERAP2) + plot_layout(ncol=4)

# fcs$fc[fcs$YPcpm == 0] <- -6
# fcs$fc[fcs$NIcpm == 0] <- 5
# fcs$fc[fcs$NIcpm == 0 & fcs$YPcpm == 0] <- 0

fcs$gene2[fcs$gene %in% c("ERAP1", "ERAP2", "LNPEP")] <- paste("1", fcs$gene[fcs$gene %in% c("ERAP1", "ERAP2", "LNPEP")], sep="_")
fcs$gene2[fcs$gene %in% c("TICAM2", "TMED7")] <- paste("2", fcs$gene[fcs$gene %in% c("TICAM2", "TMED7")], sep="_")
fcs$gene2[fcs$gene %in% c("ICOS")] <- paste("3", fcs$gene[fcs$gene %in% c("ICOS")], sep="_")
fcs$gene2[fcs$gene %in% c("CTLA4")] <- paste("4", fcs$gene[fcs$gene %in% c("CTLA4")], sep="_")

fcs <- fcs[!fcs$gene == "TMED7-TICAM2",]

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

# save(fig3a, file="~/Barreiro/pestis_aDNA/fig3B.rdata")
