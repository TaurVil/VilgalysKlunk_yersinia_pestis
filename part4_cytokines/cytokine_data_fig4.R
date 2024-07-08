setwd("C:/Users/taurv/Documents/GitHub/VilgalysKlunk_yersinia_pestis/part4_cytokines")

## Read in data
# each row contains the cytokine level for pestis-stimulated samples versus the null control samples at the same time point. we'll calculate the difference between them below, and call it `value`
###########################
data <- read.delim("./DATA/cytokines.txt")
data$value <- data$Pestis - data$Null
# Effects become more obvious over time. We therefore focus on the 24h time point rather than the 5h timepoint. 
data <- subset(data, data$time == '24h') 
# get the list of cytokines analyzed
cks <- unique(data$cytokine)
###########################

## Reshape into square format for PCA
## PCA with all samples and "rescale" option
###########################
library(ggplot2); library(reshape2)
data$selection <- gsub("2_deleterious", "Nonexpressed", data$selection)
data$id <- paste(data$sample, data$time, data$status, data$selection, sep="_")

d2 <- reshape(data[,c(11,10,7)], idvar='id', timevar='cytokine', direction='wide')
colSums(is.na(d2))
d2$`value.CXCL8/IL8`[is.na(d2$`value.CXCL8/IL8`)] <- mean(d2$`value.CXCL8/IL8`[1:25], na.rm=T)

pca <- prcomp(d2[,2:11], scale. = T)
library(plyr); strsplit(d2$id, "_") -> pc_names; ldply(pc_names) -> pc_names

tmp <- cbind(pca$x[,1], pca$x[,2], pca$x[,3], pca$x[,4], pc_names)
colnames(tmp) <- c("PC1", "PC2", "PC3", "PC4","SAMPLE", "TIME", "GT", "STATUS")
p2 <- ggplot(tmp[tmp$TIME == "24h",], aes(x=PC1, y=PC2, color = GT)) + geom_point(size=6) + theme_classic() + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
p1 <- ggplot(tmp, aes(x=GT, y=PC1, fill=GT)) + geom_boxplot() + coord_flip() +  theme_classic() + scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + xlab("genotype"); p1
# scale_x_discrete(limits = rev(levels(as.factor(tmp$GT)))) +

library(patchwork); part1 <- p2 / p1 + plot_layout(heights = c(4,1)) ; part1

summary(lm(data = tmp[tmp$TIME == "24h",], PC1 ~ as.numeric(as.character(GT))))
summary(lm(data = tmp[tmp$TIME == "24h",], PC2 ~ as.numeric(as.character(GT))))
summary(lm(data = tmp[tmp$TIME == "24h",], PC3 ~ as.numeric(as.character(GT))))
summary(lm(data = tmp[tmp$TIME == "24h",], PC4 ~ as.numeric(as.character(GT))))

###########################

## For each cytokine, measure cytokine level as a function of genotype and plot
## Normalize data for each cytokine using qqnorm
###########################
data$GT[data$genotype == "Allele 2"] <- 0
data$GT[data$genotype == "Allele 1"] <- 2
data$GT[data$genotype == "HE"] <- 1

# replace "value" with "Null" to get the p-values for unstimulated data in line 53
for (ck in cks) {
  tmp <- data[data$cytokine == ck,]
  print(ck)
  tmp2 <- tmp[tmp$time == '24h',]
  tmp2$value <- qqnorm(tmp2$value, plot.it = F)$x
  
  print(summary(lm(tmp2$value ~ as.numeric(as.character(tmp2$GT))))$coefficients)
  print(summary(lm(tmp2$value ~ as.numeric(as.character(tmp2$GT))))$r.squared)
  #print(summary(lm(tmp2$value ~ tmp2$selection))$coefficients)
}
###########################

## Plot the 4 significant cytokines
###########################
data$gt <- gsub("Allele 1", "2_C/C", data$genotype)
data$gt <- gsub("Allele 2", "0_T/T", data$gt)
data$gt <- gsub("HE", "1_T/C", data$gt)

tmp <- data[data$cytokine == "G-CSF",]; tmp$value <- qqnorm(tmp$value, plot.it = F)$x
p3 <- ggplot(data=tmp, aes(y=value, x=as.factor(gt), fill=as.factor(gt))) + 
  geom_boxplot(outlier.shape = NA) + ggtitle("G-CSF") + geom_jitter(width=0.2,col='gray45', size=2, show.legend = F) + xlab("rs249794 genotype") +
  theme_classic() + scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
tmp <- data[data$cytokine == "IL1B",]; tmp$value <- qqnorm(tmp$value, plot.it = F)$x
p4 <- ggplot(data=tmp, aes(y=value, x=as.factor(gt), fill=as.factor(gt))) + 
  geom_boxplot(outlier.shape = NA) + ggtitle("IL1B") + geom_jitter(width=0.2,col='gray45', size=2, show.legend = F) + xlab("rs249794 genotype") +
  theme_classic() + scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
tmp <- data[data$cytokine == "IL10",]; tmp$value <- qqnorm(tmp$value, plot.it = F)$x
p5 <- ggplot(data=tmp, aes(y=value, x=as.factor(gt), fill=as.factor(gt))) + 
  geom_boxplot(outlier.shape = NA) + ggtitle("IL10") + geom_jitter(width=0.2,col='gray45', size=2, show.legend = F) + xlab("rs249794 genotype") +
  theme_classic() + scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
tmp <- data[data$cytokine == "CCL3",]; tmp$value <- qqnorm(tmp$value, plot.it = F)$x
p6 <- ggplot(data=tmp, aes(y=value, x=as.factor(gt), fill=as.factor(gt))) + 
  geom_boxplot(outlier.shape = NA) + ggtitle("CCL3") + geom_jitter(width=0.2,col='gray45', size=2, show.legend = F) + xlab("rs249794 genotype") +
  theme_classic() + scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))

part2 <- p3 + p4 + p5 + p6 
part2 <- part2 + plot_annotation(tag_levels = 'A') + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
part2
###########################


## Merge elements of figure 4 in Illustrator (no idea why this didn't work smoothly)
###########################

cfu$gt <- gsub("Allele 1", "2_C/C", cfu$allele)
cfu$gt <- gsub("Allele 2", "0_T/T", cfu$gt)
cfu$gt <- gsub("Allèle 1", "2_C/C", cfu$gt)
cfu$gt <- gsub("HE", "1_T/C", cfu$gt)

pcfu <- ggplot(data=cfu, aes(y=perc_kill, x=as.factor(gt), fill=as.factor(gt))) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width=0.2,col='gray45', size=4, show.legend = F) + xlab("rs249794 genotype") + ylab("proportion of colony forming units remaining after 24h") + ggtitle("remaining colony forming units") +
  theme_classic() + scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) + theme(legend.position = 'none'); pcfu

part2 | pcfu1 + plot_annotation(tag_levels = 'A') + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')

( p2 / p1 ) + part2 + plot_annotation(tag_levels = 'A') + plot_layout(guides = 'collect', nrow = 2, ncol=3)
