## R2R: background selection 

load("./sites_with_rcr.RData")

## Get B values## these were hg18
## moved to hg19/37 using http://grch37.ensembl.org/Homo_sapiens/Tools/AssemblyConverter?db=core. From NNCBI36/hg18
library(plyr); library(dplyr); library(data.table); library(parallel)

calc_mean_B <- function(window) {
  subset(b, b$start <= tmp$end[window] & b$end > tmp$start[window]) -> tmp2
  tmp2$start[tmp2$start < tmp$start[window]] <- tmp$start[window]
  tmp2$end[tmp2$end > tmp$end[window]] <- tmp$end[window]
  tmp2$l <- tmp2$end - tmp2$start
  return(sum(tmp2$Bvalue*tmp2$l)/sum(tmp2$l))
}

all_b <- rbind(read.delim("~/my_genomes/hg37/newB1.bed", header=F), read.delim("~/my_genomes/hg37/newB2.bed", header=F), read.delim("~/my_genomes/hg37/newB3.bed", header=F))
colnames(all_b) <- c('chr', 'start', 'end', 'Bvalue')

all2 <- NULL; for (k in 1:22) {
  b <- subset(all_b, all_b$chr == paste0("chr",k)) 
  b$chr <- paste("chr",k,sep="")
  
  tmp <- subset(with_rcrs, with_rcrs$chr == paste("chr",k,sep=""))
  tmp$B <- do.call("c",mclapply(1:nrow(tmp), FUN = calc_mean_B))
  
  all2 <- rbind(all2, tmp); print(k) ; rm(tmp)
  rm(b)
}
all2 -> with_beta; rm(all2, calc_mean_B, k, all_b)

#plot(all$B2 ~ all$B); abline(a=0,b=1,col='red')
rm(empty_windows)

mean(with_beta$B[with_beta$set == "01_neutral"])
mean(with_beta$B[with_beta$set == "02_test"])
mean(with_beta$B[with_beta$set == "03_candidate"])

t.test(with_beta$B[with_rcrs$set == "01_neutral"], with_beta$B[!with_beta$set == "01_neutral"])
t.test(with_beta$B[with_rcrs$set == "01_neutral"], with_beta$B[!with_beta$set == "01_neutral"])$p.value

t.test(with_beta$B[with_rcrs$set == "02_test"], with_beta$B[with_beta$set == "03_candidate"])
mean(with_beta$B[with_beta$set == "03_candidate"])
with_beta$B[with_beta$site %in% D_replicated_in_Denmark$site]
with_beta$site[with_beta$site %in% D_replicated_in_Denmark$site]

save.image("./sites_with_rcr_with_B.RData")


# with beta if using B value files I can't seem to find anymore. 
# pull from original paper (load(../part1_aDNA/sites_with_rcr_with_B.RData))


quantile(with_beta$B[with_beta$set == '01_neutral'], c(0.025, 0.975))
with_beta$B[with_beta$site %in% D_replicated_in_Denmark$site]



with_beta$set <- gsub("candidate", "03_candidate",with_beta$set)
with_beta$set <- gsub("test", "02_test",with_beta$set)
with_beta$set <- gsub("neutral", "01_neutral",with_beta$set)

library(ggplot2)
plot_bs <- ggplot(data=with_beta, aes(x=factor(set), y=B, fill=set)) +
  geom_violin() + 
  scale_fill_brewer(palette="Blues") +
  geom_boxplot(width=0.1, fill='white', outlier.shape = NA) + 
  theme_classic() + 
  xlab("set of sites") + ylab("recombination rate")
plot_bs

library(patchwork)

left_panels <- plot_recombination / plot_bs + plot_layout(guides="collect") & theme(legend.position = 'bottom') + theme(text = element_text(size = 16, color='black'))
left_panels
