## R2R: recombination rates

dist = 50000 #window size around each snp, we'll divide this by 2 to get the start/end

## Get candidate loci 
library(data.table); library(ggplot2); library(parallel)
load("./candidate_loci.RData")
## we're going to test the 245 loci which are highly differentiated in London against the ~3k tested variants and ~700 neutral variants

test <- A_all_loci[,1:10]
neut <- info_neut[,1:10]
test$set <- "test"; test$set[test$site %in% B_L13_95th$site] <- "candidate"; neut$set <- "neutral"
for_rcr <- rbind(test, neut)

## get window around each candidate loci 
for_rcr$start <- for_rcr$pos - dist/2; for_rcr$end <- for_rcr$pos + dist/2

table(for_rcr$set)


## get mean recombination rate near each site
## hapmap rate is cM per Mb
with_rcrs <- NULL

calc_mean_rcr <- function(window, tmp) {
  end <- tmp$end[window]
  st <- tmp$start[window]
  k2 <- rcr22$start <= end & rcr22$end > st
  
  rcr22[c(k2==T),] -> tmp2
  tmp2$start[tmp2$start < tmp$start[window]] <- tmp$start[window]
  tmp2$end[tmp2$end > tmp$end[window]] <- tmp$end[window]
  tmp2$l <- tmp2$end - tmp2$start
  return(sum(tmp2$rcr*tmp2$l)/sum(tmp2$l))
}

for (k in 1:22) {
  tmp_chrom <- paste0("chr",k)
  name=paste0("~/my_genomes/hg37/HapMap_Recombination/genetic_map_GRCh37_chr",k,".txt")
  rcr22 <- fread(name)
  colnames(rcr22)[2] <- 'start'
  rcr22$end <- c(rcr22$start[-1], max(rcr22$start)+1)
  colnames(rcr22)[3] <- 'rcr'
  
  tmp_sites <- for_rcr[for_rcr$chr == tmp_chrom,]
  tmp_sites$hapmap_rcr <- do.call("c",mclapply(1:nrow(tmp_sites), FUN = calc_mean_rcr, tmp=tmp_sites))
  
  with_rcrs <- rbind(with_rcrs, tmp_sites); print(k)
  rm(name, rcr22, tmp_sites, tmp_chrom)
}; rm(k, calc_mean_rcr)
with_rcrs -> for_rcr

calc_mean_rcr <- function(window, map,tmp) {
  end <- tmp$end[window]
  st <- tmp$start[window]
  k2 <- rcr22$start <= end & rcr22$end > st
  
  rcr22[c(k2==T),] -> tmp2
  tmp2$start[tmp2$start < tmp$start[window]] <- tmp$start[window]
  tmp2$end[tmp2$end > tmp$end[window]] <- tmp$end[window]
  tmp2$l <- tmp2$end - tmp2$start
  return(sum(tmp2[[map]]*tmp2$l)/sum(tmp2$l))
}
with_rcrs <- NULL; for (k in 1:22) {
  tmp_chrom <- paste0("chr",k)
  name=paste0("E:/Backup_Dec2022/my_genomes/hg37/maps_b37/maps_chr.",k)
  rcr22 <- fread(name)
  colnames(rcr22)[1] <- 'end'
  rcr22$start <- c(0, rcr22$end[-nrow(rcr22)])
  rcr22$l <- rcr22$end - rcr22$start
  d <- rcr22[,2:7]
  d2 <- rbind(0,d[-nrow(d),], fill=T); d2[1,] <- 0; d2[,-1] -> d2
  diff <- d-d2
  rcr22[,2:7] <- diff/rcr22$l; rm(d, d2, diff)
  
  tmp_sites <- for_rcr[for_rcr$chr == tmp_chrom,]
  tmp_sites$YRI <- do.call("c",mclapply(1:nrow(tmp_sites), FUN = calc_mean_rcr, map="YRI_LD", tmp=tmp_sites))
  tmp_sites$deCODE_Icelandic <- do.call("c",mclapply(1:nrow(tmp_sites), FUN = calc_mean_rcr, map="deCODE", tmp=tmp_sites))
  tmp_sites$AA <- do.call("c",mclapply(1:nrow(tmp_sites), FUN = calc_mean_rcr, map="AA_Map", tmp=tmp_sites))
  tmp_sites$CEU <- do.call("c",mclapply(1:nrow(tmp_sites), FUN = calc_mean_rcr, map="CEU_LD", tmp=tmp_sites))
  
  
  with_rcrs <- rbind(with_rcrs, tmp_sites); print(k)
  rm(name, rcr22, tmp_sites, tmp_chrom)
}; rm(k, calc_mean_rcr)


model <- lm(with_rcrs$CEU ~ with_rcrs$set)
anova(model)


rm(splot_enrichment, test, neut, panel_D_enrichment, panel_enrichment, panel_manhattan, model, for_rcr, don)

with_rcrs$set <- gsub("candidate", "03_candidate",with_rcrs$set)
with_rcrs$set <- gsub("test", "02_test",with_rcrs$set)
with_rcrs$set <- gsub("neutral", "01_neutral",with_rcrs$set)

save.image("./sites_with_rcr.RData")



load("./sites_with_rcr.RData")

# with_rcrs$deCODE_Icelandic <- log10(with_rcrs$deCODE_Icelandic)

library(ggplot2)
plot_recombination <- ggplot(data=with_rcrs, aes(x=factor(set), y=log10(deCODE_Icelandic), fill=set)) +
  geom_violin() + 
  scale_fill_brewer(palette="Blues") +
  geom_boxplot(width=0.1, fill='white', outlier.shape = NA) + 
  theme_classic() + 
  xlab("set of sites") + ylab("recombination rate")
plot_recombination

# save.image("./sites_with_rcr.RData")

t.test(with_rcrs$deCODE_Icelandic[with_rcrs$set == "02_test"], with_rcrs$deCODE_Icelandic[with_rcrs$set == "03_candidate"])
t.test(with_rcrs$deCODE_Icelandic[with_rcrs$set == "01_neutral"], with_rcrs$deCODE_Icelandic[!with_rcrs$set == "01_neutral"])

t.test(with_rcrs$deCODE_Icelandic[with_rcrs$set == "03_candidate"], with_rcrs$deCODE_Icelandic[!with_rcrs$set == "01_neutral"])

t.test(log10(with_rcrs$deCODE_Icelandic[with_rcrs$set == "03_candidate"]+1e-7), log10(with_rcrs$deCODE_Icelandic[!with_rcrs$set == "01_neutral"]+1e-7))

t.test(log10(with_rcrs$deCODE_Icelandic[with_rcrs$set == "01_neutral"]+1e-7), log10(with_rcrs$deCODE_Icelandic[!with_rcrs$set == "01_neutral"]+1e-7))

save.image("./sites_with_rcr.RData")
