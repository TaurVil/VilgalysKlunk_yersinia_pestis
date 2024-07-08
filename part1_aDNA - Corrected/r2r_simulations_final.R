######## 
# Get study design matrix and R packages
######## 
library(data.table)

num=21288093
set.seed(21288093)

read.delim(paste0("./DATA/sims.", num, ".exon.txt"), header=F) -> exon; exon$type <- 'exon'
read.delim(paste0("./DATA/sims.", num, ".immune.txt"), header=F) -> immune; immune$type <- 'immune'
read.delim(paste0("./DATA/sims.", num, ".neutral.txt"), header=F) -> neutral; neutral$type <- 'neutral'
rbind(exon, immune, neutral) -> sites; rm(exon, immune,neutral)
  
sites$alt_af <- rowMeans(sites[,4:5])
colnames(sites)[1:5] <- c('chr', 'pos1', 'pos', 'pre', 'post')
  
########
# shift from alternate frequency to minor allele frequency
########
sites$maf <- apply(cbind(sites$alt_af, 1-sites$alt_af), 1, min)
sites$pre[sites$maf != sites$alt_af] <- 1-sites$pre[sites$maf != sites$alt_af]
sites$post[sites$maf != sites$alt_af] <- 1-sites$post[sites$maf != sites$alt_af]
########

## adjustment for sampling error, 20 individuals so 40 chromosomes
new_sites <- NULL; for (i in 1:3) {
  nindiv_pre=rpois(nrow(sites),50)
  nindiv_post=rpois(nrow(sites),50)
  
  sites -> tmp
  
  tmp$pre <- sapply(1:nrow(sites),function (x) {sum(runif(2*nindiv_pre[x]) < sites$pre[x])/(2*nindiv_pre[x])})
  tmp$post <- sapply(1:nrow(sites),function (x) {sum(runif(2*nindiv_post[x]) < sites$post[x])/(2*nindiv_post[x])})
  tmp$alt_af <- rowMeans(tmp[,4:5])
  tmp$maf <- apply(cbind(tmp$alt_af, 1-tmp$alt_af), 1, min)
  tmp$pre[tmp$maf != tmp$alt_af] <- 1-tmp$pre[tmp$maf != tmp$alt_af]
  tmp$post[tmp$maf != tmp$alt_af] <- 1-tmp$post[tmp$maf != tmp$alt_af]
  
  new_sites <- rbind(tmp, new_sites)
}
print("25 sets of sites sampled")
new_sites -> sites; rm(new_sites, i, tmp)

########
# Calculate Fst
########
sites$ehet_pre <- 2*sites$pre*(1-sites$pre)
sites$ehet_post <- 2*sites$post*(1-sites$post)
sites$ehet <- 2*sites$maf*(1-sites$maf)
sites$fst <- (sites$ehet - (sites$ehet_pre + sites$ehet_post)/2)/sites$ehet
########
sites <- sites[sites$maf > 0,]

########
# Filter based on minor allele frequency
########
min_maf=0.05
n_nearby=250 # number of neutral loci to compare candidates against


info <- sites[sites$maf >= min_maf & sites$type != "neutral",]
info_neut <- sites[sites$type == "neutral",]
info_neut$maf_rank <- rank(info_neut$maf)

########
# calculate p-values using method
########
library(parallel)
calc_pval <- function(site) {
  input_maf <- info$maf[site]; input_fst <- info$fst[site]
  n_less <- sum(info_neut$maf <= input_maf)
  
  tmp_fst <- info_neut$fst[info_neut$maf_rank <= (n_less + n_nearby/2) & info_neut$maf_rank >= n_less - n_nearby/2]
  if (length(tmp_fst) < n_nearby) {tmp_fst <- info_neut$fst[info_neut$maf_rank >= length(info_neut$maf_rank) - n_nearby]}
  
  tmp_percentiles <- ecdf(tmp_fst)
  tmp_p <- 1-tmp_percentiles(input_fst)
  return(tmp_p)
}

info$fst_pval <- do.call("c",mclapply(1:nrow(info), FUN = calc_pval))

rm(sites, calc_pval, nindiv_post, nindiv_pre)
########

########
# enrichment, binned by maf
########
library(ggplot2); library(patchwork)

## we need to get a matrix that has the maf window, the percentile, and the degree of enrichment
tmp_bin_breaks <- c(0.1,.2,.3,.4,0.5)
tmp_bins <-  matrix(ncol=2, nrow=length(tmp_bin_breaks)); tmp_bins[,1] <- c(0,tmp_bin_breaks[-length(tmp_bin_breaks)]); tmp_bins[,2] <- c(tmp_bin_breaks); rm(tmp_bin_breaks)

print("number of sites in each maf bin")
for (i in 1:nrow(tmp_bins)) {
  print(paste(100*tmp_bins[i,1], "% to ", 100*tmp_bins[i,2], "%: ", 
              nrow(info[info$maf > tmp_bins[i,1] & info$maf <= tmp_bins[i,2],]), " sites", sep=""))
}

res <- NULL; for (tmp_enrich in c(0.005, seq(0.01,0.1,0.01), seq(0.1,0.2,0.05))) {
  for (tmp_bin in 1:nrow(tmp_bins)) {
    tmp_data <- info[info$maf > tmp_bins[tmp_bin,1] & info$maf <= tmp_bins[tmp_bin,2],]
    tmp_res <- as.data.frame(matrix(ncol=6, nrow=1)); colnames(tmp_res) <- c('maf_low', 'maf_high', 'population', 'enrichment', 'observed', 'expected')
    tmp_res$maf_low <- tmp_bins[tmp_bin,1]; tmp_res$maf_high <- tmp_bins[tmp_bin,2]
    tmp_res$enrichment <- tmp_enrich
    tmp_res$observed <- sum(tmp_data$fst_pval < tmp_enrich)
    tmp_res$expected <- nrow(tmp_data)*tmp_enrich
    rbind(res, tmp_res) -> res; rm(tmp_res, tmp_data)
  }; rm(tmp_bin)
}; rm(tmp_enrich, tmp_bins)

res$maf_bin <- paste(100*res$maf_low, "% to ", 100*res$maf_high, "%", sep="")
res$fc <- log2(res$observed/res$expected)
res$perm <- num



print("binomial test: 99%")
tmp <- res[res$maf_high > 0 & res$enrichment == 0.01,]; tmp
binom.test(sum(tmp$observed), n = (1/.01)*sum(tmp$expected), p=0.01); log2(sum(tmp$observed)/sum(tmp$expected)); sum(tmp$observed)/sum(tmp$expected)


res_sim <- res; rm(res)
res_sim$fc[res_sim$fc < 0] <- 0

perm_results <- ggplot(res_sim, aes(x=1-enrichment , y=fc, color=maf_bin)) +
  geom_point(size=2) + geom_line() +
  scale_color_manual(values=c("#999999", "#D55E00", "#CC79A7", "#56B4E9", "#0072B2")) +
  theme_classic() + ylab("log2(enrichment)") +
  ggtitle("pre vs post") + coord_cartesian(ylim=c(-0.25,4)) +
  geom_abline(slope=0, intercept=0, color='dark gray') + theme(legend.position = "bottom") + 
  guides(fill = guide_legend(label.position = "bottom")) + xlab("percentile of neutral variants") + 
  ggtitle("simulated")
perm_results



load("./fig2_image.RData")

library(patchwork) 
cbPalette <- c("#999999", "#D55E00", "#CC79A7", "#56B4E9", "#0072B2")
real_results <- panel_enrichment + ggtitle('observed') + coord_cartesian(ylim=c(-0.1,6)) + scale_color_manual(values=cbPalette) 
perm_results <- perm_results + coord_cartesian(ylim=c(-0.1,6))


perm_results + real_results + plot_layout(guides='collect') & theme(legend.position = 'bottom')
  
 
