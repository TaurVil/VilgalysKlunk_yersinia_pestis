library(data.table); library(ggplot2); library(patchwork)

load("./candidate_loci.RData")
rm(method, min_maf, min_n, min_neutral_nsites, tmp_n, window_size, drop_samples_missing, mafs, bins, design2)

rsids <- candidates
sites <- D_replicated_in_Denmark$site #  "chr5_96244585",

likelihood <- function(p, data){
  gt.freq <- c((1-p)^2, 2*p*(1-p), p^2)
  ll <- sum(log(rowSums(t(t(data)*gt.freq))))
  return(ll)
}
estimate_af_ml<-function(d){
  n<-NROW(d)
  af<-rep(NA, n)
  for(i in 1:n){
    gl<-d[i,,drop=FALSE]
    gl <- gl[!is.na(gl)]
    gl <- matrix(gl, ncol=3, byrow=T)
    gl <- gl[!is.na(gl[,1]),,drop=FALSE]
    opt <- optimize(likelihood, interval=c(0,1), maximum=TRUE, data=gl)
    af[i] <- opt$maximum
  }
  return(af)
}

info <- D_replicated_in_Denmark

plot_info <- NULL
design <- expand.grid(time = gl(3, 1, labels = c("pre", "post", "during")), pop = gl(2, 1, labels = c("denmark", "london")))

sites_gwas <- D_replicated_in_Denmark$site[D_replicated_in_Denmark$type == "gwas"]
for (ttt in c(4:nrow(design),1:2)) {
  name=paste("~/GitHub/VilgalysKlunk_response_to_commentary/DATA/genoliks/genolik.gwas_",design$pop[ttt],"_", design$time[ttt],".genolik",sep="")
  n2=paste(design$pop[ttt],design$time[ttt],sep="_")
  d <- as.data.frame(fread(name))
  row.names(d) <- paste(d$V1, d$V2, sep="_")
  d <- d[,-c(1:2)]
  d <- d[row.names(d) %in% sites_gwas,]
  d[d == -1] <- NA; n <- ncol(d)/3
  
  #### bootstrap step 
  boot_data <- NULL
  for (i in 1:1000) {
    k <- sample(1:n, n, replace=T)
    new_d <- as.data.frame(matrix(ncol=1, nrow=nrow(d)))
    for (j in k) {
      new_d <- cbind(new_d, d[,c((j*3 - 2):(j*3))])
    }; rm(j)
    new_d[,1] <- row.names(d)
    colnames(new_d) <- c("site", paste0("V",1:(3*n)))
    boot_data <- rbind(boot_data, new_d); rm(new_d, k)
  }; rm(i)
  rm(d)
  
  boot_info <- as.data.frame(boot_data$site)
  colnames(boot_info)[1] <- "site"
  boot_info$freq <- estimate_af_ml(boot_data[,-1])
  
  tmp_info <- info[info$site %in% sites_gwas,c(1:5)]
  tmp_info$pop <- design$pop[ttt]
  tmp_info$time <- design$time[ttt]
  tmp_info$real_af <- info[[paste0(design$pop[ttt],".",design$time[ttt],".alternate.ML")]][info$site %in% sites_gwas]
  
  for (s in sites_gwas) {
    mean(boot_info$freq[boot_info$site == s]) -> tmp_info$boot.mean[tmp_info$site == s]
    quantile(boot_info$freq[boot_info$site == s], 0.05) -> tmp_info$boot.lower95[tmp_info$site == s]
    quantile(boot_info$freq[boot_info$site == s], 0.95) -> tmp_info$boot.upper95[tmp_info$site == s]
    sd(boot_info$freq[boot_info$site == s]) -> tmp_info$boot.sd[tmp_info$site == s]
  }; rm(s)
  print(design[ttt,])
  
  rbind(plot_info, tmp_info) -> plot_info
  
  rm(n, n2, name, tmp_info, boot_data, boot_info)
  
}

sites_exons <- D_replicated_in_Denmark$site[D_replicated_in_Denmark$type == "exon"]
if (length(sites_exons) > 0) {for (ttt in c(4:nrow(design),1:2)) {
  name=paste("~/GitHub/VilgalysKlunk_response_to_commentary/DATA/genoliks/genolik.exons_",design$pop[ttt],"_", design$time[ttt],".genolik",sep="")
  n2=paste(design$pop[ttt],design$time[ttt],sep="_")
  d <- as.data.frame(fread(name))
  row.names(d) <- paste(d$V1, d$V2, sep="_")
  d <- d[,-c(1:2)]
  d <- d[row.names(d) %in% sites_exons,]
  d[d == -1] <- NA; n <- ncol(d)/3
  
  #### bootstrap step 
  boot_data <- NULL
  for (i in 1:1000) {
    k <- sample(1:n, n, replace=T)
    new_d <- as.data.frame(matrix(ncol=1, nrow=nrow(d)))
    for (j in k) {
      new_d <- cbind(new_d, d[,c((j*3 - 2):(j*3))])
    }; rm(j)
    new_d[,1] <- row.names(d)
    colnames(new_d) <- c("site", paste0("V",1:(3*n)))
    boot_data <- rbind(boot_data, new_d); rm(new_d, k)
  }; rm(i)
  rm(d)
  
  boot_info <- as.data.frame(boot_data$site)
  colnames(boot_info)[1] <- "site"
  boot_info$freq <- estimate_af_ml(boot_data[,-1])
  
  tmp_info <- info[info$site %in% sites_exons,c(1:5)]
  tmp_info$pop <- design$pop[ttt]
  tmp_info$time <- design$time[ttt]
  tmp_info$real_af <- info[[paste0(design$pop[ttt],".",design$time[ttt],".alternate.ML")]][info$site %in% sites_exons]
  
  for (s in sites_exons) {
    mean(boot_info$freq[boot_info$site == s]) -> tmp_info$boot.mean[tmp_info$site == s]
    quantile(boot_info$freq[boot_info$site == s], 0.05) -> tmp_info$boot.lower95[tmp_info$site == s]
    quantile(boot_info$freq[boot_info$site == s], 0.95) -> tmp_info$boot.upper95[tmp_info$site == s]
    sd(boot_info$freq[boot_info$site == s]) -> tmp_info$boot.sd[tmp_info$site == s]
  }; rm(s)
  print(design[ttt,])
  
  rbind(plot_info, tmp_info) -> plot_info
  
  rm(n, n2, name, tmp_info, boot_data, boot_info)
  
}}

rm(ttt, sites_gwas, sites_exons)


plot_info$time <- factor(plot_info$time, levels=c("pre", "during", "post"))
plot_info$xpos <- as.numeric(plot_info$time)
plot_info$xpos[plot_info$pop == 'denmark'] <- plot_info$xpos[plot_info$pop == 'denmark'] + 0.1

plot_info[,c(8:11)] <- 1-plot_info[,c(8:11)]
plot(data=plot_info, boot.mean ~ real_af); abline(a=0,b=1,col='red')

plot1 <- ggplot(data=plot_info, aes(x=xpos, y=boot.mean, col=pop)) +
  facet_grid( . ~ site ) + 
  geom_point(size=4) + 
  geom_line(linetype="dashed") + 
  coord_cartesian(ylim=c(0.2,1)) +
  # geom_segment(aes(x=xpos, xend=xpos, y=boot.lower95, yend=boot.upper95), size=1.1) +
  # geom_segment(aes(x=xpos, xend=xpos, y=real_af - boot.sd, yend=real_af + boot.sd), size=1.1) +
  geom_segment(aes(x=xpos, xend=xpos, y=boot.mean - boot.sd, yend=boot.mean + boot.sd), linewidth=1.1) +
  # ggtitle("linked ERAP2 variants") +
  scale_color_manual(values=c("#56B4E9", 'red', "#56B4E9", 'red')) + 
  theme_classic() + theme(legend.position = 'inset') + xlab("time") + ylab("derived allele frequency") 
plot1


panel_ERAP_loci <- ggplot(data=plot_info, aes(x=xpos, y=boot.mean, col=pop)) +
  geom_point(size=4) + 
  geom_line(linetype="dashed") + 
  coord_cartesian(ylim=c(0.2,1)) +
  # geom_segment(aes(x=xpos, xend=xpos, y=boot.lower95, yend=boot.upper95), size=1.1) +
  # geom_segment(aes(x=xpos, xend=xpos, y=real_af - boot.sd, yend=real_af + boot.sd), size=1.1) +
  geom_segment(aes(x=xpos, xend=xpos, y=boot.mean - boot.sd, yend=boot.mean + boot.sd), linewidth=1.1) +
  # ggtitle("linked ERAP2 variants") +
  scale_color_manual(values=c("#56B4E9", 'red', "#56B4E9", 'red')) + 
  ggtitle("rs2549794") +
  theme_classic() + theme(legend.position = 'inset') + xlab("time") + ylab("derived allele frequency") ; panel_ERAP_loci


cbPalette <- c("#999999", "#D55E00", "#CC79A7", "#56B4E9", "#0072B2")
panel_enrichment <- panel_enrichment + ggtitle('London') + coord_cartesian(ylim=c(-0.1,6)) + scale_color_manual(values=cbPalette) 
panel_D_enrichment <- panel_D_enrichment + ggtitle('Denmark') + coord_cartesian(ylim=c(-0.1,6)) + scale_color_manual(values=cbPalette) 

((panel_enrichment + panel_D_enrichment + plot_layout(guides = "collect") & theme(legend.position = "bottom")) + panel_ERAP_loci) / panel_manhattan  + 
  plot_annotation(tag_levels = 'A') & 
  theme(text = element_text(color="black", family = 'sans'))

save.image("./fig2_image.RData")

