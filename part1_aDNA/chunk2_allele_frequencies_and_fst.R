######## 
# Get study design matrix and R packages
######## 
library(data.table)
# Create design matrix for which we'll pull in data
design <- expand.grid(time = gl(3, 1, labels = c("pre", "post", "during")), pop = gl(2, 1, labels = c("denmark", "london")))
design <- design[-3,]
                      #type = gl(3, 1, labels = c("exons", "gwas", "neutral"))
design2 <- expand.grid(time1 = gl(3, 1, labels = c("pre", "post", "during")), time2 = gl(3, 1, labels = c("pre", "post", "during")), pop = gl(2, 1, labels = c("denmark", "london")))
design2<-subset(design2, !(design2$time1 == design2$time2)); design2 <- design2[c(1,7:8,10),]

########

min_n=10; min_maf=0.05
drop_samples_missing=.5

method='sliding'

# for method=bins
if (method=='bins') {
  bin_breaks <- c(0.07, 0.1, 0.25, 0.5)
  bins <- matrix(ncol=2, nrow=length(bin_breaks)); bins[,1] <- c(0,bin_breaks[-length(bin_breaks)]); bins[,2] <- c(bin_breaks)
}
# for method=sliding: sliding windows centered on each percentage with at least 200 nearby variants
if (method == 'sliding') {
  min_neutral_nsites <- 200
  if (min_neutral_nsites == '200') {
    mafs <- seq(0.005,0.40, 0.01); window_size <- 0.05*rep(1, length(mafs))
    window_size[mafs>=0.135] <- 0.06; window_size[mafs>=0.145] <- 0.07; window_size[mafs>=0.155] <- 0.08
    window_size[mafs>=0.165] <- 0.09; window_size[mafs>=0.175] <- 0.11; window_size[mafs>=0.185] <- 0.12
    window_size[mafs>=0.205] <- 0.13; window_size[mafs>=0.225] <- 0.15; window_size[mafs>=0.255] <- 0.16
    window_size[mafs>=0.275] <- 0.17; window_size[mafs>=0.305] <- 0.18; window_size[mafs>=0.335] <- 0.19
    window_size[mafs>=0.345] <- 0.2; window_size[mafs>=0.365] <- 0.21; window_size[mafs>=0.375] <- 0.22
    window_size[mafs>=0.395] <- 0.24
    bins <- cbind(mafs-0.005, mafs+0.005, mafs-window_size/2, mafs+window_size/2)
    bins[nrow(bins),2] <- 0.5
  }
  # we used the following line to check this once we had the set of neutral sites
  # for (i in 1:length(mafs)) { print(paste("maf: ", mafs[i], "; nearby sites: ", sum(info_neut$maf > mafs[i] - window_size[i]/2 & info_neut$maf <= mafs[i] + window_size[i]/2), sep=""))}
}
######## 

######## 
# Get trimmed aDNA frequencies
######## 
# Info files with the number and alternate allele freq in each population
# mac >= 3; biallelic; minQ >= 30
# exons filtered for sites within bed file; rest mapped to the bed file
######## 
info_gwas <- fread("./DATA/genoliks/immune.frq")[,-c(3:4)]
info_neut <- fread("./DATA/genoliks/neutral.frq")[,-c(3:4)]
info_exon <- fread("./DATA/genoliks/exon.frq")[,-c(3:4)]

colnames(info_gwas) <- c("chr", "pos", "ref", "alt")
colnames(info_neut) <- c("chr", "pos", "ref", "alt")
colnames(info_exon) <- c("chr", "pos", "ref", "alt")

info_gwas$ref <- gsub(":.*", "",info_gwas$ref); info_gwas$alt <- gsub(":.*", "",info_gwas$alt)
info_neut$ref <- gsub(":.*", "",info_neut$ref); info_neut$alt <- gsub(":.*", "",info_neut$alt)
info_exon$ref <- gsub(":.*", "",info_exon$ref); info_exon$alt <- gsub(":.*", "",info_exon$alt)

paste(info_gwas$chr,info_gwas$pos,sep="_") -> info_gwas$site
paste(info_neut$chr,info_neut$pos,sep="_") -> info_neut$site
paste(info_exon$chr,info_exon$pos,sep="_") -> info_exon$site


for (i in c(3:nrow(design),1:2)) {
  # read in exon data
  name=paste("./DATA/genoliks/genolik.exons_",design$pop[i],"_", design$time[i],".genolik",sep="")
  n2=paste(design$pop[i],design$time[i],sep="_"); d <- as.data.frame(fread(name)); rm(name)
  snames=read.delim(paste("./DATA/SampleNames/", design$pop[i], "_", design$time[i],"_exons.012.indv",sep=""), header=F)
  ncol(d)-2 == 3*nrow(snames)
  if (i == 3) {paste(d[,1],d[,2],sep="_") -> info_exon$site_check}
  # remove site. change missing data to NA. get number of samples. 
  d <- d[,-(1:2)] ; d[d == -1] <- NA; n <- ncol(d)/3
  # report data on missingness per sample
  assign(paste("missing",n2, "exon", sep="_"), value = (colSums(is.na(d))/nrow(d))[seq(1,3*n,3)])
  # remove samples missing too much data
  keep <- which(get(paste("missing",n2, "exon", sep="_")) <= drop_samples_missing)
  assign(paste("keep",n2,"exon",sep="_"), value=snames[keep,]); write.table(snames[keep,], paste("./DATA/SampleNames/keep",n2,"exon",sep="_"), row.names=F, col.names=F, sep="\t", quote=F)
  k2 <- as.data.frame(matrix(nrow=nrow(d), ncol=1)); for (k in keep) {k2 <- cbind(k2, d[,(k*3-2):(k*3)])}
  d <- k2[,-1]; rm(keep, k2, k); n <- ncol(d)/3
  # get the allele frequency, and the number of individuals it was called from
  info_exon[[paste(n2,"alternate",sep=".")]] <- rowMeans(d[,seq(3,n*3,3)],na.rm=T) + rowMeans(d[,seq(2,n*3,3)],na.rm=T)/2
  info_exon[[paste(n2,"called",sep=".")]] <- (n*3-rowSums(is.na(d)))/3
  # Report number of samples included 
  assign(paste("samples",n2, "exon", sep="_"), n)
  # cleanup
  rm(d, n)
  
  # read in gwas data
  name=paste("./DATA/genoliks/genolik.gwas_",design$pop[i],"_", design$time[i],".genolik",sep="")
  n2=paste(design$pop[i],design$time[i],sep="_"); d <- as.data.frame(fread(name)); rm(name)
  snames=read.delim(paste("./DATA/SampleNames/", design$pop[i], "_", design$time[i],"_immune.012.indv",sep=""), header=F)
  ncol(d)-2 == 3*nrow(snames)
  if (i == 3) {paste(d[,1],d[,2],sep="_") -> info_gwas$site_check}
  # remove site. change missing data to NA. get number of samples. 
  d <- d[,-(1:2)] ; d[d == -1] <- NA; n <- ncol(d)/3
  # report data on missingness per sample
  assign(paste("missing",n2, "gwas", sep="_"), value = (colSums(is.na(d))/nrow(d))[seq(1,3*n,3)])
  # remove samples missing too much data
  keep <- which(get(paste("missing",n2, "gwas", sep="_")) <= drop_samples_missing)
  assign(paste("keep",n2,"gwas",sep="_"), value=snames[keep,]); write.table(snames[keep,], paste("./DATA/SampleNames/keep",n2,"gwas",sep="_"), row.names=F, col.names=F, sep="\t", quote=F)
  k2 <- as.data.frame(matrix(nrow=nrow(d), ncol=1)); for (k in keep) {k2 <- cbind(k2, d[,(k*3-2):(k*3)])}
  d <- k2[,-1]; rm(keep, k2, k); n <- ncol(d)/3
  # get the allele frequency, and the number of individuals it was called from
  info_gwas[[paste(n2,"alternate",sep=".")]] <- rowMeans(d[,seq(3,ncol(d),3)],na.rm=T) + rowMeans(d[,seq(2,ncol(d),3)],na.rm=T)/2
  info_gwas[[paste(n2,"called",sep=".")]] <- (ncol(d)-rowSums(is.na(d)))/3
  # Report number of samples included 
  assign(paste("samples",n2, "gwas", sep="_"), n)
  # cleanup
  rm(d, n)
  
  # read in neutral data
  name=paste("./DATA/genoliks/genolik.neutral_",design$pop[i],"_", design$time[i],".genolik",sep="")
  n2=paste(design$pop[i],design$time[i],sep="_"); d <- as.data.frame(fread(name)); rm(name)
  snames=read.delim(paste("./DATA/SampleNames/", design$pop[i], "_", design$time[i],"_neutral.012.indv",sep=""), header=F)
  ncol(d)-2 == 3*nrow(snames)
  if (i == 3) {paste(d[,1],d[,2],sep="_") -> info_neut$site_check}
  # remove site. change missing data to NA. get number of samples.
  d <- d[,-(1:2)] ; d[d == -1] <- NA; n <- ncol(d)/3
  # report data on missingness per sample
  assign(paste("missing",n2, "neutral", sep="_"), value = (colSums(is.na(d))/nrow(d))[seq(1,3*n,3)])
  # remove samples missing too much data
  keep <- which(get(paste("missing",n2, "neutral", sep="_")) <= drop_samples_missing)
  assign(paste("keep",n2,"neutral",sep="_"), value=snames[keep,]); write.table(snames[keep,], paste("./DATA/SampleNames/keep",n2,"neutral",sep="_"), row.names=F, col.names=F, sep="\t", quote=F)
  k2 <- as.data.frame(matrix(nrow=nrow(d), ncol=1)); for (k in keep) {k2 <- cbind(k2, d[,(k*3-2):(k*3)])}
  d <- k2[,-1]; rm(keep, k2, k); n <- ncol(d)/3
  # get the allele frequency, and the number of individuals it was called from
  info_neut[[paste(n2,"alternate",sep=".")]] <- rowMeans(d[,seq(3,ncol(d),3)],na.rm=T) + rowMeans(d[,seq(2,ncol(d),3)],na.rm=T)/2
  info_neut[[paste(n2,"called",sep=".")]] <- (ncol(d)-rowSums(is.na(d)))/3
  # Report number of samples included 
  assign(paste("samples",n2, "neutral", sep="_"), n)
  # cleanup
  rm(d, n)
}; rm(i,n2)
save.image("./r01.summarized_allele_frequencies.RData")
## saved for John & Matthias

sum(info_gwas$site == info_gwas$site_check) == nrow(info_gwas)
sum(info_exon$site == info_exon$site_check) == nrow(info_exon)
sum(info_neut$site == info_neut$site_check) == nrow(info_neut)

for (f in ls(pattern="missing")) {print(paste(f,sum(get(f) < 0.5), sep=": "))}; rm(f)


########

########
# calculate alternate and minor allele frequency
# calculate a "London" frequency and a "Denmark" frequency as the mean between time points within that population
# NOTE: THIS TAKES US UP TO THE FIRST 20 COLUMNS
########
# let's do gwas first
info_gwas$london.alternate <- rowMeans(cbind(info_gwas$london_pre.alternate, info_gwas$london_during.alternate, info_gwas$london_post.alternate))
info_gwas$denmark.alternate <- rowMeans(cbind(info_gwas$denmark_pre.alternate, info_gwas$denmark_post.alternate))
info_gwas$alternate <- rowMeans(cbind(info_gwas$london.alternate, info_gwas$denmark.alternate))

# let's do exons
info_exon$london.alternate <- rowMeans(cbind(info_exon$london_pre.alternate, info_exon$london_during.alternate, info_exon$london_post.alternate))
info_exon$denmark.alternate <- rowMeans(cbind(info_exon$denmark_pre.alternate, info_exon$denmark_post.alternate))
info_exon$alternate <- rowMeans(cbind(info_exon$london.alternate, info_exon$denmark.alternate))

# let's do neutral
info_neut$london.alternate <- rowMeans(cbind(info_neut$london_pre.alternate, info_neut$london_during.alternate, info_neut$london_post.alternate))
info_neut$denmark.alternate <- rowMeans(cbind(info_neut$denmark_pre.alternate, info_neut$denmark_post.alternate))
info_neut$alternate <- rowMeans(cbind(info_neut$london.alternate, info_neut$denmark.alternate))

# add in maf as well as the alternate allele frequency
info_gwas$maf <- apply(cbind(info_gwas$alternate, 1-info_gwas$alternate), 1, min)
info_exon$maf <- apply(cbind(info_exon$alternate, 1-info_exon$alternate), 1, min)
info_neut$maf <- apply(cbind(info_neut$alternate, 1-info_neut$alternate), 1, min)


########

########
# Calculate Fst
########
for (i in 1:nrow(design)) {
  n2=paste(design$pop[i],design$time[i],sep="_")
  # expected heterozygosity
  info_gwas[[paste(n2,"ehet",sep=".")]] <- 2*info_gwas[[paste(n2,"alternate",sep=".")]]*(1-info_gwas[[paste(n2,"alternate",sep=".")]])
  info_exon[[paste(n2,"ehet",sep=".")]] <- 2*info_exon[[paste(n2,"alternate",sep=".")]]*(1-info_exon[[paste(n2,"alternate",sep=".")]])
  info_neut[[paste(n2,"ehet",sep=".")]] <- 2*info_neut[[paste(n2,"alternate",sep=".")]]*(1-info_neut[[paste(n2,"alternate",sep=".")]])
}; rm(i, n2)
# pairwise means, pairwise expectations, and Fst
for (i in 1:nrow(design2)) {
  n2=paste(design2$pop[i],design2$time1[i],design2$time2[i],sep="_")
  info_gwas[[paste(n2,"mean",sep=".")]] <- rowMeans(cbind(info_gwas[[paste(design2$pop[i],"_",design2$time1[i],".alternate",sep="")]], info_gwas[[paste(design2$pop[i],"_",design2$time2[i],".alternate",sep="")]]))
  info_gwas[[paste(n2,"ehet",sep=".")]] <- 2*info_gwas[[paste(n2,"mean",sep=".")]]*(1-info_gwas[[paste(n2,"mean",sep=".")]])
  info_gwas[[paste(n2,"fst",sep=".")]] <- (info_gwas[[paste(n2,"ehet",sep=".")]] - (info_gwas[[paste(design2$pop[i],"_",design2$time1[i],".ehet",sep="")]] + info_gwas[[paste(design2$pop[i],"_",design2$time2[i],".ehet",sep="")]])/2)/info_gwas[[paste(n2,"ehet",sep=".")]]
  
  info_exon[[paste(n2,"mean",sep=".")]] <- rowMeans(cbind(info_exon[[paste(design2$pop[i],"_",design2$time1[i],".alternate",sep="")]], info_exon[[paste(design2$pop[i],"_",design2$time2[i],".alternate",sep="")]]))
  info_exon[[paste(n2,"ehet",sep=".")]] <- 2*info_exon[[paste(n2,"mean",sep=".")]]*(1-info_exon[[paste(n2,"mean",sep=".")]])
  info_exon[[paste(n2,"fst",sep=".")]] <- (info_exon[[paste(n2,"ehet",sep=".")]] - (info_exon[[paste(design2$pop[i],"_",design2$time1[i],".ehet",sep="")]] + info_exon[[paste(design2$pop[i],"_",design2$time2[i],".ehet",sep="")]])/2)/info_exon[[paste(n2,"ehet",sep=".")]]
  
  info_neut[[paste(n2,"mean",sep=".")]] <- rowMeans(cbind(info_neut[[paste(design2$pop[i],"_",design2$time1[i],".alternate",sep="")]], info_neut[[paste(design2$pop[i],"_",design2$time2[i],".alternate",sep="")]]))
  info_neut[[paste(n2,"ehet",sep=".")]] <- 2*info_neut[[paste(n2,"mean",sep=".")]]*(1-info_neut[[paste(n2,"mean",sep=".")]])
  info_neut[[paste(n2,"fst",sep=".")]] <- (info_neut[[paste(n2,"ehet",sep=".")]] - (info_neut[[paste(design2$pop[i],"_",design2$time1[i],".ehet",sep="")]] + info_neut[[paste(design2$pop[i],"_",design2$time2[i],".ehet",sep="")]])/2)/info_neut[[paste(n2,"ehet",sep=".")]]
  
}; rm(i,n2)
info_gwas -> info_gwas_trim; rm(info_gwas)
info_neut -> info_neut_trim; rm(info_neut)
info_exon -> info_exon_trim; rm(info_exon)
######## 

######## 
# Filter for 10 individual per time point
######## 
# Report number of sites before filtering
nrow(info_exon_trim); nrow(info_gwas_trim); nrow(info_neut_trim)
min_n <- 10
######## 
attach(info_exon_trim)
info_exon_trim <- subset(info_exon_trim, apply(X = cbind(london_during.called, london_pre.called, london_post.called, denmark_pre.called, denmark_post.called), 1, min) >= min_n)
detach(info_exon_trim)

attach(info_neut_trim)
info_neut_trim <- subset(info_neut_trim, apply(X = cbind(london_during.called, london_pre.called, london_post.called, denmark_pre.called, denmark_post.called), 1, min) >= min_n)
detach(info_neut_trim)

attach(info_gwas_trim)
info_gwas_trim <- subset(info_gwas_trim, apply(X = cbind(london_during.called, london_pre.called, london_post.called, denmark_pre.called, denmark_post.called), 1, min) >= min_n)
detach(info_gwas_trim)
######## 
# Report number of sites after filtering
nrow(info_exon_trim); nrow(info_gwas_trim); nrow(info_neut_trim)
info_gwas <- info_gwas_trim; rm(info_gwas_trim)
info_exon <- info_exon_trim; rm(info_exon_trim)
info_neut <- info_neut_trim; rm(info_neut_trim)
######## 

########
# Get difference between time points, then p-values and candidate loci
########
# differences between time points
info_exon$delta_L13 <- info_exon$london_post.alternate-info_exon$london_pre.alternate
info_exon$delta_L12 <- info_exon$london_during.alternate-info_exon$london_pre.alternate
info_exon$delta_L23 <- info_exon$london_post.alternate-info_exon$london_during.alternate
info_exon$delta_D13 <- info_exon$denmark_post.alternate-info_exon$denmark_pre.alternate

info_neut$delta_L13 <- info_neut$london_post.alternate-info_neut$london_pre.alternate
info_neut$delta_L12 <- info_neut$london_during.alternate-info_neut$london_pre.alternate
info_neut$delta_L23 <- info_neut$london_post.alternate-info_neut$london_during.alternate
info_neut$delta_D13 <- info_neut$denmark_post.alternate-info_neut$denmark_pre.alternate

info_gwas$delta_L13 <- info_gwas$london_post.alternate-info_gwas$london_pre.alternate
info_gwas$delta_L12 <- info_gwas$london_during.alternate-info_gwas$london_pre.alternate
info_gwas$delta_L23 <- info_gwas$london_post.alternate-info_gwas$london_during.alternate
info_gwas$delta_D13 <- info_gwas$denmark_post.alternate-info_gwas$denmark_pre.alternate

########
# trim out some of the columns: leaving just Fst, site, maf
# site information: 1:5
# alternate and minor allele frequency: 7, 9,11,13,15, 19, 20
# fst: 28, 31, 34, 37
# delta_AF: 38:41
########
info_exon <- info_exon[,c(1:5, 7, 9, 11, 13, 15, 19, 20, 28, 31, 34, 37, 38:41)]
info_neut <- info_neut[,c(1:5, 7, 9, 11, 13, 15, 19, 20, 28, 31, 34, 37, 38:41)]
info_gwas <- info_gwas[,c(1:5, 7, 9, 11, 13, 15, 19, 20, 28, 31, 34, 37, 38:41)]
########

########
# Setup to get p-values
########
# group functional sites together
########
info_neut$type <- "neutral"
info_exon$type <- "exon"
info_gwas$type <- "gwas"

rbind(info_exon, info_gwas) -> info; rm(info_exon, info_gwas)
########

########
# Filter based on minor allele frequency
########
info_neut <- subset(info_neut, info_neut$maf > min_maf)
info <- subset(info, info$maf > min_maf)
########

pvals <- info[,1:5]
neut_pvals <- info_neut[,1:5]
########
# calculate p-values using method
########
library(parallel)

calc_pval_L13 <- function(site) {
  input_maf <- info$maf[site]
  input_fst <- info$london_post_pre.fst[site]
  neut_maf <- info_neut$maf
  neut_fst <- info_neut$london_post_pre.fst
  
  if (method == "bins") {
    tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
    tmp_fst <- neut_fst[neut_maf > tmp_bin[1] & neut_maf <= tmp_bin[2]]
  }
  if (method == "sliding") {
    tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
    tmp_fst <- neut_fst[neut_maf > tmp_bin[3] & neut_maf <= tmp_bin[4]]
  }

  tmp_percentiles <- ecdf(tmp_fst)
  tmp_p <- 1-tmp_percentiles(input_fst)
  return(tmp_p)
}
pvals$L13.pval_split <- do.call("c",mclapply(1:nrow(info), FUN = calc_pval_L13))
calc_pval_L12 <- function(site) {
  input_maf <- info$maf[site]
  input_fst <- info$london_during_pre.fst[site]
  neut_maf <- info_neut$maf
  neut_fst <- info_neut$london_during_pre.fst
  
  if (method == "bins") {
    tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
    tmp_fst <- neut_fst[neut_maf > tmp_bin[1] & neut_maf <= tmp_bin[2]]
  }
  if (method == "sliding") {
    tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
    tmp_fst <- neut_fst[neut_maf > tmp_bin[3] & neut_maf <= tmp_bin[4]]
  }
  
  tmp_percentiles <- ecdf(tmp_fst)
  tmp_p <- 1-tmp_percentiles(input_fst)
  return(tmp_p)
}
pvals$L12.pval_split <- do.call("c",mclapply(1:nrow(info), FUN = calc_pval_L12))
calc_pval_L23 <- function(site) {
  input_maf <- info$maf[site]
  input_fst <- info$london_during_post.fst[site]
  neut_maf <- info_neut$maf
  neut_fst <- info_neut$london_during_post.fst
  
  if (method == "bins") {
    tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
    tmp_fst <- neut_fst[neut_maf > tmp_bin[1] & neut_maf <= tmp_bin[2]]
  }
  if (method == "sliding") {
    tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
    tmp_fst <- neut_fst[neut_maf > tmp_bin[3] & neut_maf <= tmp_bin[4]]
  }
  
  tmp_percentiles <- ecdf(tmp_fst)
  tmp_p <- 1-tmp_percentiles(input_fst)
  return(tmp_p)
}
pvals$L23.pval_split <- do.call("c",mclapply(1:nrow(info), FUN = calc_pval_L23))

calc_pval_D13 <- function(site) {
  input_maf <- info$maf[site]
  input_fst <- info$denmark_post_pre.fst[site]
  neut_maf <- info_neut$maf
  neut_fst <- info_neut$denmark_post_pre.fst
  
  if (method == "bins") {
    tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
    tmp_fst <- neut_fst[neut_maf > tmp_bin[1] & neut_maf <= tmp_bin[2]]
  }
  if (method == "sliding") {
    tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
    tmp_fst <- neut_fst[neut_maf > tmp_bin[3] & neut_maf <= tmp_bin[4]]
  }
  
  tmp_percentiles <- ecdf(tmp_fst)
  tmp_p <- 1-tmp_percentiles(input_fst)
  return(tmp_p)
}
pvals$D13.pval_split <- do.call("c",mclapply(1:nrow(info), FUN = calc_pval_D13))

{
  calc_neut_pval_L13 <- function(site) {
    input_maf <- info_neut$maf[site]
    input_fst <- info_neut$london_post_pre.fst[site]
    neut_maf <- info_neut$maf
    neut_fst <- info_neut$london_post_pre.fst
    
    if (method == "bins") {
      tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
      tmp_fst <- neut_fst[neut_maf > tmp_bin[1] & neut_maf <= tmp_bin[2]]
    }
    if (method == "sliding") {
      tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
      tmp_fst <- neut_fst[neut_maf > tmp_bin[3] & neut_maf <= tmp_bin[4]]
    }
    tmp_percentiles <- ecdf(tmp_fst)
    tmp_p <- 1-tmp_percentiles(input_fst)
    return(tmp_p)
  }
  neut_pvals$L13.pval_split <- do.call("c",mclapply(1:nrow(info_neut), FUN = calc_neut_pval_L13))
  calc_neut_pval_L12 <- function(site) {
    input_maf <- info_neut$maf[site]
    input_fst <- info_neut$london_during_pre.fst[site]
    neut_maf <- info_neut$maf
    neut_fst <- info_neut$london_during_pre.fst
    
    if (method == "bins") {
      tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
      tmp_fst <- neut_fst[neut_maf > tmp_bin[1] & neut_maf <= tmp_bin[2]]
    }
    if (method == "sliding") {
      tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
      tmp_fst <- neut_fst[neut_maf > tmp_bin[3] & neut_maf <= tmp_bin[4]]
    }
    
    tmp_percentiles <- ecdf(tmp_fst)
    tmp_p <- 1-tmp_percentiles(input_fst)
    return(tmp_p)
  }
  neut_pvals$L12.pval_split <- do.call("c",mclapply(1:nrow(info_neut), FUN = calc_neut_pval_L12))
  calc_neut_pval_L23 <- function(site) {
    input_maf <- info_neut$maf[site]
    input_fst <- info_neut$london_during_post.fst[site]
    neut_maf <- info_neut$maf
    neut_fst <- info_neut$london_during_post.fst
    
    if (method == "bins") {
      tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
      tmp_fst <- neut_fst[neut_maf > tmp_bin[1] & neut_maf <= tmp_bin[2]]
    }
    if (method == "sliding") {
      tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
      tmp_fst <- neut_fst[neut_maf > tmp_bin[3] & neut_maf <= tmp_bin[4]]
    }
    
    tmp_percentiles <- ecdf(tmp_fst)
    tmp_p <- 1-tmp_percentiles(input_fst)
    return(tmp_p)
  }
  neut_pvals$L23.pval_split <- do.call("c",mclapply(1:nrow(info_neut), FUN = calc_neut_pval_L23))
  
  calc_neut_pval_D13 <- function(site) {
    input_maf <- info_neut$maf[site]
    input_fst <- info_neut$denmark_post_pre.fst[site]
    neut_maf <- info_neut$maf
    neut_fst <- info_neut$denmark_post_pre.fst
    
    if (method == "bins") {
      tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
      tmp_fst <- neut_fst[neut_maf > tmp_bin[1] & neut_maf <= tmp_bin[2]]
    }
    if (method == "sliding") {
      tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
      tmp_fst <- neut_fst[neut_maf > tmp_bin[3] & neut_maf <= tmp_bin[4]]
    }
    
    tmp_percentiles <- ecdf(tmp_fst)
    tmp_p <- 1-tmp_percentiles(input_fst)
    return(tmp_p)
  }
  neut_pvals$D13.pval_split <- do.call("c",mclapply(1:nrow(info_neut), FUN = calc_neut_pval_D13))
  
}
rm(list=ls(pattern="calc_neut_pval")); rm(list=ls(pattern="calc_pval_"))


########
rm(list=ls(pattern="calc_gwas_pval")); rm(list=ls(pattern="calc_neut_pval")); rm(list=ls(pattern="calc_exon_pval"))
########

# We'll end this document here and switch to an R markdown file to take a closer look at the results. 
save.image("./DATA_part1.RData")
