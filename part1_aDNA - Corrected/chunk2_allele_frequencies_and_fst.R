# This file starts with the genotype likelihood files

######## Initiation, parameters, and functions ######## 
# libraries
library(data.table); library(ggplot2); library(patchwork)
# parameters
min_n=10
min_maf=0.0
drop_samples_missing=.5
# functions
#### These functions to estimate allele frequencies based on maximum likelihood is the primary difference in the correction
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
########


######## get time points to consider ########
design <- expand.grid(time = gl(3, 1, labels = c("pre", "post", "during")), pop = gl(2, 1, labels = c("denmark", "london"))); design <- design[-3,]
design2 <- expand.grid(time1 = gl(3, 1, labels = c("pre", "post", "during")), time2 = gl(3, 1, labels = c("pre", "post", "during")), pop = gl(2, 1, labels = c("denmark", "london")))
design2<-subset(design2, !(design2$time1 == design2$time2)); design2 <- design2[c(1,7:8,10),]
######## get sites to consider ########
info_gwas <- fread("~/GitHub/VilgalysKlunk_yersinia_pestis/part1_aDNA/DATA/genoliks/immune.frq")[,-c(3:4)]
info_neut <- fread("~/GitHub/VilgalysKlunk_yersinia_pestis/part1_aDNA/DATA/genoliks/neutral.frq")[,-c(3:4)]
info_exon <- fread("~/GitHub/VilgalysKlunk_yersinia_pestis/part1_aDNA/DATA/genoliks/exon.frq")[,-c(3:4)]
colnames(info_gwas) <- colnames(info_neut) <- colnames(info_exon) <- c("chr", "pos", "ref", "alt")

info_gwas$ref <- gsub(":.*", "",info_gwas$ref); info_gwas$alt <- gsub(":.*", "",info_gwas$alt)
info_neut$ref <- gsub(":.*", "",info_neut$ref); info_neut$alt <- gsub(":.*", "",info_neut$alt)
info_exon$ref <- gsub(":.*", "",info_exon$ref); info_exon$alt <- gsub(":.*", "",info_exon$alt)

paste(info_gwas$chr,info_gwas$pos,sep="_") -> info_gwas$site
paste(info_neut$chr,info_neut$pos,sep="_") -> info_neut$site
paste(info_exon$chr,info_exon$pos,sep="_") -> info_exon$site
########


######## get sample names per population and time point ########
L_snames_1 <- unique(c(
  read.delim(paste("./DATA/SampleNames/keep_london_pre_exon",sep=""), header=F)$V1, 
  read.delim(paste("./DATA/SampleNames/keep_london_pre_gwas",sep=""), header=F)$V1, 
  read.delim(paste("./DATA/SampleNames/keep_london_pre_neutral",sep=""), header=F)$V1
))
L_snames_2 <- unique(c(
  read.delim(paste("./DATA/SampleNames/keep_london_during_exon",sep=""), header=F)$V1, 
  read.delim(paste("./DATA/SampleNames/keep_london_during_gwas",sep=""), header=F)$V1, 
  read.delim(paste("./DATA/SampleNames/keep_london_during_neutral",sep=""), header=F)$V1
))
L_snames_3 <- unique(c(
  read.delim(paste("./DATA/SampleNames/keep_london_post_exon",sep=""), header=F)$V1, 
  read.delim(paste("./DATA/SampleNames/keep_london_post_gwas",sep=""), header=F)$V1, 
  read.delim(paste("./DATA/SampleNames/keep_london_post_neutral",sep=""), header=F)$V1
))
# same for denmark
D_snames_1 <- unique(c(
  read.delim(paste("./DATA/SampleNames/keep_denmark_pre_exon",sep=""), header=F)$V1, 
  read.delim(paste("./DATA/SampleNames/keep_denmark_pre_gwas",sep=""), header=F)$V1, 
  read.delim(paste("./DATA/SampleNames/keep_denmark_pre_neutral",sep=""), header=F)$V1
))
D_snames_3 <- unique(c(
  read.delim(paste("./DATA/SampleNames/keep_denmark_post_exon",sep=""), header=F)$V1, 
  read.delim(paste("./DATA/SampleNames/keep_denmark_post_gwas",sep=""), header=F)$V1, 
  read.delim(paste("./DATA/SampleNames/keep_denmark_post_neutral",sep=""), header=F)$V1
))
########


######## Get aDNA frequencies ######## 
# No filtering at this point
for (type in c("gwas", "neut", "exon")) {
  # Deal with discrepancy in intermediate file names
  if (type == "exon") {type2 <- "exons"} else {type2 <- type}
  
  for (pop in c("london", "denmark")) {
    if (pop == "london") {pop2 <- "L"}; if (pop == "denmark") {pop2 <- "D"}
    if (pop == "london") {times <- c("pre", "post", "during")}; if (pop == "denmark") {times <- c("pre", "post")}
    
    # create DATA, a p sites by 3*n samples matrix of genotype likelihoods, and the corresponding list of sample names 
    NAMES <- NULL
    DATA <- as.data.frame(matrix(nrow=nrow(get(paste0("info_",type))), ncol=1))
    for (time in times) {
      # read in sample names 
      tmp_names=read.delim(paste("./DATA/SampleNames/", pop, "_", time, "_",type2,".012.indv",sep=""), header=F)
      # get genolik file 
      d <- as.data.frame(fread(paste("./DATA/genoliks/genolik.",type2,"_",pop,"_", time, ".genolik",sep="")))
      # remove site. change missing data to NA. get number of samples (tmp_n).
      row.names(d) <- paste(d[,1], d[,2], sep="_"); d <- d[,-(1:2)] ; d[d == -1] <- NA; tmp_n <- ncol(d)/3
      # remove samples missing too much data
      tmp_missing=(colSums(is.na(d))/nrow(d))[seq(1,3*tmp_n,3)]
      keep <- which(tmp_missing <= drop_samples_missing); k2 <- as.data.frame(matrix(nrow=nrow(d), ncol=1)); for (k in keep) {k2 <- cbind(k2, d[,(k*3-2):(k*3)])}; d <- k2[,-1];
      tmp_names <- tmp_names$V1[keep]; n_pre <- ncol(d)/3
      
      # SAVE TO FILE FOR pop and type 
      NAMES <- c(NAMES, tmp_names)
      DATA <- cbind(DATA, d)
      rm(k2, k, tmp_missing,keep, tmp_names)
    }
    # Remove buffer row from data
    DATA <- DATA[,-1]
    
    # estimate allele frequencies for each time point, based on their assignments given above ([L,D]_snames_[1:3], which can be permuted)
    info <- get(paste0("info_", type))
    for (time in times) {
      if (time == "pre") {t2 <- "1"}; if (time == "during") {t2 <- "2"}; if (time == "post") {t2 <- "3"}
      # keep individuals given by "[L,D]_snames_[1:3]"
      keep <- get(paste0(pop2, "_snames_", t2)); keep <- which(NAMES %in% keep)
      # get data for just those individuals
      k2 <- as.data.frame(matrix(nrow=nrow(DATA), ncol=1)); for (k in keep) {k2 <- cbind(k2, DATA[,(k*3-2):(k*3)])}; tmp_data <- k2[,-1]; rm(k, keep, k2)
      n <- ncol(tmp_data)/3
      
      # save number of samples called and the genotype likelihoods
      info[[paste(pop, time, "called", sep=".")]] <- (ncol(tmp_data)-rowSums(is.na(tmp_data)))/3
      # info[[paste(pop,time,"alternate.GL",sep=".")]] <- rowMeans(tmp_data[,seq(3,n*3,3)],na.rm=T) + rowMeans(tmp_data[,seq(2,n*3,3)],na.rm=T)/2
      info[[paste(pop,time,"alternate.ML",sep=".")]] <- estimate_af_ml(tmp_data)
    }; rm(time)
    assign(x = paste0("info_", type), value = info); rm(info, pop2)
  }; rm(type2)
}; rm(type, pop, times)
rm(n, n_pre, NAMES, t2)
rm(L_snames_1, L_snames_2, D_snames_1, D_snames_3, L_snames_3, tmp_data,d,DATA)
rm(estimate_af_ml, likelihood)
########  

########  Merge into single "info" file ########  
info_gwas$type <- "gwas"; info_exon$type <- "exon"; info_neut$type <- "neut"
info <- rbind(info_gwas, info_exon, info_neut)
rm(info_gwas, info_exon, info_neut)
info <- info[!duplicated(info$site),]
########  

####### filter for sites with 10 samples per population in each time point #######
attach(info)
info <- subset(info, apply(X = cbind(london.during.called, london.pre.called, london.post.called, denmark.pre.called, denmark.post.called), 1, min) >= min_n)
detach(info)
#######

####### calculate alternate and minor allele frequency #######
info$london.alternate.ML <- rowMeans(cbind(info$london.pre.alternate.ML, info$london.post.alternate.ML, info$london.during.alternate.ML))
info$denmark.alternate.ML <- rowMeans(cbind(info$denmark.pre.alternate.ML, info$denmark.post.alternate.ML))
info$alternate.ML <- rowMeans(cbind(info$london.alternate.ML, info$denmark.alternate.ML))
info$maf.ML <- apply(cbind(info$alternate.ML, 1-info$alternate.ML), 1, min)
#######

####### calculate fst #######
# calculate expected heterozygosity
for (pop in c("london", "denmark")) {
  if (pop == "london") {times <- c("pre", "post", "during")}
  if (pop == "denmark") {times <- c("pre", "post")}
  
  for (time in times) {
    # expected heterozygosity
    info[[paste(pop,time,"ehet.ML",sep=".")]] <- 2*info[[paste(pop,time,"alternate.ML",sep=".")]]*(1-info[[paste(pop,time,"alternate.ML",sep=".")]])
  }
  rm(time, times)
}; rm(pop)
# pairwise means, pairwise expectations, and Fst
for (i in 1:3) {
  n2=paste(design2$pop[i],design2$time1[i],design2$time2[i],sep=".")
  info[[paste(n2,"mean.ML",sep=".")]] <- rowMeans(cbind(info[[paste(design2$pop[i],".",design2$time1[i],".alternate.ML",sep="")]], info[[paste(design2$pop[i],".",design2$time2[i],".alternate.ML",sep="")]]))
  info[[paste(n2,"ehet.ML",sep=".")]] <- 2*info[[paste(n2,"mean.ML",sep=".")]]*(1-info[[paste(n2,"mean.ML",sep=".")]])
  info[[paste(n2,"fst.ML",sep=".")]] <- (info[[paste(n2,"ehet.ML",sep=".")]] - (info[[paste(design2$pop[i],".",design2$time1[i],".ehet.ML",sep="")]] + info[[paste(design2$pop[i],".",design2$time2[i],".ehet.ML",sep="")]])/2)/info[[paste(n2,"ehet.ML",sep=".")]]
}; rm(i,n2)
########

######## Get difference between time points ########
info$delta_L13.ML <- info$london.post.alternate.ML-info$london.pre.alternate.ML
info$delta_L12.ML <- info$london.during.alternate.ML-info$london.pre.alternate.ML
info$delta_D13.ML <- info$denmark.post.alternate.ML-info$denmark.pre.alternate.ML
########



######## Assign rsid for each variant ########
# Download rsids from NCBI's dbSNP:  https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz & https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi
# bcftools annotate --rename-chrs rename_chroms.txt --threads 8 -o ncbi.renamedchrs.vcf.gz -Oz GCF_000001405.25.gz #rename to chr# rather than ncbi chromosome ids
# zcat ncbi.renamedchrs.vcf.gz | grep -v '^#' | grep '1000Genomes' | awk -F "\t" -v OFS="\t" '{print $1, $2, $2, $3, $4, $5}' > ncbi_1000G.bed
# bedtools intersect -a ncbi_1000G.bed -b targetted_regions.bed -wa -wb > captured.ncbi1k_rsid.txt

rsid <- fread("./DATA/captured.ncbi1k_rsid.txt")
rsid <- rsid[!duplicated(paste(rsid$V3, rsid$V1, rsid$V4, sep="_")),]
rsid$site <- paste(rsid$V1, "_", rsid$V3, sep="")
colnames(rsid)[1:6] <- c("chr", "pos", "pos2", "rsid", "ref", "alt")  # columns 7-9 carry the targeted regions and don't matter. 

info$rsid <- info$rsid_number <- NA
info$rsid_number <- as.numeric(info$rsid_number)
rsid <- as.data.frame(rsid)

# info$rsid <- as.character(info$rsid)
for (i in 1:nrow(info)) {
  # get rsid at that site (i.e., same start site)
  tmp <- rsid[rsid$site == info$site[i],]
  # get rsid that say the same reference allele
  if (nrow(tmp) > 0) {
    tmp <- tmp[tmp$ref == info$ref[i],]
  }
  # get rsid with the same alternate allele
  if (nrow(tmp) > 0) {
    # using strsplit is necessary rather than 'like' to prevent partial matches. Like would allow "C" to match "ACG" like it matches "A,C"
    keep <- NULL; for (j in 1:nrow(tmp)) {
      tmp2 <- strsplit(tmp$alt[j],",")[[1]]
      if (info$alt[i] %in% tmp2) {keep <- cbind(keep,tmp$rsid[j])}
    }; rm(j, tmp2)
    
    length(keep) -> info$rsid_number[i]
    tmp <- tmp[tmp$rsid %in% keep,]
  }
  if (nrow(tmp) > 0) {
    info$rsid[i] <- paste(as.character(tmp$rsid), collapse=", ")
  }
  if (i/1000 == round(i/1000)) {print(i)}
}; rm(keep, rsid, tmp, tmp_n,i)
########


# sites only filtered for 10 individuals and known rsID. No maf filter. 
info <- as.data.frame(info[info$rsid_number == 1,])
info <- info[,-which(colnames(info) == "rsid_number")]
# save.image("./fst_estimates.RData")


# To obtain the 5,884 variants mentioned in the text, run with min_n = 0. 