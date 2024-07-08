# Identifying selection using changes in allele frequency during the black plague

## Genotype calling
Start with the bam files provided by Jennifer Klunk. 

## Uniquely mapped reads
Files were originally mapped to only the regions of the genome targetted for enrichment. This produced weird patterns of heterozygosity due to calling as genetic variants regions which would have actually mapped to another part of the genome. To correct this, we remap reads to the entire genome, and then retain reads which mapped uniquely within the target region. 

First, we'll get the list of bam files to consider. We'll extract just the name of the sample, which we'll use as the sample name from here on. This file is included here `00_bam.list`. 

```console 
# Get list of individuals and genomes
ls */*bam > 00_bam.list; sed -i 's/.min24.MQ30.merged.RG.bam//g' 00_bam.list; sed -i 's/_20.*//g' 00_bam.list; sed -i 's/Final//g' 00_bam.list; sed -i 's/_redesigned//g' 00_bam.list; sed -i 's/bams\///g' 00_bam.list; sed -i 's/_seqs//g' 00_bam.list
```

Next map each bam file to the full human genome using BWA bam2bam mapping after removing singleton reads whose pairmate failed to map. At this point the files are stored in `./bams/`. 

For each file, we can check how many reads mapped uniquely out of the klunk bam file. We can also use samtools flagstat to get the proportion of reads that mapped to the full genome period, but that's done below. the "XT:A:U" tag is used by bwa to flag uniquely mapped reads

`samtools view ./klunk_bams/$f*.bam | wc -l`
`samtools view ./bams/$f.remapped.sort.bam | wc -l`
`samtools view ./bams/$f.remapped.sort.bam | grep "XT:A:U" | wc -l` 

```console
# index genome
# indexing needs to be done with same bwa version as mapping: ../Programs/bwa_bam2bam/network-aware-bwa-master/bwa index ../my_genomes/hg19/hg19.fa

# Get working path for bam2bam mapping
	module load samtools; module load bedtools; module load java
	export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:../Programs/zeromq/libzmq-master/build
	export LIBRARY_PATH=$LIBRARY_PATH:../Programs/zeromq/libzmq-master/build/lib:../Programs/luarocks/lib/lua/5.3/
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../Programs/zeromq/libzmq-master/build/lib:../Programs/luarocks/lib/lua/5.3/
	export CPATH=$CPATH:../Programs/zeromq/libzmq-master/include/
	module load gcc; module load cmake/3.6.2; module load lua; module load bwa

# For each bam file, sort, convert to fastq, remove any PE read information (failure to collapse reads), convert to bam
# there are a total of 1003 files to run through
# the fastest way would be to do this all in parallele on the cluster, but due to cueuing times it's actually faster to do a mix of interactive jobs and queued jobs unless you want to run it overnight or something
# sbatch --array=1-503%225 --mem=8G --account=pi-lbarreiro --partition=lbarreiro run.01.bwa_map.sh
for index in `seq 335 336`; do f=`head -$index 00_bam.list | tail -1`;
	# sort original bam by queryname, and convert to fastq files
	samtools sort -n klunk_bams/$f*min24*.bam -o $f.sortname.bam
	bedtools bamtofastq -i $f.sortname.bam -fq $f.R.fq
	bedtools bamtofastq -i $f.sortname.bam -fq $f.R1.fq -fq2 $f.R2.fq

	# remove PE reads from SE dataset 
	awk 'NR%4==1 {print substr($1,2)}' $f.R1.fq | sed 's/\/1//g'> $f.names #extract read name from fastq
	../Programs/bbmap/filterbyname.sh in=$f.R.fq names=$f.names out=$f.remaining.fq exclude qin=33 overwrite=T
	
	# convert fq to bam; merge the two bams
	java -jar ../Programs/picard.jar FastqToSam FASTQ=$f.R1.fq FASTQ2=$f.R2.fq OUTPUT=tmp5.$f.umap_fq.bam READ_GROUP_NAME=$f SAMPLE_NAME=$f LIBRARY_NAME=ILL
	java -jar /project2/lbarreiro/users/tauras/Programs/picard.jar FastqToSam FASTQ=$f.remaining.fq OUTPUT=tmp5.$f.umap_fq2.bam READ_GROUP_NAME=$f SAMPLE_NAME=$f LIBRARY_NAME=ILL
	samtools merge tmp6.$f.tomap.bam tmp5.$f.umap_fq2.bam tmp5.$f.umap_fq.bam

	# map 
	../Programs/bwa_bam2bam/network-aware-bwa-master/bwa bam2bam -g ../my_genomes/hg19/hg19.fa -n 0.01 -o 2 -l 16500 -f bams/$f.remapped.bam tmp6.$f.tomap.bam
	samtools sort -o bams/$f.remapped.sort.bam bams/$f.remapped.bam
	samtools index bams/$f.remapped.sort.bam; samtools flagstat bams/$f.remapped.sort.bam
	rm bams/$f.remapped.bam; rm tmp*$f.*bam; rm $f.*fq; rm $f.names; rm $f.sortname.bam
	echo $index
done

# let's get information such as the number of reads per sample, and the number of uniquely mapped reads
touch 01.reads_per_bam.txt; for index in `seq 1 1003`; do f=`head -$index 00_bam.list | tail -1`; samtools view ./bams/$f.remapped.sort.bam | wc -l >> 01.reads_per_bam.txt; done
touch 01.unique_reads_per_bam.txt; for index in `seq 1 1003`; do f=`head -$index 00_bam.list | tail -1`; samtools view ./bams/$f.remapped.sort.bam | grep "XT:A:U" | wc -l >> 01.unique_reads_per_bam.txt; done

# let's also get a single bam per sample 
sed -e 's/_.*//g' 00_bam.list | sed -e 's/a$//g' | sed -e 's/R$//g' |sort | uniq > 00_sample.list
module load samtools
for index in `seq 1 363`; do f=`head -$index 00_sample.list | tail -1`; samtools merge ./bams/merged.$f.bam ./bams/$f[a,R,_]*.remapped.sort.bam ; samtools index ./bams/merged.$f.bam; echo $f; done
## check that we don't have anything greater than 3: for index in `seq 1 363`; do f=`head -$index 00_sample.list | tail -1`; ls ./bams/$f[a,R,_]*.remapped.sort.bam | wc -l | grep -v 3 | grep -v 2; done
```

# Clean up aDNA damage

```console
# For each sample, run mapDamage (n=364)
sbatch --array=111-120 --mem=2G run.02.mapdamage.sh
# Again, we'll mix this between batch submitted and interactive jobs
	module load gcc; module load python/cpython-3.8.5; module load R
	export PATH=$PATH:$HOME/.local/bin
	cd ./bams/
	for index in `seq 111 120`; do f=`head -$index ../00_sample.list | tail -1`; mapDamage -i merged.$f.bam -r ../my_genomes/hg19/hg19.fa --rescale; echo $index $f; done 

# move the bam files to a new directory rather than the sub-directories created by mapDamage 
mkdir bams_mapDamage; mv ./bams/merged.*.mapDamage/*bam ./bams_mapDamage/
ls bams_mapDamage/ | sed -e 's/merged.//g' | sed -e 's/.rescaled.bam//g' > 00_mapDamage_bams.list

# Trim each read
# -L and -R commands indicate to trim the first/last 4bp of each read
mkdir bams_trim4
for f in `cat 00_mapDamage_bams.list`; do ../Programs/bamUtil/bam trimBam bams_mapDamage/merged.$f.rescaled.bam bams_trim4/$f.bam -L 4 -R 4; echo $f; done

```


# Get genotype calls
Now we need to call genotypes for each dataset, starting with the gVCF files then assembling the combined call sets. 
```console
mkdir gVCF
## make gVCF file for each individual
	sbatch --array=341-361%80 --mem=16G --account=pi-lbarreiro --partition=lbarreiro ./run.03.gVCF_clipped.sh
	sbatch --array=201-300%100 --mem=16G ./run.03.gVCF_clipped.sh
  
ls gVCF/*gz > 02_cohortmap; awk 'BEGIN { FS="/t"; OFS="/t" } { print $1 $1}'  02_cohortmap > tmp4 ; sed 's/gzgVCF/gz \t gVCF/g' tmp4 > tmp2; sed 's/ gVCF/gVCF/g' tmp2 > 02_cohortmap; sed -i 's/^gVCF\///g' 02_cohortmap; sed -i 's/.g.vcf.gz //g' 02_cohortmap; rm tmp2; rm tmp4

ls gVCF/*gz > 03_gvcf.list; sbatch --array=1-22 --mem=16G --account=pi-lbarreiro --partition=lbarreiro run.04.merge_gVCF.sh

mv vcf/*g.vcf.gz* gVCF/ # remove *.g.vcf.gz files
module load bcftools; bcftools concat ./vcf/chr*.vcf.gz -a > merged.vcf.gz
module load htslib; tabix merged.vcf.gz

module load vcftools
vcftools --gzvcf merged.vcf.gz --mac 1 --max-alleles 2 --remove-indels --minQ 30 --freq --out merged
vcftools --gzvcf merged.vcf.gz --mac 1 --max-alleles 2 --remove-indels --minQ 30 --singletons --out merged

vcftools --gzvcf merged.vcf.gz --mac 1 --max-alleles 2 --remove-indels --bed ./neutral.bed --minQ 30 --recode --out ./vcf/neutral_v1
vcftools --vcf ./vcf/neutral_v1.recode.vcf --exclude-positions merged.singletons --recode --out ./vcf/neutral

vcftools --gzvcf merged.vcf.gz --mac 1 --max-alleles 2 --remove-indels --bed ./exon.bed --minQ 30 --recode --out ./vcf/exon_v1
vcftools --vcf ./vcf/exon_v1.recode.vcf --exclude-positions merged.singletons --recode --out ./vcf/exon

vcftools --gzvcf merged.vcf.gz --mac 1 --max-alleles 2 --remove-indels --bed ./immune.bed --minQ 30 --recode --out ./vcf/immune_v1
vcftools --vcf ./vcf/immune_v1.recode.vcf --exclude-positions merged.singletons --recode --out ./vcf/immune


```

# Filter genotype calls for target regions
```console

bcftools annotate -x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ,^FORMAT/PL ./vcf/immune.recode.vcf | bcftools +setGT -- -ta -nu > joint.immune.forgenolik.vcf
bcftools annotate -x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ,^FORMAT/PL ./vcf/exon.recode.vcf | bcftools +setGT -- -ta -nu > joint.exon.forgenolik.vcf
bcftools annotate -x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ,^FORMAT/PL ./vcf/neutral.recode.vcf | bcftools +setGT -- -ta -nu > joint.neutral.forgenolik.vcf
```

# Get genotype likelihoods for each subset
```console
module load vcftools
## Exons
vcftools --vcf joint.exon.forgenolik.vcf --keep 05.pre_exon_london.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out london_pre_exons
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../Programs/LCLAE/filtbaboon1b 38 > genolik.exons_london_pre.genolik
vcftools --vcf joint.exon.forgenolik.vcf --keep 05.post_exon_london.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out london_post_exons
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../Programs/LCLAE/filtbaboon1b 63 > genolik.exons_london_post.genolik
vcftools --vcf joint.exon.forgenolik.vcf --keep 05.BD_exon_london.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out london_during_exons
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../Programs/LCLAE/filtbaboon1b 41 > genolik.exons_london_during.genolik

	vcftools --vcf joint.exon.forgenolik.vcf --keep 05.pre_exon_denmark.txt --recode --out temp2; vcftools --vcf temp2.recode.vcf --012 --out denmark_pre_exons 
	sed '/^#/d' temp2.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../Programs/LCLAE/filtbaboon1b 42 > genolik.exons_denmark_pre.genolik
	vcftools --vcf joint.exon.forgenolik.vcf --keep 05.post_exon_denmark.txt --recode --out temp2; vcftools --vcf temp2.recode.vcf --012 --out denmark_post_exons
	sed '/^#/d' temp2.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../Programs/LCLAE/filtbaboon1b 58 > genolik.exons_denmark_post.genolik
	vcftools --vcf joint.exon.forgenolik.vcf --keep 05.BD_exon_denmark.txt --recode --out temp2; vcftools --vcf temp2.recode.vcf --012 --out  denmark_during_exons
	sed '/^#/d' temp2.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../Programs/LCLAE/filtbaboon1b 24 > genolik.exons_denmark_during.genolik


## Immune
vcftools --vcf joint.immune.forgenolik.vcf --keep 00_London_pre_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out london_pre_immune; sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../Programs/LCLAE/filtbaboon1b 65 > genolik.gwas_london_pre.genolik
vcftools --vcf joint.immune.forgenolik.vcf --keep 00_London_post_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out london_post_immune; sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../Programs/LCLAE/filtbaboon1b 98 > genolik.gwas_london_post.genolik
vcftools --vcf joint.immune.forgenolik.vcf --keep 00_London_during_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out london_during_immune; sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../Programs/LCLAE/filtbaboon1b 62 > genolik.gwas_london_during.genolik

vcftools --vcf joint.immune.forgenolik.vcf --keep 00_Denmark_pre_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out denmark_pre_immune
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../Programs/LCLAE/filtbaboon1b 42 > genolik.gwas_denmark_pre.genolik

vcftools --vcf joint.immune.forgenolik.vcf --keep 00_Denmark_post_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out denmark_post_immune
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../Programs/LCLAE/filtbaboon1b 57 > genolik.gwas_denmark_post.genolik

vcftools --vcf joint.immune.forgenolik.vcf --keep 00_Denmark_during_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out denmark_during_immune
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../Programs/LCLAE/filtbaboon1b 24 > genolik.gwas_denmark_during.genolik

## neutral
vcftools --vcf joint.neutral.forgenolik.vcf --keep 00_London_pre_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out london_pre_neutral; sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../Programs/LCLAE/filtbaboon1b 65 > genolik.neutral_london_pre.genolik
vcftools --vcf joint.neutral.forgenolik.vcf --keep 00_London_post_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out london_post_neutral; sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../Programs/LCLAE/filtbaboon1b 98 > genolik.neutral_london_post.genolik
vcftools --vcf joint.neutral.forgenolik.vcf --keep 00_London_during_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out london_during_neutral; sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../Programs/LCLAE/filtbaboon1b 62 > genolik.neutral_london_during.genolik

vcftools --vcf joint.neutral.forgenolik.vcf --keep 00_Denmark_pre_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out denmark_pre_neutral
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../Programs/LCLAE/filtbaboon1b 42 > genolik.neutral_denmark_pre.genolik

vcftools --vcf joint.neutral.forgenolik.vcf --keep 00_Denmark_post_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out denmark_post_neutral
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../LCLAE/filtbaboon1b 57 > genolik.neutral_denmark_post.genolik

vcftools --vcf joint.neutral.forgenolik.vcf --keep 00_Denmark_during_neutral.txt --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out denmark_during_neutral
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ../Programs/LCLAE/filtbaboon1b 24 > genolik.neutral_denmark_during.genolik


## get site information and alternate variants for table
vcftools --vcf joint.neutral.forgenolik.vcf --freq --out neutral
vcftools --vcf joint.exon.forgenolik.vcf --freq --out exon
vcftools --vcf joint.immune.forgenolik.vcf --freq --out immune

```
Exons, London: pre=38, post=63, BD=41. Denmark: pre=42, post=58, BD=24 (n=264)

Immune & Neutral, London: pre=65, post=98 & 100, BD=64 & 63. Denmark: pre=42, post=57, BD=24 (n=350 & 352)

The resulting genolikelihood files feed into part2. 

