# We sequenced sets of multiplexed samples with separate libraries for infected and uninfected cells. 

##------------------------------------------
# Raw data were first mapped to the hg19 genome using cellranger. 
##------------------------------------------
# create genome to use with cellranger


# set sample or multiplexed library in place of NAME
INPUT=NAME 
# set directories where input files can be found (one director per lane of sequencing)
FILE_DIR1=loc1
FILE_DIR2=loc2
FILE_DIR3=loc3
FILE_DIR4=loc4

# run cellranger
srun -n 1 -c 12 cellranger count --localcores=12 --localmem=160 --id=${INPUT}_OUT 
--fastqs=$FILE_DIR1/$INPUT,$FILE_DIR2/$INPUT,$FILE_DIR3,$FILE_DIR4 
--sample=$INPUT --transcriptome=$GENOME  --expect-cells=12000 > logs/joblog-${INPUT}.out
#submit with >45G of memory

##------------------------------------------


##------------------------------------------
# Use souporcell to cluster cells by sample and then use those genotypes to assign samples to a previously genotyped individual via hierarchical clustering
##------------------------------------------

BAM_DIR=./cellranger/${INPUT}_OUT/outs
BARCODES_DIR=./cellranger/${BATCH}_OUT/outs/filtered_feature_bc_matrix
REF_DIR=./SOUPORCELL/references

cp ~/my_genomes/hg38_pestis/merged.fa ${BAM_DIR}/genome.fa
cp ${REF_DIR}/filtered_2p_1kgenomes_GRCh38.vcf ${BAM_DIR}

gunzip ${BARCODES_DIR}/barcodes.tsv.gz
singularity exec ~/Programs/SOUPORCELL/souporcell_latest.sif souporcell_pipeline.py \
  -i ${BAM_DIR}/possorted_genome_bam.bam \
                -b ${BARCODES_DIR}/barcodes.tsv \
                --fasta ${BAM_DIR}/genome.fa \
                --common_variants ${BAM_DIR}/filtered_2p_1kgenomes_GRCh38.vcf \
                --skip_remap SKIP_REMAP \
                -t 16 \
                -o ./v1_skip_remap_OUTS/${SAMPLE}_SOC_OUTS_bashTest \
                -k 6
bgzip ${BARCODES_DIR}/barcodes.tsv
# souporcell_pipeline.py is included in this github repository. souporcell_latest.sif is available along with the installation of souporcell
##------------------------------------------

# the outputs at this point can be read into R using Seurat. 
