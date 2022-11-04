#!/bin/bash

module load java; module load samtools

g=${SLURM_ARRAY_TASK_ID}
f=`head -$g 00_mapDamage_bams.list | tail -1`;

samtools addreplacerg -r ID:$f -r PL:Illumina -r SM:$f -r LB:$f final_bams/$f.bam -o final_bams/rg.$f.bam
samtools index final_bams/rg.$f.bam

Programs/gatk-4.1.4.1/gatk HaplotypeCaller -ERC GVCF -I final_bams/rg.$f.bam  -R my_genomes/hg19/hg19.fa  -O gVCF/$f.g.vcf.gz -mbq 20 --sample-name $f

echo $f $g
