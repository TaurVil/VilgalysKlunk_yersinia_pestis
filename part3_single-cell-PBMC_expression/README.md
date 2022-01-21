In this section, we further explore the effect of our candidate variant on chromosome 5 on ERAP2 expression, using single-cell RNA sequencing. Using genotypes previously called for these samples (Randolph et al. 2021, Science), we focus on 5 individuals who are homozygous for the protective allele and 5 individuals homozygous for the susceptible variant. 

**chunk 1**: Get matrix of single-cell counts the fastq files, using cellranger and souporcell.

**chunk 2**: Using Seurat in R, merge matrixes across individuals, assign cell types, and integrate. Then calculate pseudobulk and save that output file (`./DATA/pbmc_n10_pseudobulk_means.RData`). Intermediate files including cell-level data are too large to store on GitHub (1-2G), but available upon request (please contact Tauras).  

**chunk 3**: Output Nebulosa plot (fig3) from Seurat objects created in chunk 2. 

**chunk 4**: Analyze pseudobulk gene expression returned from chunk 2, for effects of stimulation and genotype. 



