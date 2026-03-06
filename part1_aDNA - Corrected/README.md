In this section, we describe: 

**chunk 1**: We map targetted enrichment for ancient genomes, and then call genotypes using GATK. Genotype likelihoods are extracted using LCLAE (Wall et al. 2016, Molecular Ecology) which are stored in `DATA` and used as the input for chunk 2. 

**chunk  2**: In R, we use the genotype likelihoods from earlier to calculate allele frequencies in each population and time point. We then calculate Fst, and a p-value for candidate variants based upon the distribution of Fst in neutral regions of the genome. We save an R environment at the end of this step labeled `DATA_part1.RData`. 

**chunk 3**: Using the results from chunk 2, identify an enrichment of positive selection in immune-related regions and identify candidate variants for positive selection (Figure 2, Table S1). 

These files differ from the originally published version by using a maximum likelihood based estimate of allele frequency, rather than the sum of genotype likelihoods. This avoids a bias where low frequency variants were overestimated. 