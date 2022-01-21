In this section, we describe the identificaiton of genes year our candidates for positive selection which response with differential expression upon _Yersinia pestis_ infection and how our candidate variants affect the expression of those genes. 

**chunk 1**: Read in expression data and model the response to infection. Here we start with the raw count data (also available on GEO) which we filter for genes expressed at a mean of at least 1 cpm in either the infected or non-infected condition. We then text for differences in expression between infected and uninfected cells, controlling for individual identity. 

**chunk  2**: We combine the previous data frame with genotypes for these samples to test whether variants under positive selection have an effect on gene expression in normal conditions and when responding to pestis. We also ask whether these two responses differ (which they don't for the candidate loci we focused on). 
