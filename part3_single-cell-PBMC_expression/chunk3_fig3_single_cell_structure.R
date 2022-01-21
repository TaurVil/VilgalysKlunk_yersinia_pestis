##### Post processing, actual analysis and figures
# cd ./pestis_single_cell/
# salloc --mem=200G --partition=bigmem2 srun --pty bash
# module load R/4.0.4; module load hdf5; module load hdf5_hl; R
library(Seurat); library(SeuratDisk)
library(ggplot2); library(patchwork)
library(data.table); library(dplyr)
# BiocManager::install("scran")
library(scater); library(scran)
library(ggridges); library(ggthemes)
library(pbapply)
logCPM_filter = 0.5
options(future.globals.maxSize= 2*891289600)

##-----------------------------------------------
## get cells which we could assign to a known cell type
##-----------------------------------------------
logCPM_filter = 0.5
options(future.globals.maxSize= 2*891289600)

my_data <- readRDS("./pbmc_integrated_by_genotype-condition_n10_EU_hg38pestis.rds")
cell_types <- c("NK", "CD8 T", "CD4 T", "B", "Mono")
data <- subset(my_data, subset = predicted.celltype.l1 %in% cell_types)
rm(my_data)
##-----------------------------------------------


##-----------------------------------------------
## Plot UMAP with cell type clusters (fig 3d)
##-----------------------------------------------
fig3d <- DimPlot(data, reduction = "umap", group.by = "predicted.celltype.l1", label = T, label.size = 7, repel = TRUE) + ggtitle('') + theme(legend.position = 'none'); fig3d
DimPlot(data, reduction = "umap", group.by = "predicted.celltype.l2", label = T, label.size = 4, repel = TRUE) + theme(legend.position = 'none')
save(fig3d, file="~/Barreiro/pestis_aDNA/fig3D.rdata")
##-----------------------------------------------


##-----------------------------------------------
# UMAP colored by ERAP2
##-----------------------------------------------
DefaultAssay(data) <- 'RNA'

data1 <- subset(data, subset = stim == "UNSTIM" & genotype == 2)
data2 <- subset(data, subset = stim == "UNSTIM" & genotype == 0)

p1 <- FeaturePlot(data1, features= "ERAP2", reduction='umap', cols=c('lightgrey', 'darkred')) + ggtitle("unstim, genotype T/T (NMD)") & theme(plot.title = element_text(size = 10))
p2 <- FeaturePlot(data2, features= "ERAP2", reduction='umap', cols=c('lightgrey', 'darkred')) + ggtitle("unstim, genotype C/C (functional)") & theme(plot.title = element_text(size = 10))

p1 + p2 + plot_layout(ncol=2) +
  plot_layout(guides="collect") & theme(legend.position = 'bottom') &
  scale_color_gradient(low = "lightgrey", high = "darkred", limits = range(c(0, 1.6)))
##-----------------------------------------------

##-----------------------------------------------
# Nebulosa plot: fig 3E
##-----------------------------------------------
library(Nebulosa)
# get just the stimulated condition for the main figure, just the unstimulated condition for the SI figure. 
subset(data, subset = stim == "STIM") -> shifted_geno2

# move cells containing the protective genotype so they're plotted side by side
shifted_geno2[["umap"]]@cell.embeddings[shifted_geno2@meta.data$genotype == 0,1] <- shifted_geno2[["umap"]]@cell.embeddings[shifted_geno2@meta.data$genotype == 0,1] + 25
plot_density(shifted_geno2, "ERAP2") & theme(legend.position = 'bottom') &
  scale_color_gradient(low = "lightgrey", high = "darkred")

# separate out T cells and other cells so they have their own color gradient
Tcells <- subset(shifted_geno2, subset = predicted.celltype.l1 %in% c("CD8 T", "CD4 T"))
othercells <- subset(shifted_geno2, subset = predicted.celltype.l1 %in% c("NK", "B", "Mono"))

t <- plot_density(Tcells, "ERAP2") & theme(legend.position = 'bottom') &
  scale_color_gradient(low = "lightgrey", high = "darkred")
oth <- plot_density(othercells, "ERAP2") & theme(legend.position = 'bottom') &
  scale_color_gradient(low = "lightgrey", high = "darkred")

# scale colors
t$data$feature <- t$data$feature * (1/max(t$data$feature))
oth$data$feature <- oth$data$feature * (1/max(oth$data$feature))

# merge back together, and plot
tmp <- rbind(t$data, oth$data)
fig3e <- ggplot(data = tmp, aes(x = UMAP_1, y=UMAP_2, col = feature)) + 
  geom_point() + 
  scale_color_gradient(low = "lightgrey", high = "darkred") + 
  geom_vline(xintercept = 12, linetype = 'dashed', col='gray75') +
  theme_classic() & 
  theme(legend.position = 'none'); fig3e
save(fig3e, file="./fig3E.rdata")
##-----------------------------------------------

# repeat with unstimulated cells for Fig S5

