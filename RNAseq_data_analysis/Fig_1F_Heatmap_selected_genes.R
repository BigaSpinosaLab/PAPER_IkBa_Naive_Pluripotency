################################################################################
##        Title : RNAseq data analysis - mESCs and EBs
##  Description : This script is for plotting  a heatmap of selected genes - 
#                 Reviewed version Fig 1F
##   Researcher : Luis Galán
##         Date : 6th Jun 2025
################################################################################

# NOTE: Since Fig1F only includes WT samples, we can normalize the data only 
# considering these samples. Method: VST implemented in DESeq2 package

################################################################################
## 1. Load packages
################################################################################

require(ggplot2)
require(ComplexHeatmap)
require(openxlsx)
require(DESeq2)

# Import counts and samples data ------------------------------------------

################################################################################
## 2. Import data: Raw values and annotations
################################################################################

# Raw counts to be normalized: Exprs matrix from GEO
raw.data <- read.table(file = "Downloads/GSE239563_RNAseq_All_samples_Raw_counts.txt",
                     header = TRUE)
rownames(raw.data) <- raw.data$gene_ID
raw.data <- raw.data[,-1]

# Corresponding metadata: this is just a dataframe with the basic description of the data
# It can be downloaded from Github repo
metadata <- read.table(file="Downloads/Metadata.txt", 
                       header=TRUE, sep="\t")
head(metadata)
#   Sample_id Timepoint Group Clone.n
# 1     D4_48   EBs_48h    KO       1
# 2     D4_96   EBs_96h    KO       1
# 3    E14_48   EBs_48h    WT       1
# 4    E14_96   EBs_96h    WT       1
# 5     E2_48   EBs_48h    KO       2
# 6     E2_96   EBs_96h    KO       2

# Complete gene annotation: we will use this just to select genes by name instead of ENSEMBL id
gene_annotations <- readRDS("~/Documents/Genomes/MusMusculus/Ensembl_GRCm38_102/Mus_musculus.GRCm38.102.ONLY_GENES.RDS")

################################################################################
## 3. Normalize raw data - all WT samples
################################################################################

# Just select the WT samples
metadata.wt <- metadata[which(metadata$Group %in% "WT"),]
raw.data.wt <- raw.data[ ,which(colnames(raw.data) %in% metadata.wt$Sample_id)]

# VST normalization NOT blind to experimental deseing
data.dds <- DESeq2::DESeqDataSetFromMatrix(countData = raw.data.wt, #counts.subset
                                           colData = metadata.wt, #metadata.subset
                                           design =  ~ Timepoint, 
                                           tidy=FALSE)  # We've already add the rownames when importing data
norm.data <- vst(data.dds, blind=FALSE) 
norm.mat <- assay(norm.data)

################################################################################
## 4. Define genes of interest per category
################################################################################

# Genes of interest split into categories
pluri.genes <- c("Pou5f1", "Dppa4", "Tbx3","Prdm14","Klf4","Zfp42", # Fig1F
                 "Sox2", "Klf2", "Gbx2", "Nanog", # Fig1F
                 "Gab1", "Esrrb", "Tfcp2l1") #additional form Fig2D
diff.genes <- c("Lef1", "Mesp1", "Gata6", "Lhx1", "Gata4","Wnt3a",
                "Fgf5", "Axin2", "T", "Mixl1", "Wnt3", "Eomes", "Foxa2", "Nes")
nfkb.genes <- c("Nfkbid", "Nfkbia", "Relb", "Rel", "Nfkb2", "Rela", "Nfkbie",
                "Nfkb1", "Nfkbiz", "Nfkbib")
polycomb.genes <- c("Ezh1", "Ezh2", "Eed", "Epop", "Jarid2", "Suz12", "Mtf2", "Aebp2")

# Put them in the same df
GoI <- data.frame("Gene" = c(pluri.genes, diff.genes, nfkb.genes, polycomb.genes),
                  "Category" = c(rep("Pluripotency", length(pluri.genes)),
                                 rep("Differentiation", length(diff.genes)),
                                 rep("NF-kB", length(nfkb.genes)),
                                 rep("Polycomb", length(polycomb.genes))))

# Check if all these symbols are included in the annotations list
setdiff(GoI$Gene, gene_annotations$gene_name)
#character(0) # Correct

# Add the ensembl ID in this dataframe
GoI.final <- merge.data.frame(x = GoI, y = gene_annotations, 
                 by.x = "Gene", by.y= "gene_name", all=FALSE)
GoI.final <- GoI.final[,c("Gene","Category","gene_id")]

################################################################################
## 5. Subset Genes of Interest from the complete matrix
################################################################################
                      
# Subset them from the matrix
mat.interest <- norm.mat[which(rownames(norm.mat) %in% GoI.final$gene_id),] 

# Sort them in the same order
mat.interest <- mat.interest[match(GoI.final$gene_id, rownames(mat.interest)),]

# Now Change row names to gene names (instead of Ensembl)
rownames(mat.interest) <- GoI.final$Gene

################################################################################
## 6. Construct and plot the heatmap
################################################################################

# Matrix should be genes x samples. Scale data
mat.interest.scaled <- t(scale(t(mat.interest)))

# Change the order of the columns
mat.interest.scaled <- mat.interest.scaled[,c("WT_1","WT_2","WT_3",
                                              "E14_48", "H7_48", "H9_48",
                                              "E14_96", "H7_96", "H9_96")]

timepoint = c(rep("mESCs",3), rep("EBs 48h",3), rep("EBs 96h",3))

# Value range in the matrix
quantile(mat.interest.scaled, c(0.01, 0.99))
# 1%       99% 
#   -1.607247  2.017222 

# Let's not consider outlier values for color grading
col_fun = circlize::colorRamp2(c(-1.6, 0, 1.8), c("skyblue", "white", "indianred3"))

htmp<- Heatmap(t(mat.interest.scaled), 
               name = "Scaled Norm Exprs",  
               cluster_rows = FALSE, 
               row_title_gp = gpar(fontsize = 10, fontface="bold"),
               row_gap = unit(0.2, "mm"),
               column_split = factor(GoI.final$Category),
               column_gap = unit(0.9, "mm"),
               cluster_columns = TRUE,
               show_column_dend = TRUE,
               col = col_fun,
               show_row_names = TRUE,
               row_names_side = "right", 
               column_names_gp=gpar(fontsize=6, fontface="bold"),
               column_names_centered = TRUE,
               column_names_rot = 45,
               show_column_names = TRUE,
               use_raster = TRUE,
               heatmap_legend_param = list(legend_direction = "horizontal", border = "black"),
               border_gp = gpar(col = "black", lty = 1),
               show_parent_dend_line = FALSE,
               raster_quality = 4)

pdf("Fig1F_Heatmap_FINAL_OPTION.pdf", width=6, height=3)
draw(htmp, heatmap_legend_side="bottom")
dev.off()
