################################################################################
##        Title : RNAseq data analysis - mESCs and EBs (48h and 96h) - 
##                Time Course analysis
##  Description : This script is for plotting Xderm terms along time for KO and
##                WT conditions
##   Researcher : Luis Galán
##       Author : María Maqueda
##         Date : 13th Jul 2023
################################################################################


################################################################################
## 1. Load packages
################################################################################

require(DESeq2)
require(ggplot2)
require(GSVA)
require(org.Mm.eg.db)

################################################################################
## 2. Import data: normalized counts and scale data gene annotations and pathways of interest
################################################################################

# Normalized counts: all samples together since we are going to show all timepoints
# Norm expression mat can be obtained from GEO
#======

norm.mat <- read.table(file = "Downloads/GSE239563_RNAseq_All_samples_Norm_counts.txt", header=TRUE)
rownames(norm.mat) <- norm.mat$gene_ID
norm.mat <- norm.mat[,-1]

# Import metadata info to change colnames - can be downloaded from GitHub repo (or created)

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

metadata <- metadata %>%
  dplyr::mutate(SampleName = paste(Timepoint, paste(Group,Clone.n, sep=""), sep="_"))

colnames(norm.mat) <- metadata$SampleName

# Counts are already normalized, we do scaling
expr.mat.scaled <- t(apply(norm.mat, 1, scale))
colnames(expr.mat.scaled) <- colnames(norm.mat)

################################################################################
## 3. Import data: gene annotations and pathways of interest (GO BP terms)
################################################################################

# General annotations
#======
gene_annotations <- readRDS("/Users/mmaqueda/Documents/Genomes/MusMusculus/Ensembl_GRCm38_102/Mus_musculus.GRCm38.102.ONLY_GENES.RDS")

# colnames(gene_annotations)
# [1] "chrom"            "source"           "feature"          "start"            "end"              "score"            "strand"          
# [8] "frame"            "attribute"        "gene_id"          "gene_type"        "gene_name"        "Entrez"           "gene_biotype"    
# [15] "gene_description"

# nrow(gene_annotations)
# [1] 55487


# GO terms of interest: Obtain Entrez id genes vectors
#======
# Endoderm, mesoderm and ectoderm related signatures

GoI <- c("GO:0001706",  "GO:0001707",  "GO:0001705")
names(GoI) <- c("Endo_form",  "Meso_form","Ecto_form")

xderm.gs <- lapply(GoI, function(term) {
  genes <- AnnotationDbi::select(org.Mm.eg.db, keytype="GOALL", keys=term, columns="ENTREZID")
  return(genes$ENTREZID) })

names(xderm.gs) <- names(GoI)

xderm.gs.ensembl <- lapply(xderm.gs, function(genes) {
  a <- clusterProfiler::bitr(geneID= genes,toType = "ENSEMBL", fromType = "ENTREZID",OrgDb = org.Mm.eg.db, drop=TRUE)
  return(a$ENSEMBL)
})

################################################################################
## 4. Compute an score for each signature of interest. Option A: GSVA
################################################################################

# NOTE: Results are exactly the same if the expression matrix is scaled or not

# GSVA Z-score - default parameters
gsva.es <- gsva(expr = expr.mat.scaled, 
                method = "zscore", 
                gset.idx.list =xderm.gs.ensembl, 
                verbose=FALSE)

################################################################################
## 5. Generate plot showing including trend (loess curve)
################################################################################

toplot <- metadata
toplot <- cbind(toplot, t(gsva.es))

# head(toplot)
# Sample_id Timepoint Group Clone.n  SampleName  Endo_form    Meso_form   Ecto_form
# EBs_48h_KO1     D4_48   EBs_48h    KO       1 EBs_48h_KO1  0.6423981 -0.380328520  0.19178054
# EBs_96h_KO1     D4_96   EBs_96h    KO       1 EBs_96h_KO1 -0.7716598 -0.001268292  0.03351503
# EBs_48h_WT1    E14_48   EBs_48h    WT       1 EBs_48h_WT1 -1.1609376  0.426216018 -0.06669277
# EBs_96h_WT1    E14_96   EBs_96h    WT       1 EBs_96h_WT1  6.1928271 11.884678968  2.30698412
# EBs_48h_KO2     E2_48   EBs_48h    KO       2 EBs_48h_KO2  0.4574388 -1.165898175 -0.08825768
# EBs_96h_KO2     E2_96   EBs_96h    KO       2 EBs_96h_KO2 -1.2188309  0.401597330 -0.12415166


# Example with mesoderm formation

pdf(file = "Mesoderm_Formation_Score.all_timepoints.pdf",
    width=4, height=4)

ggplot(toplot, aes(x=as.numeric(Timepoint), y=Meso_form, colour=Group)) +
  scale_x_continuous(expand = expansion(mult=0.01, add=0)) +
  geom_jitter(
    aes(color = Group), 
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.1),
    size = 1.2) +
  geom_smooth(aes(color = Group),se = TRUE, method = "loess", alpha=0.1,level=0.9) +  #Or try 'lm'
  theme_classic() +
  theme(legend.position="right",
        plot.title = element_text(size=8),
        axis.text.x = element_text(size=4, face="bold", colour = "black"),
        axis.text.y = element_text(size=4, face="bold", colour = "black"),
        axis.ticks.y.left =  element_line(linewidth=0.2),
        axis.ticks.x.bottom =  element_blank(),
        axis.line.x.bottom=element_line(linewidth = 0.2),
        axis.line.y.left = element_line(linewidth = 0.2)) +
  ggtitle("Mesoderm formation - Score") +
  ylab("Z-score") +
  xlab("Timepoint")

dev.off()

