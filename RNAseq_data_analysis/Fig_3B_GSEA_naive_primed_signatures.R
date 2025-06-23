################################################################################
##        Title : RNAseq data analysis - mESCs and EBs (48h and 96h) ikba KO
##  Description : This script is for exploring the naïve/primed or ground state
##                of our mESCs cells. Based on publisched signatures of 
##                corresponding states.
##   Researcher : Luis Galán
##        Author: M.Maqueda
##         Date : 21st Nov 2022
################################################################################

## THIS SCRIPTS REPRODUCES FIGURE 3B

# Public dataset from: Ghimire et al. Nature Scientific Reports 2018
# https://www.nature.com/articles/s41598-018-24051-5

# Set up paths ------------------------------------------------------------

################################################################################
## 1. Set up paths
################################################################################

# Load packages -----------------------------------------------------------

################################################################################
## 2. Load packages
################################################################################

require(readxl)
require(fgsea)
require(ggplot2)
require(ggvenn)

################################################################################
## 3. Import data: naive, ground and primed signatures
################################################################################

# NOTE: GSEA is usually restricted to gene set of maximum 500 genes (original GSEA, min 15)

# Data from Ghimire et al. : Remove rows wo gene id
naive.signature <- read_excel(path = "data/Naive_signature/Naive_Ground_Primed_signatures.xls",
                              sheet = "Signature LS Up")
# 706 entries
naive.signature <- naive.signature[!is.na(naive.signature$geneid),]
naive.signature <- naive.signature[order(naive.signature$adj.P.Val),]
naive.signature <- naive.signature[!duplicated(naive.signature$geneid),]
# 509 entries

ground.signature <- read_excel(path = "data/Naive_signature/Naive_Ground_Primed_signatures.xls",
                               sheet = "Signature 2i Up")
# 997 entries
ground.signature <- ground.signature[!is.na(ground.signature$geneid),]
ground.signature <- ground.signature[order(ground.signature$adj.P.Val),]
ground.signature <- ground.signature[!duplicated(ground.signature$geneid),]
# 633 genes

# For directly using fGSEA
# Let's include all and increased the max size to 600
states.gs.list.all <- list("Naive" = unlist(naive.signature[,"geneid"]),
                           "Ground" = unlist(ground.signature[,"geneid"]))

################################################################################
## 4. Import data: our RNAseq data - Results from differential expression analysis
################################################################################

# DEGs: this RDS file includes a list with three elements, one per timepoint
# Each element is the result of the DEA (check corresponding script for generating it)
degs.all <- readRDS(file = "Documents/Projects/Luis_Galán/mESCs_EBs_Ikba/RNAseq_mESCs_EBs_Ikba/RData/KO_vs_WT_comparisons/KO_vs_WT_all_comparisons.rds")

# No need to select DEGs, for GSEA we need the complete list of genes tested in DESeq2
# We will sort them by shrunken logFC
degs.ordered <- lapply(degs.all,function(df) {
  df.order <-  df$Shrunken_lFC  
  #df.order <-  df$stat  # ¿Order by stat?
  names(df.order) <- df$gene_name # This is the gene symbol
  # Decreasing order
  df.order = sort(df.order, decreasing = TRUE)
  # Remove duplicates
  df.order <- df.order[!duplicated(names(df.order))]
  return(df.order)
})


################################################################################
## 5. Compute GSEA and compute some results
# NOTE: WE ONLY TEST THE mESCs (0h)
################################################################################

# WARNING: There are ties in the preranked stats (0.5%)

# About score type check: https://github.com/ctlab/fgsea/issues/87

set.seed(123)
fgseaRes <- fgseaMultilevel(pathways = states.gs.list, # This is for Hu et al. dataset
                                stats = degs.ordered$`0h`,
                                scoreType = "std",
                                eps=0,
                            maxSize = 600)

# Random walk plots
plotEnrichment(pathway = states.gs.list$Naive,
               stats = degs.ordered$`0h`) +
  labs(title="Naive")


plotEnrichment(pathway = states.gs.list$Ground,
               stats = degs.ordered$`0h`) +
  labs(title="Ground")
