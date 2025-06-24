################################################################################
##        Title : RNAseq data analysis - mESCs and EBs (48h and 96h) ikba KO
##  Description : This script is for exploring the NfKB KEGG Pathway (mmu04064)
##                in terms of the expression of its genes components
##   Researcher : Luis Galán
##         Date : 27th March 2023
################################################################################


# Set up paths ------------------------------------------------------------

################################################################################
## 1. Set up paths
################################################################################

# Load packages -----------------------------------------------------------

################################################################################
## 2. Load packages
################################################################################

require(ggplot2)
require(KEGGREST)
require(dplyr)
require(reshape2)
require(ggpubr)

################################################################################
## 3. Import data: geneset from KEGG mmu04064 (NfkB signaling pathway)
################################################################################

# Get the list of numbers, gene symbols and gene description
names <- keggGet("mmu04064")[[1]]$GENE
# Entrez ids
nfkb.entrez <- names[seq(1,length(names),2)]
#Delete the gene number by deleting every other line
nfkb.symbol <-  names[seq(0,length(names),2)]
nfkb.symbol <- unlist(sapply(nfkb.symbol, function(n) unlist(strsplit(n, split=";"))[1]))

################################################################################
## 4. Import data: our RNAseq data: Normalized expression matrices (for plots)
## DEGs lists for GSEA
################################################################################

# Normalized counts: 
# Option A: normalized matrix - independently per timepoint as done in DEA - [original figure]
norm.counts <- list("mESCs" = readRDS(file = "RData/KO_vs_WT_comparisons/Normalized_counts_DataSubset_mEScs.rds"),
                    "EBs_48h" = readRDS(file = "RData/KO_vs_WT_comparisons/Normalized_counts_DataSubset_EBs_48h.rds"),
                    "EBs_96h" = readRDS(file = "RData/KO_vs_WT_comparisons/Normalized_counts_DataSubset_EBs_96h.rds"))

norm.counts <- lapply(norm.counts, function(dds) as.data.frame(assay(dds)))

# Option B: all samples normalized together (i.e. norm exprs matrix from GEO)
norm.counts <- read.table(file = "Downloads/GSE239563_RNAseq_All_samples_Norm_counts.txt", header=TRUE)
rownames(norm.counts) <- norm.counts$gene_ID
norm.counts <- norm.counts[,-1]

# Import metadata info to change colnames - can be downloaded from GitHub repo (or created)
                      
metadata <- read.table(file="Downloads/Metadata.txt", 
                       header=TRUE, sep="\t")
#metadata$Clone.n <- as.factor(metadata$Clone.n)  # Optional

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
                      
colnames(norm.counts) <- metadata$SampleName

  # Split them in a list: one element per timepoint
norm.counts <- list("mESCs" = norm.counts[,grep("mEScs", colnames(norm.counts))], 
                    "EBs_48h" = norm.counts[,grep("48h", colnames(norm.counts))],
                    "EBs_96h" = norm.counts[,grep("96h", colnames(norm.counts))])


# NOTE: Both options show the same overall conclusion -> no significant differences

################################################################################
## 5. Select the genes of interest in our normalized counts matrices
################################################################################

# Import gene annotations (complete set)
gene_annot <- readRDS("/Volumes/projectscomput/cancer/db_files/Genomes/Ensembl/mouse/mm10/release-102/Mus_musculus.GRCm38.102.ONLY_GENES.RDS")
# colnames(gene_annot)
# [1] "chrom"            "source"           "feature"          "start"            "end"              "score"            "strand"          
# [8] "frame"            "attribute"        "gene_id"          "gene_type"        "gene_name"        "Entrez"           "gene_biotype"    
# [15] "gene_description"

# NOTE: this is just to associate ENSEMBL ids. It can be done in different ways

norm.counts.nfkb <- lapply(norm.counts, function(mat){
  mat$GENE_ID <- rownames(mat)
  mat.match <- left_join(x=mat, y=gene_annot, by = c("GENE_ID" = "gene_id"))
  res <- mat.match[which(mat.match$Entrez %in% nfkb.entrez),]
  return(res)
})

# Put the data in long format: easy to plot with ggplots

norm.counts.nfkb.long <- lapply(norm.counts.nfkb, function(data){
  data <- data[,grep("KO|WT|GENE_ID", colnames(data))] # Select just the exprs data
  data <- melt(data,id.vars="GENE_ID")
  data$Condition <- ifelse((data$variable == "mEScs_WT1" | data$variable == "mEScs_WT2" | data$variable == "mEScs_WT3"|
                            data$variable == "EBs_48h_WT1" | data$variable == "EBs_48h_WT2" | data$variable == "EBs_48h_WT3"|
                            data$variable == "EBs_96h_WT1" | data$variable == "EBs_96h_WT2" | data$variable == "EBs_96h_WT3"),"WT", "KO")
  return(data)
})

################################################################################
## 6. Make some plots: Violin plot with boxplots
################################################################################

for(i in 1:length(norm.counts.nfkb.long))
{
  p <-ggplot(norm.counts.nfkb.long[[i]], aes(Condition, value, fill=Condition)) +
    geom_violin(alpha=0.5, lwd=0.4)  +
    geom_boxplot(width=0.2,lwd=0.3,outlier.size = 0.15) +
    geom_jitter(width=0.25, alpha=0.8, size=0.4, shape=21, color="black", stroke=0.1) +
    stat_compare_means(method = "wilcox.test", label.x = 0.6, label.y = max(norm.counts.nfkb.long[[1]]$value) +0.5, label="p.format", size=4) +
    theme_classic() +
    theme(legend.position="none",
          plot.title = element_text(size=12),
          axis.title.y=element_text(size=8, colour="black"),
          axis.text.x = element_text(size=8, face="bold", colour = "black"),
          axis.text.y = element_text(size=7, face="bold", colour = "black"),
          axis.ticks.y.left =  element_line(linewidth=0.2),
          axis.ticks.x.bottom =  element_blank(),
          axis.line.x.bottom=element_line(linewidth = 0.2),
          axis.line.y.left = element_line(linewidth = 0.2)) +
    ggtitle(names(norm.counts.nfkb.long)[i]) +
    ylab("Normalized Expression") +
    xlab("")
  
  pdf(file=paste("NfKB_GeneSet_ViolinPlots_", names(norm.counts.nfkb.long)[i],".pdf",sep=""), width=5, height=5)
  print(p)
  dev.off()
}

