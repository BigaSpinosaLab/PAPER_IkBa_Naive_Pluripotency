################################################################################
##        Title : RNAseq data analysis - mESCs and EBs (48h and 96h) ikba KO
##  Description : This script is for creating the raw and normalized (VST) 
##                DESeqDataSet objects and conducting DEA analysis KO vs WT for
##                Timepoints independently
##   Researcher : Luis Galán
##       Author : María Maqueda
##         Date : 21st July 2022
################################################################################

# NOTE: KO vs WT commparison is independently conducted per timepoint given that
# samples dispersion at 96h is notably different from mESCs, so it was decided to 
# do it separately.

# Set up paths ------------------------------------------------------------

################################################################################
## 1. Set up paths
################################################################################

.wd <- getwd()
.res <- file.path(.wd, "results")
.fig <- file.path(.res,"figures")
.tab <- file.path(.res,"tables")
.dat <- file.path(.wd,"data")
.RData <- file.path(.wd, "RData")

# Load packages -----------------------------------------------------------

################################################################################
## 2. Load packages
################################################################################

require("openxlsx")
require("DESeq2")
require("dplyr")
require("clusterProfiler")
require("ggvenn")
require("EnhancedVolcano")

# Import counts and samples data ------------------------------------------

################################################################################
## 3. Import data: counts, samples information and complete gene annotation
################################################################################

# 1.Import count data file: 
  # OPTION A: Expression matrix from featureCounts
counts <- read.delim(file=file.path(.dat,"IKB_EBs_gene_quantification.txt"), 
                     header=TRUE, sep="\t", skip=1)
counts <- counts[,-which(colnames(counts) %in% c("Chr","Start","End","Strand","Length"))]

# The Gene ID column is converted to the row names and removed
rownames(counts) <- counts$Geneid
counts <- counts[-(which(colnames(counts) %in% "Geneid"))]

# Simplify samples names to get the real ID
long.names <- colnames(counts)
long.names <- sapply(long.names, function(name) unlist(strsplit(name, split= "X.projects.cancer.IKB_EBs.STAR_align.BAM."))[2])
long.names <- sapply(long.names, function(name) paste(unlist(strsplit(name, split= "_"))[1:2], collapse="_"))

colnames(counts) <- long.names

  # OPTION B: Import raw counts table directly from GEO
counts <- read.table(file = "Downloads/GSE239563_RNAseq_All_samples_Raw_counts.txt",
                     header = TRUE)
rownames(counts) <- counts$gene_ID
counts <- counts[,-1]

# 2. Import metadata files (can be found at GitHub repo)
metadata <- read.table(file="Desktop/Scripts_IkBa_eLife_paper/RNAseq/Metadata.txt", 
                       header=TRUE, sep="\t")
metadata$Clone.n <- as.factor(metadata$Clone.n)

head(metadata)
#   Sample_id Timepoint Group Clone.n
# 1     D4_48   EBs_48h    KO       1
# 2     D4_96   EBs_96h    KO       1
# 3    E14_48   EBs_48h    WT       1
# 4    E14_96   EBs_96h    WT       1
# 5     E2_48   EBs_48h    KO       2
# 6     E2_96   EBs_96h    KO       2

# Optional (recommended): Construct a new sample name based on targets and store them in the counts matrix instead
metadata <- metadata %>%
  dplyr::mutate(SampleName = paste(Timepoint, paste(Group,Clone.n, sep=""), sep="_"))
colnames(counts) <- metadata$SampleName

# Define the reference levels
metadata$Group <- relevel(as.factor(metadata$Group),ref="WT")
metadata$Timepoint <- relevel(as.factor(metadata$Timepoint),ref="mEScs")

# 3. Import the complete gene annotation
gene_annotations <- readRDS("~/Documents/Genomes/MusMusculus/Ensembl_GRCm38_102/Mus_musculus.GRCm38.102.ONLY_GENES.RDS")

# colnames(gene_annotations)
# [1] "chrom"            "source"           "feature"          "start"            "end"              "score"            "strand"          
# [8] "frame"            "attribute"        "gene_id"          "gene_type"        "gene_name"        "Entrez"           "gene_biotype"    
# [15] "gene_description"

# nrow(gene_annotations)
# [1] 55487

# Sanity check: columns line up with your counts matrix
all(rownames(counts) == gene_annotations$gene_id)
# [1] TRUE

# If not, you would need to reorder
# gene_annotations <- gene_annotations[match(rownames(counts), gene_annotations$gene_id),]

################################################################################
## 2. Function to conduct DEA per subset (per timepoint, KO vs WT comparison)
## This function will return:
## a) RDS objects: DESeqDataSet already filtered by non (or very low) -expressed genes
## b) Excel files: DEGs (KO vs WT per timepoint)
## c) Volcano per KO vs WT comparison
################################################################################

conduct_DEA_on_subset <- function(complete_counts, 
                               complete_metadata,
                               complete_annotations,
                               timepoint,
                               path_to_store_RDS_objects,
                               path_to_store_excels,
                               path_to_store_figures)
{
  # 1. Subset the samples of interest to the selected timepoint
  subset <- which(complete_metadata$Timepoint %in% timepoint)
  metadata.subset <-  complete_metadata[subset,]
  counts.subset <- complete_counts[,subset] 
  
  # 2. Construct DESeqDataSet object with the subset
  data.dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts.subset, #counts.subset
                                             colData = metadata.subset, #metadata.subset
                                             design =  ~ Group, 
                                             tidy=FALSE)  # We've already add the rownames when importing data
  
  # 3. Remove low-count genes: those with less than 10 counts accross all samples 
  data.dds.red <- data.dds[rowSums(counts(data.dds)) >=10]
  keep <- which(rowSums(counts(data.dds)) >=10)  # This will be used for latter gene annotations
  
  # 4. Normalize data for later visualization
  vstData.dds<- vst(data.dds.red, blind=FALSE) 
  
  # 5. Conduct the test and get results for comparison of interest
  dea.dds <- DESeq(data.dds.red)
  results <- results(dea.dds,
                     name= "Group_KO_vs_WT",
                     alpha=0.05)
  
  # 6. LogFC shrinkage
  res.shr <- lfcShrink(dds = dea.dds,
                       coef= "Group_KO_vs_WT",
                       res = results,
                       type = "apeglm")
  
  # 7. Print summary of results in console to have a quick view
  summary(results)
  # 
  # 8. Store results (complete list) + lfcShrinkage + annotations
  KO_vs_WT <- cbind(as.data.frame(results(dea.dds, contrast=list(c("Group_KO_vs_WT")))),
                    "Shrunken_lFC" = as.data.frame(res.shr)$log2FoldChange,
                    "Shrunken_lFCSE" = as.data.frame(res.shr)$lfcSE,
                    gene_annotations[keep,c("chrom", "source", "start", "end",
                                            "gene_id", "gene_name", "Entrez", "gene_biotype")])
  # 
  # # 9. Print number of DEGs if applying abs(log2FoldChange) >1
  print(paste("UP:", nrow(as.data.frame(KO_vs_WT) %>%
                            filter(padj< 0.05, Shrunken_lFC >1)), sep=""))
  print(paste("DOWN:", nrow(as.data.frame(KO_vs_WT) %>%
                              filter(padj< 0.05, Shrunken_lFC < -1)), sep=""))
  
  ########
  # Volcano plotting: 
  #######

  p <- EnhancedVolcano(toptable = KO_vs_WT,
                       lab= KO_vs_WT$gene_name,
                       x ='Shrunken_lFC',
                       y= 'padj',
                       xlab = "Shrunken fold change",
                       ylab = "-Log10(adj pval)",
                       title = paste("KO_vs_WT: ",timepoint,sep=""),
                       pCutoff = 0.05,
                       FCcutoff = 1,
                       pointSize = 1.5,
                       labSize = 4.5,
                       legendLabels=c("Not sign.", "Shrunken FC", "adj pval","adj pval & shrunken FC"),
                       legendPosition = "right",
                       legendIconSize = 5.0,
                       legendLabSize = 10)
  
  pdf(file = file.path(path_to_store_figures,paste("Volcano_KO_vs_WT_",timepoint,".pdf",sep="")), width=14, height=9)
  print(p)
  dev.off()

  # ########
  # Store VST object and results table
  # #######
  
  saveRDS(vstData.dds,
          file=file.path(path_to_store_RDS_objects,paste("Normalized_counts_DataSubset_", timepoint,".rds", sep="")))
  # 
  hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12, fontName = "Arial Narrow", fgFill = "#4F80BD")
  write.xlsx(KO_vs_WT, file = file.path(path_to_store_excels, paste("DEA_results_KO_vs_WT_", timepoint, ".xlsx", sep="" )),
             rowNames=TRUE,
             headerStyle= hs)
  # 
  # # Return the dea results 
  return(KO_vs_WT)
}


# Conduct_Analysis --------------------------------------------------------

################################################################################
## 3. Conduct DEA analysis per timepoint of interest
################################################################################

# Subset: 0h

timepoint_subsets <- c("mEScs", "EBs_48h", "EBs_96h")

KO_vs_WT_timepoints <- lapply(timepoint_subsets, function(case){
  
  conduct_DEA_on_subset(complete_counts = counts, 
                        complete_metadata = metadata,
                        complete_annotations = gene_annotations,
                        timepoint = case,
                        path_to_store_RDS_objects = file.path(.RData,"KO_vs_WT_comparisons"),
                        path_to_store_excels = file.path(.tab,"KO_vs_WT_comparisons"),
                        path_to_store_figures = file.path(.fig,"KO_vs_WT_comparisons"))
})

names(KO_vs_WT_timepoints) <- timepoint_subsets

# Store all results in a RDS file
saveRDS(KO_vs_WT_timepoints, file=file.path(.RData,"KO_vs_WT_comparisons/KO_vs_WT_all_comparisons.RDS"))

################################################################################
## AUX. Venn Diagram to compare DEGs
################################################################################

A <- list("0h" = (KO_vs_WT_timepoints$mEScs %>% filter(padj <0.05, abs(Shrunken_lFC)>1))$gene_id,
          "48h" = (KO_vs_WT_timepoints$EBs_48h %>% filter(padj <0.05, abs(Shrunken_lFC)>1))$gene_id,
          "96h" = (KO_vs_WT_timepoints$EBs_96h %>% filter(padj <0.05, abs(Shrunken_lFC)>1))$gene_id)

pdf(file = file.path(.fig,"KO_vs_WT_comparisons/Venn_Diagram_DEGs_KO_vs_WT_different_timepoints.pdf"), width=6, height=6)
ggvenn(A,
       text_size = 3,
       set_name_size = 3)
dev.off()
