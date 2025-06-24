################################################################################
##        Title : RNAseq data analysis - mESCs and EBs (48h and 96h) ikba KO
##  Description : This script is for plotting PCA (all samples)
##   Researcher : Luis Galán
##       Author : M.Maqueda
##         Date : 21st July 2022
################################################################################


################################################################################
## 1. Load packages
################################################################################

require("openxlsx")
require("DESeq2")
require("dplyr")

################################################################################
## 2. Import data: normalized counts and samples information
################################################################################

# All samples normalized together (i.e. norm exprs matrix from GEO)
norm.counts <- read.table(file = "Downloads/GSE239563_RNAseq_All_samples_Norm_counts.txt", header=TRUE)
rownames(norm.counts) <- norm.counts$gene_ID
norm.counts <- norm.counts[,-1]

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

colnames(norm.counts) <- metadata$SampleName

################################################################################
## 3. Conduct PCA and plot first two principal components
################################################################################

# Manual computation with prcomp
mat <- norm.counts 
pca <-  prcomp(t(mat), scale=TRUE, center=TRUE)
PCs <- as.data.frame(pca$x)
PCs$Sample <- rownames(PCs)
PCs$Group <- metadata$Group
PCs$Timepoint <- metadata$Timepoint
Variance <- round(summary(pca)$importance[2,]*100, digits=1)

pdf(file = "PCA_PC1PC2_Complete.pdf", width=6, height=6)

ggplot(PCs, aes(PC1, PC2, fill=Group,label=Sample, col= Group, shape=Timepoint)) +
  geom_point(size=3,alpha=0.8,) +
  geom_text(aes(label=Sample),vjust= 0.5,hjust=-0.1,size=2, color="black")+
  xlab(paste0("PC1: ",Variance[1],"% variance")) +
  ylab(paste0("PC2: ",Variance[2],"% variance")) +
  scale_color_manual(values=c("lightskyblue3","darkred")) +
  scale_fill_manual(values=c("lightskyblue3","darkred"))+
  scale_shape_manual(values=c(24,22,21)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color="black"),
        legend.title = element_text(size = 12,face="italic"),
       legend.text = element_text(size = 12),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(size=14, hjust = 1),
        axis.text.y = element_text(size=14, hjust = 1)) +
  coord_fixed() +
  xlim(-110,210) +
  ggtitle("Transcriptome dynamics during mESC differentiation") 

dev.off()

