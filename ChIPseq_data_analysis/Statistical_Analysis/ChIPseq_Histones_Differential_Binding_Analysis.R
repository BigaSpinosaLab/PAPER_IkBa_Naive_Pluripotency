################################################################################
##        Title : ChIPseq H3K27ac/H3K4me1/H3K4me3/H3K27me3 histone marks on mESCs 
##                - IbBa KO and WT conditions
##  Description : This script is for conducting differential binding analysis with
##                diffBind between identified peaks in KO and WT. Annotation with
##                ChIPseeker is also included
##   Researcher : Luis Galan
##       Author : María Maqueda
##         Date : 7th July 2022 / 23rd May 2023
################################################################################

# NOTE 1: Same script is used for the four histone marks although there are slight
# differences upon a narrow mark (H3K4me3) and broad mark (the rest).
# Code is referred to H3K4me1 as an example. The different parameters for H3K4me3
# are specified. Partial results are also shown for better reading.

# NOTE 2: For differential peaks annotation, release 102 from Ensembl is used (mm10)

# NOTE 3: For complete reproducibility of this script, access to corresponding 
# BAM files is necessary (including unique mapped reads, marked duplicates and 
# already filtered for blacklisted regions). On the other hand, corresponding
# called peaks are available at GEO repo

################################################################################
## 1. Load packages
################################################################################

require(DiffBind)
require(openxlsx)
require(rtracklayer)
require(profileplyr)

require(ChIPseeker)  # For annotation of diff binded peaks
require(GenomicFeatures)  # For creating a TxDb object from gtf file
require(org.Mm.eg.db)

require(dplyr)
require(EnhancedVolcano)
require(grid)


########################
# 2. Import Data: samples sheet referring to BAM IP and input files and called peaks
########################

# Example for H3K4me1. Same for the rest
samples <- read.xlsx(xlsxFile = "Desktop/Scripts_IkBa_eLife_paper/ChIPseq/Statistical_Analysis/ChIPseq_Histones_SampleSheets_DBA.xlsx",
                     sheet= "H3K4me1")
#samples[1,]
#   SampleID Tissue  Factor Condition Replicate
# 1 K4me1_WT1  mESCs H3K4me1        WT         1

# bamReads
# 1 WT1K4M_03061AAF_GTCATTGCTC-GATCGTAACT_R1_001_merged_trimmed.sorted.unique.markdup.filtered_blacklisted.bam

# bamControl
# 1 I_WT_03057AAF_TGGCTTACAT-GAGGCCGATG_R1_001_merged_trimmed.sorted.unique.markdup.filtered_blacklisted.bam

# Peaks
# 1 WT1K4M_03061AAF_GTCATTGCTC-GATCGTAACT.epic2_results.txt

########################
# 3. Differential Binding Analysis: Complete analysis
########################

# Example of a complete direct analysis. However, it will be performed step by step
# dbObj <- dba(sampleSheet=samples) %>%
#   dba.blacklist() %>%
#   dba.count()     %>%
#   dba.normalize() %>%
#   dba.contrast()  %>%
#   dba.analyze()

# 1. Reading in Peaksets: 
  # For peaks called with epic2 (broad i.e. H3K4me1)
dbObj <- dba(sampleSheet=samples,
             peakCaller = "raw",  # text file with peak score in 4th column
             peakFormat = "raw",  # "sicer" with scoreCol=7 (ChIPcount)
             scoreCol = 4,  
             bLowerScoreBetter= TRUE) # Since score is a p-value

  # For peaks called with MACS2 (narrow i.e. H3K4me3)
# dbObj <- dba(sampleSheet=samples,
#              peakCaller = "narrow")  # For K4 peakset

# 2. Counting on regions of interest
# Summit values were selected to have interval widths between the minimum and 
# first quartile peak width values for each histone mark

dbObj <- dba.count(dbObj, 
                   bRemoveDuplicates=FALSE, # By default
                   bUseSummarizeOverlaps=TRUE,
                   summits = 500)  # To test intervals of 1000bp (H3K27ac, H3K4me1 and H3K27me3)
                  #summits = 150) # To test intervals of 300bp (H3K4me3)
#dbObj

# 6 Samples, 53290 sites in matrix:
#   ID Tissue  Factor Condition Replicate    Reads FRiP
# 1 K4me1_WT1  mESCs H3K4me1        WT         1 79866966 0.08
# 2 K4me1_WT2  mESCs H3K4me1        WT         2 81871116 0.08
# 3 K4me1_WT3  mESCs H3K4me1        WT         3 80063436 0.08
# 4 K4me1_KO1  mESCs H3K4me1        KO         1 77775393 0.08
# 5 K4me1_KO2  mESCs H3K4me1        KO         2 75045111 0.08
# 6 K4me1_KO3  mESCs H3K4me1        KO         3 73892339 0.08
 
# 3. Plot PCA and correlation heatmap just to check how samples cluster

pdf(file = "results/DiffBind/figures/PCA_H3K4me1.pdf", height=6, width=6)
dba.plotPCA(dbObj,  
            attributes=DBA_CONDITION, 
            label=DBA_ID)
dev.off()

  # Correlation Heatmap
pdf(file.path("results/DiffBind/figures/Corr_Heatmap_H3K4me1.pdf"), height=6, width=6)
plot(dbObj)
dev.off()

# 4. Define the constrast of interest KO vs WT
dbObj <- dba.contrast(dbObj,
                      contrast=c("Condition","KO","WT"),
                      minMembers = 3)

# 5. Conduct the diff analysis with DESeq2 and EdgeR
# NOTE: DESeq2 results will be kept for information, EdgeR ARE THE ONES TO BE CONSIDERED

# (RLE normalization for DESeq2 and TMM normalization for EdgeR by default if not performed earlier)
dbObj <- dba.analyze(dbObj, 
                     bBlacklist = FALSE, bGreylist = FALSE,  # We've applied blacklisting before to BAM
                     method=DBA_ALL_METHODS)  # To use DESeq2 and EdgeR

# 6. Show the results of the testing (FDR 5% by default)
dba.show(dbObj, bContrasts=T)	

# H3K4me1
# Factor Group Samples Group2 Samples2 DB.edgeR DB.DESeq2
# 1 Condition    KO       3     WT        3     7989      7647

# 7. (Optional) PCA only with the diff peaks
dba.plotPCA(dbObj, contrast=1, 
            method=DBA_EDGER, 
            attributes=DBA_CONDITION, 
            label=DBA_ID)

# 8. Diff Peaks obtained from DESeq2 or EdgeR
pdf(file = "results/DiffBind/figures/Venn_DESeq2_EdgeR_H3K4me1.pdf", height=6, width=6)
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
dev.off()

# 9. Plot log concentrations KO against WT
pdf(file = "results/DiffBind/figures/MAplot_EdgeR_H3K4me1.pdf", height=6, width=6)
dba.plotMA(dbObj, method=DBA_EDGER)
dev.off()

pdf(file= "results/DiffBind/figures/Conc_KOvsWT_EdgeR_H3K4me1.pdf", height=6, width=6)
dba.plotMA(dbObj, bXY=TRUE,method=DBA_EDGER)   # Red points are the significant
dev.off()

# FINAL: Store the results in an RDS file
saveRDS(dbObj, file = "results/DiffBind/RDS/DiffBind_H3K4me1_results.RDS")

########################
# 4. Annotation of diff binding analysis results with ChIPseeker
# PREAMBLE: prepare data into correct format
########################

# a) Extracting Final results in a dataframe
# =====
db <- dba.report(dbObj, 
                 bDB=TRUE, 
                 bGain=TRUE, bLoss=TRUE,
                 method = DBA_EDGER)
db

#H3K4me1
# 3 Samples, 7989 sites in matrix:
#   Contrast Direction DB Method Intervals
# 1 KO vs. WT       All DB  edgeR      7989
# 2 KO vs. WT      Gain DB  edgeR      3362
# 3 KO vs. WT      Loss DB  edgeR      4627

res <- dba.report(dbObj, method=DBA_EDGER, 
                  contrast = 1, th=1) # Store all intervals
out <- as.data.frame(res)
#head(out)

# b) Create bed files for the gained and lost peaks 
# =====

KO_enrich <- out %>%
   filter(FDR < 0.05 & Fold > 0) 

write.table(x = KO_enrich[,c("seqnames","start","end")], file= "results/DiffBind/BED/H3K4me1_KO_gained_intervals_adjpval.bed", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")

KO_loss <- out %>%
  filter(FDR < 0.05 & Fold < 0) 

write.table(x = KO_loss[,c("seqnames","start","end")], file= "results/DiffBind/BED/H3K4me1_KO_lost_intervals_adjpval.bed", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep="\t")


# c) Create your own TxDb object i.e. from Ensembl if you are not using
#     a built one (i.e. from UCSC). Stick to genome version used during alignment
#     NOTE: This is common for any histone mark included in the paper
# =====

metadata <- data.frame(name="Resource URL",
                       value=paste0("ftp://ftp.ensemblgenomes.org/pub/","release-102/gtf/mus_musculus/"))

txdb <- makeTxDbFromGFF(file = "/Volumes/projectscomput/cancer/db_files/Genomes/Ensembl/mouse/mm10/release-102/Mus_musculus.GRCm38.102.gtf",
                        format="gtf",
                        dataSource="Ensembl_FTP_repository",
                        organism="Mus Musculus",
                        taxonomyId=10090,
                        circ_seqs=NULL,
                        chrominfo=NULL,
                        miRBaseBuild=NA,
                        metadata=metadata,
                        dbxrefTag="gene_id")

# Just to check is all correct
# columns(txdb)
# transcripts(txdb)

# d) Transform the diff bind peaks into a GR object
# =====
Sites.gr <- makeGRangesFromDataFrame(df = out[,-which(colnames(out) %in% c("width","strand"))],
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqnames.field="seqnames",
                                     start.field="start",
                                     end.field="end",
                                     starts.in.df.are.0based=FALSE)  # Ensembl uses a 1-based coordinate


########################
# 5. Peak Annotation with 'annotatePeak' function
########################

peakAnno <- annotatePeak(peak = Sites.gr, 
                         tssRegion = c(-5000,100), 
                         TxDb = txdb,
                         genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                                       "Downstream", "Intergenic"),  # Default genomic Annotation
                         annoDb = "org.Mm.eg.db")

# Let's split the info into different sheets
final <- as.data.frame(peakAnno)

tostore <- list("All_Peaks_H3K4me1" = final,
                "Gained_KO_H3K4me1" = final %>% filter(FDR<0.05,Fold >0),
                "Lost_KO_H3K4me1" = final %>% filter(FDR<0.05,Fold <0))

# Store annotations
hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12, fontName = "Arial Narrow", fgFill = "#4F80BD")
write.xlsx(x= tostore, file = "results/DiffBind/DiffBind_All_Peaks_H3K4me1_ANNOTATED.xlsx", headerStyle=hs)

########################
# 6. Volcano plot with results diffBind
########################

p <- EnhancedVolcano(toptable = final,
                     lab= final$SYMBOL,
                     x ='Fold',
                     y= 'FDR',
                     xlab= 'Log2 Fold change',
                     ylab = "-Log10(adj pval)",
                     title = "Differential Binding H3K4me1",
                     subtitle="KO vs WT IkBa",
                     #xlim=c(-5,5),
                     pCutoff = 0.05,
                     FCcutoff = 0,
                     #caption = "FC cutoff: 0; adj p-val cutoff: 0.05",
                     pointSize = 1.2,
                     labSize = 2.5,
                     max.overlaps = 20,
                     drawConnectors = TRUE,
                     legendLabels=c("Not sign.", "log2FC", "adj pval","adj pval & log2FC"),
                     legendPosition = "right",
                     legendIconSize = 5.0,
                     gridlines.major = FALSE,
                     gridlines.minor = FALSE)

grob <- grobTree(textGrob(paste0("GAINED (adj pval<0.05)\n N=", length(Sites.gr.gained)), x=0.7,  y=0.95, hjust=0,
                          gp=gpar(col="gray43", 
                                  fontsize=7, 
                                  fontface="bold")))

grob2 <- grobTree(textGrob(paste0("LOST (adj pval<0.05)\n N=", length(Sites.gr.lost)), x=0.1,  y=0.95, hjust=0,
                           gp=gpar(col="gray43", 
                                   fontsize=7, 
                                   fontface="bold")))

p <- p + annotation_custom(grob) + annotation_custom(grob2) 


pdf(file = "results/DiffBind/figures/Volcano_DiffBind_H3K4me1.pdf", width=7, height =6)
print(p)
dev.off()

