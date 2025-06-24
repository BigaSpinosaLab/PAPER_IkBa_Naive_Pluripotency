################################################################################
##        Title : mESCs IkBa - ChIPseq Histone marks
##  Description : This script is for creating a consensus PeakSet among the three 
##                replicates per condition (WT and KO in this case) per histone
##   Researcher : Luis Galán
##       Author : María Maqueda
##         Date : 21st March 2024
################################################################################

# Objective: Create consensus peakset (peaks in at least two out of the three
# replicates) and annotate
# Annotation will also be conducted for the peaks individually called (per sample)

# REMARK: Code refers to H3K4me1 but it can be applied to the rest of marks

################################################################################
## 0. Load packages
################################################################################

require(ChIPpeakAnno)
library(rtracklayer)

# Required for annotation
require(ChIPseeker)  
require(GenomicFeatures)  # For creating a TxDb object from gtf file
require(org.Mm.eg.db)

################################################################################
## 1. Import files with called peaks
################################################################################

results.path <- "/Volumes/cancer/IkBa_mESCs_ABigas/ChIPseq_H3K27ac_H3K4me1/epic2_peak_calling"
files <- list.files(path = results.path,full.names = TRUE)

# K4me1 samples
# ---------------
k4.files <- files[grep("K4M_",files)]

k4.files <- lapply(k4.files, function(f){
  read.table(file = f, header=TRUE,comment.char = "")})

filenames <-  list.files(path = results.path,full.names = FALSE)
filenames <- filenames[grep("K4M_",filenames)]
filenames <- unlist(sapply(filenames, function(fn) unlist(strsplit(fn, split="_"))[1]))

names(k4.files) <- filenames


# DEFINE THE SET OF PEAKS TO BE ANALYZED
peaks.set <- k4.files

################################################################################
## 2. Number of peaks per sample
################################################################################

# Number of peaks only based on adjpval < 0.05
lapply(peaks.set, nrow)

# $KO1K4M
# [1] 70180
# 
# $KO2K4M
# [1] 73373
# 
# $KO3K4M
# [1] 76172
# 
# $WT1K4M
# [1] 72952
# 
# $WT2K4M
# [1] 72126
# 
# $WT3K4M
# [1] 72485

# Or alternatively Number of peaks FDR adj pval 5% and log2FC >1 (or Score >100). 
lapply(peaks.set, function(list) nrow(list %>% filter(log2FoldChange >1)))

# $KO1K4M
# [1] 35096
# 
# $KO2K4M
# [1] 38394
# 
# $KO3K4M
# [1] 41292
# 
# $WT1K4M
# [1] 36236
# 
# $WT2K4M
# [1] 34892
# 
# $WT3K4M
# [1] 35038

# If more restrictive filter wants to be applied
#peaks.set <- lapply(peaks.set, function(list) list %>% filter(log2FoldChange >1))

################################################################################
## 3. Create GRanges objects from those peak files
################################################################################

# FIRST, peaks list must be converted GRanges objects
peaks.set.gr <- lapply(peaks.set, function(res) 
  makeGRangesFromDataFrame(df = res,
                           keep.extra.columns=TRUE,
                           seqnames.field="X.Chromosome",
                           start.field="Start",
                           end.field="End",
                           strand.field = "Strand",
                           starts.in.df.are.0based=FALSE)) 
peaks.set.gr <- GRangesList(peaks.set.gr)

# Let's add some names to the peaks 
samples <- names(peaks.set.gr)
peaks.set.gr <-GRangesList(lapply(samples, function(regions) 
  {
  names(peaks.set.gr[[regions]]) <- paste(regions,seq(1,length(peaks.set.gr[[regions]])), sep="_")
  return(peaks.set.gr[[regions]])
}))
names(peaks.set.gr) <- samples

################################################################################
## 4. Find overlapping peaks (by at least 1bp). These will be merged
################################################################################

# Overlapping peaks
overlapping.ko <- findOverlapsOfPeaks(peaks.set.gr$KO1K4M,peaks.set.gr$KO2K4M, peaks.set.gr$KO3K4M,
                                    ignore.strand = TRUE, 
                                    connectedPeaks = "keepAll",
                                    maxgap = -1)  # Maxgap default value: 1base in common

overlapping.wt <- findOverlapsOfPeaks(peaks.set.gr$WT1K4M, peaks.set.gr$WT2K4M, peaks.set.gr$WT3K4M,
                                      ignore.strand = TRUE, 
                                      connectedPeaks = "keepAll",
                                      maxgap = -1)  # Maxgap default value: 1base in common

write.table(x = as.table(overlapping.ko$venn_cnt), 
            file = "results/Consensus_PeakSet/other/Venn_Counts_K4me1KO_Replicates.txt",
            row.names = FALSE, quote = FALSE)

write.table(x = as.table(overlapping.wt$venn_cnt), 
            file = "results/Consensus_PeakSet/other/Venn_Counts_K4me1WT_Replicates.txt",
            row.names = FALSE, quote = FALSE)

################################################################################
## 5. Create txdb for Annotate the consensus PeakSet
################################################################################

# Create the txdb according to the same used during alignment
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


################################################################################
## 6. Annotate the consensus PeakSet and the individual peaks
################################################################################

# Consensus
Consensus.ko <- annotatePeak(peak = overlapping.ko$mergedPeaks, 
                 tssRegion = c(-5000,100), 
                 TxDb = txdb,
                 genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                               "Downstream", "Intergenic"),  # Default genomic Annotation
                 annoDb = "org.Mm.eg.db")
######
Consensus.wt <- annotatePeak(peak = overlapping.wt$mergedPeaks, 
                             tssRegion = c(-5000,100), 
                             TxDb = txdb,
                             genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                                           "Downstream", "Intergenic"),  # Default genomic Annotation
                             annoDb = "org.Mm.eg.db")

# Add annotations to all individual peak sets
peakAnno <- lapply(peaks.set.gr, function(set){
  p <-annotatePeak(peak = set, 
               tssRegion = c(-5000,100), 
               TxDb = txdb,
               genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                                         "Downstream", "Intergenic"),  # Default genomic Annotation
              annoDb = "org.Mm.eg.db")
  p <- cbind(names(set), as.data.frame(p))
  return(as.data.frame(p))
})

# Store a BED file for corresponding consensus WT or KO
rtracklayer::export.bed(object = overlapping.wt$mergedPeaks, 
                        con= "results/Consensus_PeakSet/BED/ChIP_H3K27ac_IkBa_WT_Consensus_Peakset.adjpval.bed")
rtracklayer::export.bed(object = overlapping.ko$mergedPeaks, 
                        con= "results/Consensus_PeakSet/BED/ChIP_H3K27ac_IkBa_KO_Consensus_Peakset.adjpval.bed")

################################################################################
## 7. Store final results
################################################################################

# Store annotations: Consensus + Individual peaks + consensus of consensus (& exclusive)
tostore <- peakAnno
tostore$Consensus.WT <- as.data.frame(Consensus.wt)
tostore$Consensus.KO <- as.data.frame(Consensus.ko)

hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12, fontName = "Arial Narrow", fgFill = "#4F80BD")
write.xlsx(x= tostore, file = "results/Consensus_PeakSet/ChIP_H3K4me1_IkBa_WT_KO_Consensus_PeakSet.xlsx", headerStyle=hs)

################################################################################
## ANNEX: Peak width distribution -> REQUIRED for summits value for DiffBind
################################################################################

lapply(peaks.set.gr, function(l) summary(width(l)))

# $KO1K4M
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 400    2400    4400    8387    9600  246400 
# 
# $KO2K4M
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 400    2400    4200    7938    9000  246800 
# 
# $KO3K4M
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 400    2200    4000    7520    8600  270200 
# 
# $WT1K4M
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 400    2400    4400    8044    9200  247200 
# 
# $WT2K4M
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 400    2600    4600    8214    9400  257600 
# 
# $WT3K4M
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 400    2600    4600    8175    9400  240000 

