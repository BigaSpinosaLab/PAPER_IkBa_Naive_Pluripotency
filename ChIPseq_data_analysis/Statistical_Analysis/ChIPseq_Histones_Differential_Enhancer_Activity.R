
# ==================
# SCOPE:
# This script identifies putative enhancers based on specific histone marks 
# (H3K4me1 and H3K27ac), specifically those diff binded regions identified in 
# our own experiments when comparing KO vs WT peaks 

# NOTE: We will not apply any filter to identified peaks in terms of distance 
# to TSS regions. No criteria in this sense is applied.

# Author: M.Maqueda
# Date: 22nd March 2024
# ==================

# Brief description of enhancers identification approach based on diff regions

# 1. Poised enhancers gained in KO: those ones with gained H3K4me1 in KO with no
# H3K27ac peaks originally (from the Consensus KO set)

# 2. Poised enhancers lost in KO: those ones with lost H3K4me1 in KO with no
# H3K27ac peaks originally (from the Consensus KO set)

# 3. Active enhancers wo K4me1 gained in KO: those ones with gained H3K27ac in KO 
# with no H3K4me1 peaks originally (from the Consensus KO set)

# 4. Active enhancers wo K4me1 lost in KO: those ones with lost H3K27ac in KO 
# with no H3K4me1 peaks originally (from the Consensus KO set)

# 5. Active enhancers w K4me1 gained in KO: those ones with gained H3K27ac in KO 
# with H3K4me1 peaks originally (from the Consensus KO set)

# 6. Active enhancers w K4me1 lost in KO: those ones with lost H3K27ac in KO 
# with H3K4me1 peaks originally (from the Consensus KO set)


# REMARK: As consensus set, we will consider diff binded regions with adj pval < 5% 

################################################################################
# Load libraries
################################################################################

require(ChIPpeakAnno)
require(rtracklayer)
require(dplyr)
require(openxlsx)

# Required for annotation
require(ChIPseeker)  
require(GenomicFeatures)  # For creating a TxDb object from gtf file
require(org.Mm.eg.db)


################################################################################
# Import data
################################################################################

# BED path DBA directory should contain all BED files independently for histone
# mark and gained or lost for KO samples. Check 'ChIPseq_Histones_DBA.R' script
# for reproducibility
bedpath.diffbind <- "../ChIPseq_K27ac_K4me1/results/DiffBind/BED"

# Consider all diff binded regions (adjpval < 5%)
bedfiles <- as.list(list.files(path = bedpath.diffbind, pattern = "_adjpval.bed",full.names = TRUE))

names(bedfiles) <- sapply(list.files(path = bedpath.diffbind, 
                                     pattern = "_adjpval.bed",full.names = FALSE), function(name) 
  paste(unlist(strsplit(name, split="_"))[1:3], collapse="_"))

gr.diff <- lapply(bedfiles, function(bed) rtracklayer::import(con = bed, format="bed"))

# And now import the consensus KO peaks: ONLY those with adjpval < 5% 
# Check 'Consensus_Peakset.R' script for obtaining BED file with this consensus

bedpath.cosensus <- "../ChIPseq_K27ac_K4me1/results/Consensus_PeakSet/BED"
bedfiles <- as.list(list.files(path = bedpath.cosensus, 
                               pattern = "KO_Consensus_Peakset.adjpval.bed",full.names = TRUE))

names(bedfiles) <- sapply(list.files(path = bedpath.cosensus, 
                                     pattern = "KO_Consensus_Peakset.adjpval.bed",
                                     full.names = FALSE), function(name) 
                                       paste(unlist(strsplit(name, split="_"))[c(2,4:5)], collapse="_"))

gr.consensus.loose <- lapply(bedfiles, function(bed) rtracklayer::import(con = bed, format="bed"))

# Define max gap allowed. We have to be loose since regions from DiffBind are just
# 1kb width
# Distance allowed (max gap): +/- X bp between gained peaks and consensus regions
max_gap = 1000

################################################################################
# 1. Poised enhancers gained in KO: those ones with gained H3K4me1 in KO with no
# H3K27ac peaks originally (from the Consensus KO set)

poised.enh.gained.KO <- IRanges::subsetByOverlaps(gr.diff$H3K4me1_KO_gained, 
                                              gr.consensus.loose$H3K27ac_KO_Consensus,
                                                        ignore.strand = TRUE, 
                                                        maxgap = max_gap, 
                                                        invert=TRUE)  
# 810 regions

################################################################################
# 2. Poised enhancers lost in KO: those ones with lost H3K4me1 in KO with no
# H3K27ac peaks originally (from the Consensus KO set)

poised.enh.lost.KO <- IRanges::subsetByOverlaps(gr.diff$H3K4me1_KO_lost, 
                                                  gr.consensus.loose$H3K27ac_KO_Consensus,
                                                  ignore.strand = TRUE, 
                                                  maxgap = max_gap, 
                                                  invert=TRUE)  
# 2,764 regions


################################################################################
# 3. Active enhancers wo K4me1 gained in KO: those ones with gained H3K27ac in KO 
# with no H3K4me1 peaks originally (from the Consensus KO set)

active.enh.gained.KO.woK4 <- IRanges::subsetByOverlaps(gr.diff$H3K27ac_KO_gained, 
                                                gr.consensus.loose$H3K4me1_KO_Consensus,
                                                ignore.strand = TRUE, 
                                                maxgap = max_gap, 
                                                invert=TRUE)  
# 78 regions


################################################################################
# 4. Active enhancers wo K4me1 lost in KO: those ones with lost H3K27ac in KO 
# with no H3K4me1 peaks originally (from the Consensus KO set)

active.enh.lost.KO.woK4 <- IRanges::subsetByOverlaps(gr.diff$H3K27ac_KO_lost, 
                                                       gr.consensus.loose$H3K4me1_KO_Consensus,
                                                       ignore.strand = TRUE, 
                                                       maxgap = max_gap, 
                                                       invert=TRUE)  
# 349 regions


################################################################################
# 5. Active enhancers w K4me1 gained in KO: those ones with gained H3K27ac in KO 
# with H3K4me1 peaks originally (from the Consensus KO set)

active.enh.gained.KO.wK4 <- IRanges::subsetByOverlaps(gr.diff$H3K27ac_KO_gained, 
                                                     gr.consensus.loose$H3K4me1_KO_Consensus,
                                                     ignore.strand = TRUE, 
                                                     maxgap = max_gap, 
                                                     invert=FALSE)  
# 2,468 regions

################################################################################
# 6. Active enhancers w K4me1 lost in KO: those ones with lost H3K27ac in KO 
# with H3K4me1 peaks originally (from the Consensus KO set)

active.enh.lost.KO.wK4 <- IRanges::subsetByOverlaps(gr.diff$H3K27ac_KO_lost, 
                                                      gr.consensus.loose$H3K4me1_KO_Consensus,
                                                      ignore.strand = TRUE, 
                                                      maxgap = max_gap, 
                                                      invert=FALSE)  
# 1,549 regions



################################################################################
# Combine all results and annotate
################################################################################

enhancers.database.gr = list("Poised.Gained.KO" = poised.enh.gained.KO,
                             "Poised.Lost.KO" = poised.enh.lost.KO,
                             "Active.Gained.KO.wo.K4" = active.enh.gained.KO.woK4,
                             "Active.Lost.KO.wo.K4" = active.enh.lost.KO.woK4,
                             "Active.Gained.KO.wK4" = active.enh.gained.KO.wK4,
                             "Active.Lost.KO.wK4" = active.enh.lost.KO.wK4)

# Create the txdb according to the same used during alignment
metadata <- data.frame(name="Resource URL - Filtered",
                       value=paste0("ftp://ftp.ensemblgenomes.org/pub/","release-102/gtf/mus_musculus/"))

txdb <- makeTxDbFromGFF(file = "/Volumes/projectscomput/cancer/db_files/Genomes/Ensembl/mouse/mm10/release-102/Mus_musculus.GRCm38.102_only_ProteinCoding_genes.gtf",
                        format="gtf",
                        dataSource="Ensembl_FTP_repository",
                        organism="Mus Musculus",
                        taxonomyId=10090,
                        circ_seqs=NULL,
                        chrominfo=NULL,
                        miRBaseBuild=NA,
                        metadata=metadata,
                        dbxrefTag="gene_id")

# Add annotations to all individual peak sets
peakAnno <- lapply(enhancers.database.gr, function(set){
  p <-annotatePeak(peak = set, 
                   tssRegion = c(-5000,100), 
                   TxDb = txdb,
                   genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron",
                                                 "Downstream", "Intergenic"),  # Default genomic Annotation
                   annoDb = "org.Mm.eg.db")
  return(p)
})

# Store annotations in a file
peakAnno.df <- lapply(peakAnno, function(k) as.data.frame(k))

hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12, fontName = "Arial Narrow", fgFill = "#4F80BD")
write.xlsx(x= peakAnno.df, 
           file = "results/Differential_Enhancers_in_IkBa_KO_samples.xlsx", headerStyle=hs)


