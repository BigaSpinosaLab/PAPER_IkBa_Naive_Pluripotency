# Chromatin activity of IkBa mediates the exit from naïve pluripotency

This repository includes scripts required for the bulk RNAseq and ChIPseq data analysis included in L. Galán-Palma et al. eLife 2025 manuscript. All scripts include comments so they are self-explanatory.

The repository is organized in the following subfolders:

## RNAseq data analysis folder

Scripts required to reproduce the complete RNAseq data analysis, specifically:

- Data preprocessing: to obtain a raw expression matrix from FASTQ files. Check repository as the lab standard pre-procesing workflow: https://github.com/BigaSpinosaLab/LAB_RNAseq_Data_Analysis
 
- Downstream analysis: script named as *1_DEA_KO_vs_WT.R* to conduct differential expression analysis between IkBa KO and WT conditions at different timepoints.
  
- Related paper figures, specifically for:
  
  - Figure 1F: Heatmap of selected genes across all samples.
  - Figure 2C: First two principal components including all samples.
  - Figure 2E: GSVA Z-scores showing the evolution of Meso/Ecto/Endo formation signatures (GO BP terms).
  - Figure 3B: Random walk plot from GSEA (mESCs RNAseq) against naïve and ground signatures.
  - Figure 5A: Violin plots NfKB related pathway.
  
To conduct data preprocessing, original FASTQ files are required. Please check GEO accession no. GSE239563 (RNAseq data) [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239563]. Preprocessing scripts were executed using Singularity images (v3.8.3).

It is possible to directly conduct the downstream analysis if the corresponding raw expression matrix from GEO is downloaded. In case of the mentioned figures, it is possible to reproduce them if the corresponding normalized (VST) expression matrix from GEO is downloaded.
NOTE: Gene annotation files are required. For this analysis, these were retrieved from Ensembl (release 102, mm10).

## ChIPseq data analysis folder

Scripts required to reproduce the complete ChIPseq data analysis, specifically:

- Data preprocessing: prepare corresponding BAM files ready for peak calling. Check repository as the lab standard pre-procesing workflow: https://github.com/BigaSpinosaLab/LAB_ChIPseq_Data_Analysis
 
- Downstream analysis: scripts in this subfolder include (i) Peak calling with MACS2 (narrow peaks) or epic2 (broad peaks), (ii) peak annotation and consensus peakset per condition,  (iii) Differential Binding analysis based on lists of called peak in different conditions and (iv) Differential Enhancer Activity (based on results of histone marks) and corresponding functional analysis.
  
To conduct data preprocessing, original FASTQ files are required. Please check GEO accession no. GSE239564 (ChIPseq data) [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239564]. 

It is possible to directly conduct differential binding or enhancer analysis if called peak are directly downloaded from GEO dataset. 
