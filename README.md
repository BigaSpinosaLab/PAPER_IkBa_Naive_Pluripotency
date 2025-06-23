# Chromatin activity of IkBa mediates the exit from naïve pluripotency

This repository includes scripts required for the bulk RNAseq and ChIPseq data analysis included in L. Galán-Palma et al. eLife 2025 manuscript. All scripts include comments so they are self-explanatory.

The repository is organized in the following subfolders:

## RNAseq data analysis folder

Scripts required to reproduce the complete RNAseq data analysis, specifically:

- Data preprocessing: to obtain a raw expression matrix from FASTQ files. Check repository as the lab standard pre-procesing workflow: https://github.com/BigaSpinosaLab/LAB_RNAseq_Data_Analysis
 
- Downstream analysis: to conduct differential expression analysis and functional analysis (GSEA against ground and naïve states signatures from public data).
  
- Other paper figures: to reproduce paper figures, specifically: Figure 2E (evolution Xderm signatures) and Figure 1F (heatmap relevant genes).
  
To conduct data preprocessing, original FASTQ files are required. Please check GEO accession no. GSE239563 (RNAseq data) [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239563]. Preprocessing scripts were executed using Singularity images (v3.8.3).

It is possible to directly conduct the downstream analysis if the corresponding raw expression matrix from GEO is downloaded. NOTE: Gene annotation files are required. For this analysis, these were retrieved from Ensembl (release 102, mm10).

## ChIPseq data analysis folder
