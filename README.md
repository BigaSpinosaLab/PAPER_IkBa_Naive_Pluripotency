# IkBa_Naive_Pluripotency

This repository includes scripts required for the bulk RNAseq and ChIPseq data analysis included in L. Galán-Palma et al. eLife 2025 manuscript. All scripts include comments so they are self-explanatory.

The repository is organized in the following subfolders:

## RNAseq data analysis folder

Scripts required to reproduce the complete RNAseq data analysis, specifically:

- Data preprocessing: to obtain a raw expression matrix from FASTQ files. Check repository as the lab standard pre-procesing workflow: https://github.com/BigaSpinosaLab/LAB_RNAseq_Data_Analysis
 
- Downstream analysis: to conduct differential expression analysis and functional analysis (GSEA and Overrepresentation analysis). Scripts 8 and 12. NOTE: Although these scripts are particularized for RNAseq inducible IkBa dataset, they have also been applied for the RNAseq Knock-in dataset. Experimental design has been adapted for the corresponding comparisons.
To conduct data preprocessing, original FASTQ files are required. Please check GEO accession no. GSE206515 (Inducible IKBa) [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206515] or GSE292288 (Knock-in IKBa) [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292288]. All required scripts were executed using Singularity images (v3.8.3) per required tool.

It is possible to directly conduct the downstream analysis if the corresponding raw expression matrix from GEO is downloaded. NOTE: Gene annotation files are required. For this analysis, these were retrieved from Ensembl (release 106, hg38).

## ChIPseq data analysis folder
