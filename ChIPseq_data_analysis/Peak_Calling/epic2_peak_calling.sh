#!/bin/bash

#SBATCH --job-name=epic2
#SBATCH --partition=long
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --nodes=1  
#SBATCH --output=logs/epic2.out
#SBATCH --error=logs/epic2.err
#SBATCH --array=1-12%6

# REMARK!!!! Adapt the number of array tasks to the number of samples i.e. if you have
# 12 samples to analyze you need to specify 12 as indicated. %X means refers to the number 
# of tasks will be sent to cluster execution simultaneously. 
# Each task meaning one sample to be analyzed

#=========================
# User defined parameters: relevant paths and analysis type
#=========================

# SPECIFY Root directory in the cluster (usually /projects/cancer)
ROOTDIR="/projects/cancer"

# SPECIFY your project working directory. 
WKD=$ROOTDIR'/IkBa_mESCs_ABigas/ChIPseq_H3K27ac_H3K4me1'

# SPECIFY the file name where the sample;input;est_fragment_size is included. 
# This is required to know which input corresponds to which sample
SAMPLESHEET=$WKD"/epic2_peak_calling/Samples_Input_epic2.txt"

#=========================
# General configuration
#=========================
START=$(date +%s)
# Enable Singularity image to look into the general path (equivalent to -B)
export SINGULARITY_BIND=$ROOTDIR 
# Path to images folder in cluster
IMAGES_PATH=$ROOTDIR"/images"
# Path to databases folder in cluster
DB_PATH=$ROOTDIR"/db_files"

################################################################################
##       Peak calling with epic2 - (ultraperformance implementation of SICER)
################################################################################

# NOTE: epic2 (and SICER) were especially developed for broad signals 
# Additional information of epic2
# https://academic.oup.com/bioinformatics/article/35/21/4392/5421513
# https://github.com/biocore-ntnu/epic2

###########################
## 1. Other relevant paths
###########################

# Folder where input BAM files are available
DATA=$WKD'/Bowtie_align/BAM_Markdup'

# Folder where to store epic2 output 
OUT=$WKD'/epic2_peak_calling'

#################################################
## 2. Singularity image and Tool Parametrization
#################################################

# Specify image/s name to be used (tool-related)
EPIC='epic2_v0.0.52.sif'

# Specify any particular tool parameters
# Recommended to use these ones adapted to read length 
# https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
GSIZE=2308125349 # mm10 50bp read length 

# Number of processors
T=4

# Genome to be considered. Epic2 already has information about several genomes
GENOME="mm10"

# FDR to consider; default is 5%
FDR=0.05

# Define the chromosome sizes file (OPTIONAL, epic2 has its own by default ==> It's the same with the only 
# difference that it does not include the non-canonical chromosomes)
CHROMSIZES=$DB_PATH'/Genomes/Ensembl/mouse/mm10/release-102/mm10_Ensembl_r102_chrom_only_canonical.sizes'


# About the effective genome size. We will use the same as considered for deepTools 
# https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
# Effective genome size mm10 50bp read length: 2308125349 # ADAPT IN FUNCTION OF READ LENGTH
# Genome length: 2,730,871,774 http://nov2020.archive.ensembl.org/Mus_musculus/Info/Annotation

GS_FRACTION=0.8451973 #2308125349/2730871774

################################################################################
## 3. Command file preparation: to execute batch array
################################################################################

while IFS=";" read -r sample input fragment; 
do
  # Sample name
  NAME=${sample%_R1_001_merged_trimmed.sorted.unique.markdup.filtered_blacklisted.bam}

    # epic2 execution
    echo "singularity exec $IMAGES_PATH/$EPIC epic2 --treatment $DATA/$sample \
    --control $DATA/$input \
    --chromsizes $CHROMSIZES \
    --effective-genome-fraction $GS_FRACTION \
    --fragment-size $fragment \
    --false-discovery-rate-cutoff $FDR \
    --output $OUT/$NAME.epic2_results.txt"
  
done < $SAMPLESHEET > $WKD'/scripts/cmds/epic2_peak_calling.cmd'

################################################################################
## 4. Peak calling
################################################################################

echo "-----------------------------------------"
echo "Starting Peack Calling with Epic2"
echo "-----------------------------------------"

DATE=$(date +%m-%d-%Y--%T)
echo "  Peak calling in array mode: $DATE"
echo " "

SEEDFILE=$WKD'/scripts/cmds/epic2_peak_calling.cmd'
SEED=$(sed -n ${SLURM_ARRAY_TASK_ID}p $SEEDFILE)
eval $SEED

DATE=$(date +%m-%d-%Y--%T)
echo "  Peaks called with Epic2: $DATE"


################################################################################
## 5. End
################################################################################

END=$(date +%s)
DIFF=$(( $END - $START ))
echo 'Peak calling completed' 
echo "Processing Time: $DIFF seconds"
