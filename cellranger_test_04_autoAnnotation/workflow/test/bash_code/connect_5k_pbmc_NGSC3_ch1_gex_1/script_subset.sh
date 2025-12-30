#!/bin/bash
#SBATCH --job-name=seqtk_3_subsets
#SBATCH --account pedrini.edoardo
#SBATCH --mem=128GB
#SBATCH --time=INFINITE 
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="subsample_log.err"
#SBATCH --output="subsample_log.out"

echo "Job started: $(date)" > process.log

. /home/pedrini.edoardo/miniconda3/bin/activate;
conda activate bioinfo;

# --- Define Input Files (The concatenated files you already have) ---
R1_BIG="connect_5k_pbmc_NGSC3_ch1_gex_1_S1_L003_R1_001.fastq.gz"
R2_BIG="connect_5k_pbmc_NGSC3_ch1_gex_1_S1_L003_R2_001.fastq.gz"

# This file will track all reads we have already used
BLACKLIST="used_read_ids.txt"
touch $BLACKLIST

# --- SUBSET 1 (5% - L004) ---
echo "Creating Subset 1 (L004)..." >> process.log
seqtk sample -s42 $R1_BIG 0.05 | gzip > connect_5k_pbmc_NGSC3_ch1_gex_1_S1_L004_R1_001.fastq.gz
seqtk sample -s42 $R2_BIG 0.05 | gzip > connect_5k_pbmc_NGSC3_ch1_gex_1_S1_L004_R2_001.fastq.gz

# Update blacklist with IDs from the new R1 file
zcat connect_5k_pbmc_NGSC3_ch1_gex_1_S1_L004_R1_001.fastq.gz | awk 'NR%4==1 {print $1}' | sed 's/^@//' >> $BLACKLIST

# --- SUBSET 2 (5% - L005) ---
echo "Creating Subset 2 (L005)..." >> process.log
# Filter out used reads from the BIG files
zgrep -v -F -f $BLACKLIST $R1_BIG | gzip > temp_R1.fastq.gz
zgrep -v -F -f $BLACKLIST $R2_BIG | gzip > temp_R2.fastq.gz

# Sample from the remainders
seqtk sample -s42 temp_R1.fastq.gz 0.05 | gzip > connect_5k_pbmc_NGSC3_ch1_gex_1_S1_L005_R1_001.fastq.gz
seqtk sample -s42 temp_R2.fastq.gz 0.05 | gzip > connect_5k_pbmc_NGSC3_ch1_gex_1_S1_L005_R2_001.fastq.gz

# Add Subset 2 IDs to the blacklist
zcat connect_5k_pbmc_NGSC3_ch1_gex_1_S1_L005_R1_001.fastq.gz | awk 'NR%4==1 {print $1}' | sed 's/^@//' >> $BLACKLIST

# --- SUBSET 3 (5% - L006) ---
echo "Creating Subset 3 (L006)..." >> process.log
# Note: $BLACKLIST now contains IDs from L004 AND L005
zgrep -v -F -f $BLACKLIST $R1_BIG | gzip > temp_R1.fastq.gz
zgrep -v -F -f $BLACKLIST $R2_BIG | gzip > temp_R2.fastq.gz

seqtk sample -s42 temp_R1.fastq.gz 0.05 | gzip > connect_5k_pbmc_NGSC3_ch1_gex_1_S1_L006_R1_001.fastq.gz
seqtk sample -s42 temp_R2.fastq.gz 0.05 | gzip > connect_5k_pbmc_NGSC3_ch1_gex_1_S1_L006_R2_001.fastq.gz

# --- Cleanup ---
rm temp_R1.fastq.gz temp_R2.fastq.gz
echo "Job finished: $(date)" >> process.log