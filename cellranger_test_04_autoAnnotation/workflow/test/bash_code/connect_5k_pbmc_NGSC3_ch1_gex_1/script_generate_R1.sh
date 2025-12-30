#!/bin/bash
#SBATCH --job-name=seqtk
#SBATCH --account pedrini.edoardo
#SBATCH --mem=128GB  # amout of RAM in MB required (and max ram available).
#SBATCH --time=INFINITE  ## OR #SBATCH --time=10:00 means 10 minutes OR --time=01:00:00 means 1 hour
#SBATCH --ntasks=4  # number of required cores
#SBATCH --nodes=1  # not really useful for not mpi jobs
#SBATCH --mail-type=FAIL ## BEGIN, END, FAIL or ALL
#SBATCH --mail-user=pedrini.edoardo@hsr.it
#SBATCH --error="script_generate_R1.err"
#SBATCH --output="script_generate_R1.out"

echo "my job strart now" > script_generate_R1.log;

date >> script_generate_R1.log;

. /home/pedrini.edoardo/miniconda3/bin/activate;
conda activate bioinfo;

# 00 concatenate the fastqs per read to make it easier the processing
zcat connect_5k_pbmc_NGSC3_ch1_gex_1_S1_L001_R1_001.fastq.gz connect_5k_pbmc_NGSC3_ch1_gex_1_S1_L002_R1_001.fastq.gz | gzip > connect_5k_pbmc_NGSC3_ch1_gex_1_S1_L003_R1_001.fastq.gz

date >> script_generate_R1.log;