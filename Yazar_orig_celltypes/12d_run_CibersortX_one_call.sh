#!/bin/bash
#SBATCH --job-name Testing_CibersortX
#SBATCH -N 1
#SBATCH --mem=196G
#SBATCH -t 3-10:30:00
#SBATCH --mail-type=END,ALL
#SBATCH --mail-user=btuftin@unc.edu # send-to address
#SBATCH --output=logs/RNAseq_CibersortX_%j.log
#SBATCH --error=logs/RNAseq_CibersortX_%j.err

ref=$1
bulk=$2
perm=$3

mkdir /proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/Yazar_orig_celltypes/CX_output_2/${ref}_${bulk}_${perm}

singularity exec -B /proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/Yazar_orig_celltypes/datasets_for_comparison/:/src/data \
  -B /proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/Yazar_orig_celltypes/CX_output_2/${ref}_${bulk}_${perm}/:/src/outdir \
  --contain /work/users/b/t/btuftin/CibersortX/fractions_latest.sif /src/CIBERSORTxFractions \
  --verbose TRUE \
  --username btuftin@unc.edu --token d2cce77ff2fa380376c31c66a3543401 \
  --single_cell TRUE --refsample ${ref}_for_cibersortX.txt \
  --mixture ${bulk}_samples_for_cibersortX.txt \
  --perm $perm --rmbatchSmode TRUE

mv /proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/Yazar_orig_celltypes/CX_output_2/${ref}_${bulk}_${perm}/CIBERSORTx_Adjusted.txt \
  /proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/Yazar_orig_celltypes/CX_output/CIBERSORTx_${ref}_${bulk}_${perm}batchcorr_permutations_one_step_run.txt

