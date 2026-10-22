#!/bin/bash
#SBATCH --job-name Testing_CibersortX
#SBATCH -N 1
#SBATCH --mem=32G
#SBATCH -t 2:30:00
#SBATCH --mail-type=END,ALL
#SBATCH --mail-user=btuftin@unc.edu # send-to address
#SBATCH --output=logs/RNAseq_CibersortX_%j.log
#SBATCH --error=logs/RNAseq_CibersortX_%j.err

ref="LM22"
bulk="400_JHS"
perm=5

mkdir /proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/Yazar_orig_celltypes/CX_output_2/${ref}_${bulk}_${perm}

singularity exec -B /proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/Yazar_orig_celltypes/Cibersort/:/src/data \
  -B /proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/Yazar_orig_celltypes/CX_output_2/${ref}_${bulk}_${perm}/:/src/outdir \
  --contain /work/users/b/t/btuftin/CibersortX/fractions_latest.sif /src/CIBERSORTxFractions \
  --verbose TRUE \
  --username btuftin@unc.edu --token d2cce77ff2fa380376c31c66a3543401 \
  --sigmatrix LM22.txt \
  --mixture ${bulk}_samples_for_LM22cibersortX.txt \
  --perm $perm  --rmbatchBmode TRUE

mv /proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/Yazar_orig_celltypes/CX_output_2/${ref}_${bulk}_${perm}/CIBERSORTx_Adjusted.txt \
  /proj/raffieldlab/users/btuftin/7702_RNAseq_JHS/Yazar_orig_celltypes/CX_output/CIBERSORTx_${ref}_${bulk}_${perm}batchcorr_permutations_one_step_run.txt

