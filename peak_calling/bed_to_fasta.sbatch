#!/bin/bash
#
#SBATCH --job-name=fasta
#SBATCH -c 8
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -o ./log/fasta-%A_%a.out.o%j
#SBATCH --mail-user=sis186@pitt.edu
#SBATCH --mail-type=END,FAIL    # email me when the job fails or ends
#SBATCH --array=0-3 # job array index
#SBATCH --mem=64g #Add memory

module purge
module load bedtools/2.31.0

SCRIPTPATH=/ix/cigcore/proj/delgoffe/code/Sissi/Delgoffe_ATAC_1_25_23_Alok/MACS2
names=($(cat $SCRIPTPATH/macs2.txt))

echo ${names[${SLURM_ARRAY_TASK_ID}]}

INDEXPATH=/ix/cigcore/utils/mouse_indexed/bowtie2/Archive
OUTPATH=/ix/cigcore/proj/delgoffe/result/Delgoffe_ATAC_1_25_23_Alok/MACS2/fasta

bedtools getfasta -fi \
$INDEXPATH/mm10.fa \
-bed $OUTPATH/${names[${SLURM_ARRAY_TASK_ID}]}_summits.bed \
-fo $OUTPATH/${names[${SLURM_ARRAY_TASK_ID}]}_summits.fasta