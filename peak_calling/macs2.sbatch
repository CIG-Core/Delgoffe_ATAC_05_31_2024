#!/bin/bash
#
#SBATCH --job-name=macs2
#SBATCH -c 8
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -o ./log/macs2-%A_%a.out.o%j 
#SBATCH --mail-user=sis186@pitt.edu
#SBATCH --mail-type=END,FAIL    # email me when the job fails or ends
#SBATCH --array=0-14 # job array index
#SBATCH --mem=64g #Add memory

module load python/anaconda2.7-5.2.0 gcc/8.2.0 macs/2.1.1.20160309

SCRIPTPATH=/ix/cigcore/proj/delgoffe/code/Sissi/Delgoffe_ATAC_05_31_2024/bowtie2
names=($(cat $SCRIPTPATH/mouse_bowtie2_listFastq.txt))

echo ${names[${SLURM_ARRAY_TASK_ID}]}

DATAPATH=/ix/cigcore/proj/delgoffe/result/Delgoffe_ATAC_05_31_2024/ATACshift
OUTPATH=/ix/cigcore/proj/delgoffe/result//Delgoffe_ATAC_05_31_2024/peak_calling_MACS2

mkdir -p $OUTPATH

macs2 callpeak \
--treatment $DATAPATH/${names[${SLURM_ARRAY_TASK_ID}]}_shifted.bam \
--name ${names[${SLURM_ARRAY_TASK_ID}]} \
--format BAMPE \
-g mm \
--keep-dup all \
--outdir $OUTPATH/${names[${SLURM_ARRAY_TASK_ID}]}