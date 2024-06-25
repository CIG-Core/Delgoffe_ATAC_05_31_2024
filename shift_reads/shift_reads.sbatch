#!/bin/bash
#
#SBATCH --job-name=shift
#SBATCH -c 8
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -o ./log/shift-%A_%a.out.o%j 
#SBATCH --mail-user=sis186@pitt.edu
#SBATCH --mail-type=END,FAIL    # email me when the job fails or ends
#SBATCH --array=0-14 # job array index
#SBATCH --mem=64g #Add memory

module purge
module load bedtools/2.31.0

SCRIPTPATH=/ix/cigcore/proj/delgoffe/code/Sissi/Delgoffe_ATAC_05_31_2024/bowtie2
names=($(cat $SCRIPTPATH/mouse_bowtie2_listFastq.txt))

echo ${names[${SLURM_ARRAY_TASK_ID}]}

DATAPATH=/ix/cigcore/proj/delgoffe/result/Delgoffe_ATAC_05_31_2024/post_align_filter
OUTPATH=/ix/cigcore/proj/delgoffe/result/Delgoffe_ATAC_05_31_2024/shift_reads

mkdir -p $OUTPATH
bamfile=$OUTPATH/${names[${SLURM_ARRAY_TASK_ID}]}_shifted.bam

bedtools bamtobed -i $DATAPATH/${names[${SLURM_ARRAY_TASK_ID}]}_mapq30_noDup_noMT_sort.bam -bedpe | awk -v OFS="\t" '{if($9=="+"){print $1,$2+4,$6+4}else if($9=="-"){print $1,$2-5,$6-5}}' > $bamfile