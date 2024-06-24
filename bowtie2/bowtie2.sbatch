#!/bin/bash
#
#SBATCH --job-name=sort_n_bowtie2
#SBATCH -o ./log/bowtie2-%A_%a.out.o%j 
#SBATCH -c 10
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH --mail-user=sis186@pitt.edu
#SBATCH --mail-type=END,FAIL    # email me when the job fails or ends
#SBATCH --array=0-14 # job array index
#SBATCH --mem=60g #Add memory

module load bowtie2/2.4.5
module load gcc/8.2.0 samtools/1.10

DATAPATH=/ix/cigcore/proj/delgoffe/data/Delgoffe_ATAC_05_31_2024
SCRIPTPATH=/ix/cigcore/proj/delgoffe/code/Sissi/Delgoffe_ATAC_05_31_2024/bowtie2
names=($(cat $SCRIPTPATH/mouse_bowtie2_listFastq.txt))

echo ${names[${SLURM_ARRAY_TASK_ID}]}

OUTPATH=/ix/cigcore/proj/delgoffe/result/Delgoffe_ATAC_05_31_2024/bowtie2
mkdir -p $OUTPATH

# better alignment results are frequently achieved with --very-sensitive
# use -X 1000 to allow larger fragment size (default is 500)

bowtie2 --very-sensitive -X 500 --threads 8 \
  -x /ix/cigcore/utils/mouse_indexed/bowtie2/mm10 \
  -1 $DATAPATH/${names[${SLURM_ARRAY_TASK_ID}]}_R1_001.fastq.gz \
  -2 $DATAPATH/${names[${SLURM_ARRAY_TASK_ID}]}_R2_001.fastq.gz \
  | samtools view -u - \
  | samtools sort -n -o $OUTPATH/${names[${SLURM_ARRAY_TASK_ID}]}_bowtie2_sort_n.bam -