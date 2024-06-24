#!/bin/bash
#
#SBATCH -J trimmomatic    # job name
#SBATCH -o /ix/cigcore/proj/delgoffe/code/Sissi/Delgoffe_ATAC_05_31_2024/trimmomatic/log/trimmomatic-%A_%a.out.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 1                # total number of mpi tasks requested
#SBATCH -N 1                # single node use
##SBATCH -p normal          # queue (partition) -- normal, development, etc.
#SBATCH -t 20:00:00          # run time (hh:mm:ss)
#SBATCH --mail-user=sis186@pitt.edu
#SBATCH --mail-type=END,FAIL     # email me when the job fails or ends
#SBATCH --array=0-29 # job array index
########################################

## trimmomatic set-up
module purge 
module load trimmomatic/0.38

set -x
#########################
FASTQPATH=/ix/cigcore/proj/delgoffe/data/Delgoffe_ATAC_05_31_2024
CLEANPATH=/ix/cigcore/proj/delgoffe/data/Delgoffe_ATAC_05_31_2024/cleaned
INTERIMPATH=/ix/cigcore/proj/delgoffe/data/Delgoffe_ATAC_05_31_2024/intermediate
SCRIPTPATH=/ix/cigcore/proj/delgoffe/code/Sissi/Delgoffe_ATAC_05_31_2024/trimmomatic
FASTQNAMES=($(cat $SCRIPTPATH/mouse_trim_listFastq.txt))

echo "Doing fastqc analysis on:"
echo ${FASTQNAMES[${SLURM_ARRAY_TASK_ID}]}

##########################
mkdir -p $CLEANPATH
mkdir -p $INTERIMPATH

################
echo "Trimming started at:"
echo $(date)

##Running trimmomatic on the .fastq file

java -jar /ihome/crc/install/trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar PE \
	-threads 63 -phred33 \
	$FASTQPATH/${FASTQNAMES[${SLURM_ARRAY_TASK_ID}]}_R1_001.fastq.gz $FASTQPATH/${FASTQNAMES[${SLURM_ARRAY_TASK_ID}]}_R2_001.fastq.gz \
	$CLEANPATH/${FASTQNAMES[${SLURM_ARRAY_TASK_ID}]}_trimmomatic_R1.fastq \
	$INTERIMPATH/${FASTQNAMES[${SLURM_ARRAY_TASK_ID}]}_trimmomatic_unpaired_R1.fastq \
	$CLEANPATH/${FASTQNAMES[${SLURM_ARRAY_TASK_ID}]}_trimmomatic_R2.fastq \
	$INTERIMPATH/${FASTQNAMES[${SLURM_ARRAY_TASK_ID}]}_trimmomatic_unpaired_R2.fastq \
	ILLUMINACLIP:/ihome/crc/install/trimmomatic/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:50

echo "Finished at:"
echo $(date)
