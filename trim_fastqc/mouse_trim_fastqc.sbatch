#!/bin/bash
#
#SBATCH -J qc    # job name
#SBATCH -o /ix/cigcore/proj/delgoffe/code/Sissi/Delgoffe_ATAC_05_31_2024/trim_fastqc/log/fastqc-%A_%a.out.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 1                # total number of mpi tasks requested
#SBATCH -N 1                # single node use
##SBATCH -p normal          # queue (partition) -- normal, development, etc.
#SBATCH -t 20:00:00          # run time (hh:mm:ss)
#SBATCH --mail-user=sis186@pitt.edu
#SBATCH --mail-type=END,FAIL    # email me when the job fails or ends
#SBATCH --array=0-29 # job array index
#SBATCH --mem=60g #Add memory

echo "Started at:"
echo $(date)

module purge
module load fastqc/0.11.9

DATAPATH=/ix/cigcore/proj/delgoffe/data/Delgoffe_ATAC_05_31_2024/cleaned
SCRIPTPATH=/ix/cigcore/proj/delgoffe/code/Sissi/Delgoffe_ATAC_05_31_2024/trim_fastqc
FASTQNAMES=($(cat $SCRIPTPATH/mouse_trim_fastqc_listFastq.txt))

echo "Doing fastqc analysis on:"
echo ${FASTQNAMES[${SLURM_ARRAY_TASK_ID}]}

QCPATH=/ix/cigcore/proj/delgoffe/data/Delgoffe_ATAC_05_31_2024/qc
mkdir -p ${QCPATH}/trim_fastqc
## The job command:
fastqc -o ${QCPATH}/trim_fastqc \
        $DATAPATH/${FASTQNAMES[${SLURM_ARRAY_TASK_ID}]}

echo "Finished at:"
echo $(date)
