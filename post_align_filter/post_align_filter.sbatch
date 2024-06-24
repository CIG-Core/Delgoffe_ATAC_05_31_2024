#!/bin/bash
#
#SBATCH --job-name=filter
#SBATCH -c 8
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -o ./log/filter-%A_%a.out.o%j 
#SBATCH --mail-user=sis186@pitt.edu
#SBATCH --mail-type=END,FAIL    # email me when the job fails or ends
#SBATCH --array=0-14 # job array index
#SBATCH --mem=64g #Add memory

module load gcc/8.2.0 samtools/1.9

SCRIPTPATH=/ix/cigcore/proj/delgoffe/code/Sissi/Delgoffe_ATAC_05_31_2024/bowtie2
names=($(cat $SCRIPTPATH/mouse_bowtie2_listFastq.txt))

echo ${names[${SLURM_ARRAY_TASK_ID}]}

DATAPATH=/ix/cigcore/proj/delgoffe/result/Delgoffe_ATAC_05_31_2024/bowtie2
OUTPATH=/ix/cigcore/proj/delgoffe/result/Delgoffe_ATAC_05_31_2024/post_align_filter

## set up file names
align_bam=$DATAPATH/${names[${SLURM_ARRAY_TASK_ID}]}_bowtie2_sort_n.bam
remove_MT=$OUTPATH/${names[${SLURM_ARRAY_TASK_ID}]}_noMT.bam
remove_Dup=$OUTPATH/${names[${SLURM_ARRAY_TASK_ID}]}_noDup_noMT.bam
MAPQ30=$OUTPATH/${names[${SLURM_ARRAY_TASK_ID}]}_mapq30_noDup_noMT.bam
sort_index=$OUTPATH/${names[${SLURM_ARRAY_TASK_ID}]}_mapq30_noDup_noMT_sort.bam

# Remove the mitochondrial reads after alignment
samtools view -@ 8 -h $align_bam | grep -v chrM | samtools sort -@ 8 -O bam -o $remove_MT

# Remove PCR duplicates 
samtools view -h -b -F 1024 $remove_MT -o $remove_Dup

# Remove multi-mapped reads (i.e. those with MAPQ < 30, using -q in SAMtools)
samtools view -h -b -q 30 $remove_Dup -o $MAPQ30

# Output bam file needs to be sorted/index again
samtools sort -@ 8 -O bam -o $sort_index $MAPQ30
samtools index -@ 8 $sort_index

# Remove intermediate files
rm $remove_MT
rm $remove_Dup
rm $MAPQ30