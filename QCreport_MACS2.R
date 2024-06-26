library(ChIPQC)
library(DiffBind)
library(tidyverse)
library(ChIPseeker)

proj_path <- "/ix/cigcore/proj/delgoffe"
data_path <- file.path(proj_path, "data", "Delgoffe_ATAC_05_31_2024")
resPath <- gsub(pattern="data", replacement="result", data_path)
after_filter_directory <- file.path(resPath, "post_align_filter")
MACS2_path <- file.path(resPath, "peak_calling_MACS2")

coldata <- data.frame(
  names = c("0hr_Fasted1_S1", "0hr_Fasted2_S2", "0hr_Fasted3_S3", "0hr_Fed1_S4", "0hr_Fed2_S5", "0hr_Fed3_S6", "24hr_Fasted1_S7", "24hr_Fasted2_S8", "24hr_Fasted3_S9", "24hr_Fed1_S10", "24hr_Fed2_S11", "24hr_Fed3_S12", "24hr_Fed_SREBP1i-1_S13", "24hr_Fed_SREBP1i-2_S14", "24hr_Fed_SREBP1i-3_S15"),
  condition = c("0hr", "0hr", "0hr", "0hr", "0hr", "0hr", "24hr", "24hr", "24hr", "24hr", "24hr", "24hr", "24hr", "24hr", "24hr"),
  factor = c("Fasted", "Fasted", "Fasted", "Fed", "Fed", "Fed", "Fasted", "Fasted", "Fasted", "Fed", "Fed", "Fed", "Fed_SREBP1", "Fed_SREBP1", "Fed_SREBP1")
)

coldata$sample_id <- sub(".*_(S[0-9]+)$", "\\1", coldata$names)

samples <- data.frame(
  SampleID = c("0hr_Fasted1_S1", "0hr_Fasted2_S2", "0hr_Fasted3_S3", "0hr_Fed1_S4", "0hr_Fed2_S5", 
               "0hr_Fed3_S6", "24hr_Fasted1_S7", "24hr_Fasted2_S8", "24hr_Fasted3_S9", "24hr_Fed1_S10", 
               "24hr_Fed2_S11", "24hr_Fed3_S12", "24hr_Fed_SREBP1i.1_S13", "24hr_Fed_SREBP1i.2_S14", 
               "24hr_Fed_SREBP1i.3_S15"),
  Tissue = rep("NA", 15),
  Factor = c("Fasted", "Fasted", "Fasted", "Fed", "Fed", "Fed", "Fasted", "Fasted", "Fasted", 
             "Fed", "Fed", "Fed", "Fed_SREBP1", "Fed_SREBP1", "Fed_SREBP1"),
  Condition = c("0hr", "0hr", "0hr", "0hr", "0hr", "0hr", "24hr", "24hr", "24hr", 
                "24hr", "24hr", "24hr", "24hr", "24hr", "24hr"),
  Replicate = c("1", "2", "3", "1", "2", "3", "1", "2", "3", "1", "2", "3", "1", "2", "3"),
  bamReads = file.path(after_filter_directory, paste0(coldata$names, "_mapq30_noDup_noMT_sort.bam")),
  Peaks = file.path(MACS2_path, coldata$names, paste0(coldata$names, "_peaks.narrowPeak")),
  PeakCaller = rep("bed", 15),
  ControlID = rep("NA", 15),
  bamControl = rep("NA", 15)
)
samples$SampleID <- make.names(samples$SampleID)
metadata_file <- "/ix/cigcore/proj/delgoffe/code/Sissi/Delgoffe_ATAC_05_31_2024/metadata.csv"
write.csv(samples, file = metadata_file, row.names = FALSE)

samples <- read.csv("/ix/cigcore/proj/delgoffe/code/Sissi/Delgoffe_ATAC_05_31_2024/metadata.csv")


aug_samples <- read.csv("/ix/cigcore/proj/delgoffe/data/Delgoffe_Cut_Tag_9_15_22/Delgoffe_Cut_Tag_Metadata.csv")

if (!all(file.exists(samples$bamReads))) {
  stop("Some bamReads files do not exist.")
}
if (!all(file.exists(samples$Peaks))) {
  stop("Some Peaks files do not exist.")
}

BiocParallel::register(BiocParallel::SerialParam())
atacObj <- ChIPQC(samples, annotation="mm10", chromosomes = "chr19", blacklist = "/ix/cigcore/proj/delgoffe/code/Sissi/Delgoffe_ATAC_05_31_2024/mm10.blacklist.bed")
#CPT1A chr19 chr19:3,321,326-3,386,735
ChIPQCreport(atacObj, reportName="ATAC QC report: Delgoffe_ATAC_05_31_2024", reportFolder="QCreport_MACS2")
