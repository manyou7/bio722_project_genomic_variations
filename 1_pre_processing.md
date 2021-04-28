---
Title: "Pre-processing"
Author: "Man You"
Date: "10/04/2021"
Output: html_document
---
Back to [README](README.md)

Illumina raw reads of two parents and their progeny were in fastq files. Before using these reads, we need to go through pre-processing steps, including quality control (QC) by ```FastQC``` and read trimming by ```Trimommatic```. Both tools are commonly used for pre-processing steps.

## Step 1: QC of raw reads

The first thing after receiving the raw reads is to check the quality and make sure the quality is good for further analysis.

```{bash}
cd raw_reads/
mkdir fastqc_out_dir/
fastqc *.fastq -o fastqc_out_dir/
```

## Step 2: Read trimming

The choice of adapters was based on the provided sequencing information by the sequencing lab. In my case, they used the TruSeq DNA PCR-FREE and TRuSeq DNA CD Indexes for library prep.

```{bash}
cd ..
mkdir reads ; cd reads
nano trim.sh
bash trim.sh
```

```{bash}
#trim.sh script

#!/bin/bash

path_to_data=/scratch/youman7/raw_reads/

samples="JF109
CDC15
YM165"

for individual in $samples

do
    # echo which individual we are on 
    
    echo ${individual}
    
    # use trimmomatic to trim the raw reads
    
    java -Xmx2G -jar trimmomatic-0.39.jar PE -phred33 -trimlog logfile 
${path_to_data}/${individual}_R1.fastq ${path_to_data}/${individual}_R2.fastq out.${individual}_1.fastq out.trim.${individual}_1.fastq out.${individual}_2.fastq out.trim.${individual}_2.fastq ILLUMINACLIP:/home/youman7/biosoft/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:50

done
```

## Step 3: QC of trimmed reads

After trimming, we expect the reads are in good quality, but we still need to double check. 

```{bash}
mkdir fastqc_trimmed_dir/
fastqc *.fastq -o fastqc_trimmed_dir/
```

Good. They are all good!

Now, all reads are ready to use for aligning. Please click [HERE](2_mapping.md)
