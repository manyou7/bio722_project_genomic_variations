---
Title: "Mapping and Polishing"
Author: "Man You"
Date: "12/04/2021"
Output: html_document
---
Back to [pre-precessing](https://github.com/manyou7/bio722_project_genomic_variations/blob/3361a74fb74382bfb4f92567513722da75cf1b3d/1_pre_processing.md) or [README](https://github.com/manyou7/bio722_project_genomic_variations/blob/01e4810d1d0725e7a8bde9d12f16d0368fbc6666/README.md)

## Read mapping

Now let's map the trimmed reads to the reference genome (R265), which was assembled to chromosome level previously.

```BWA-MEM``` was used to get the alignments and ```SAMtools``` was used to convert and sort the output bam files. Then ```Picard``` was used to remove the duplicates from the alignments. 

```{bash}
cd ..
mkdir mapping ; cd mapping
nano bwa.sh
batch bwa.sh
```

```{bash}
# bwa.sh script

#!/bin/bash

path_to_data=/scratch/youman7/reads/

individuals="JF109
CDC15
YM165"

for each_individual in $individuals

do
    # echo which individual we are on                                           
                        
    echo ${each_individual}

    # align the reads to the reference gemome                                    
                        
    bwa mem -t 16 -R "@RG\tID:LANE\tSM:${each_individual}\tPL:illumina" /scratch/youman7/reference/R265_VGII_NCBI_genomic.fna ${path_to_data}/${each_individual}_1.fastq ${path_to_data}/${each_individual}_2.fastq | samtools view -bSh - > ${each_individual}.bam

    # sort the bam file                                                         
                        
    samtools sort ${each_individual}.bam -o ${each_individual}_sorted.bam

    # index the sorted bam file                                                 
                        
    samtools index ${each_individual}_sorted.bam
    
    # counts the number of alignments for each FLAG type
    
    samtools flagstat ${each_individual}_sorted.bam
    
    # remove duplicates (rmdup)
    
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${path_to_data}/${each_individual}_sorted.bam O=${each_individual}_sorted.markdup.bam M=${each_individual}_markdup_metrics.txt

    # index the sorted rmdup bam files
    
    samtools index ${each_individual}_sorted.markdup.bam
    
    # counts the number of alignments for each FLAG type
    
    samtools flagstat ${each_individual}_sorted.markdup.bam > ${each_individual}_mappingstats.txt
              
done
```

## Alignment polishing

To find _de novo_ variations between parents and progeny, I thought using parents as references to call variants would be easier and more straightforward. So in addition, I used ```Pilon``` to polish the alignments of two parents. 

```Pilon``` uses read alignment analysis to identify inconsistencies between the input genome and the evidence in the reads. It can attempt to improve the input genome, including single-base differences, small indels, larger indel or block substitution events, gap filling, and identification of local misassemblies. It then generates a FASTA file containing an improved representation of the genome from the read data. The FASTA files were then used as reference genome to call variants later.

```{bash}
mkdir polish ; cd polish

# JF109 (Parent #1)

java -Xmx2G -jar $EBROOTPILON/pilon-1.23.jar --genome /scratch/youman7/reference/R265_VGII_NCBI_genomic.fna -- frags/scratch/youman7/mapping/JF109_sorted.markdup.bam --fix snps,indels --output JF109_pilon_polished --vcf &> JF109_sorted_pilon.log

# CDC15 (Parent #2)

java -Xmx2G -jar $EBROOTPILON/pilon-1.23.jar --genome /scratch/youman7/reference/R265_VGII_NCBI_genomic.fna --frags /scratch/youman7/mapping/CDC15_sorted.markdup.bam --fix snps,indels --output CDC15_pilon_polished --vcf &> CDC15_sorted_pilon.log

#change the file labels

cat JF109_pilon_polished.fasta | sed -e 's/>tig.*/ >genome_sequence/g' > JF109.fasta

cat CDC15_pilon_polished.fasta | sed -e 's/>tig.*/ >genome_sequence/g' > CDC15.fasta
```

## Annotation

```AUGUSTUS``` is a program that predicts genes in eukaryotic genomic sequences. Also, both _Cryptococcus gattii_ and _Cryptococcus neoformans_ are included in its database. So, I used ```AUGUSTUS``` to annotate the two parental polished alignments.

```{bash}
# JF109 (Parent #1) belongs to _Cryptococcus gattii_

augustus --progress=true --strand=both --species=cryptococcus_neoformans_gattii --protein=on --cds=on --codingseq=on JF109.fasta --gff3=on > JF109.gff

# CDC15 (Parent #2) belongs to _Cryptococcus neoformans_

augustus --progress=true --strand=both --species=cryptococcus_neoformans_neoformans_B --protein=on --cds=on --codingseq=on CDC15.fasta --gff3=on > CDC15.gff

#NOTE: the option of species is different. 
```

Now we have the alignments of each sample, we can next use them to call SNPs. Please click [HERE](https://github.com/manyou7/bio722_project_genomic_variations/blob/01e4810d1d0725e7a8bde9d12f16d0368fbc6666/3_variant_calling.md)


