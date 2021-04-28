---
Title: "Variant calling"
Author: "Man You"
Date: "15/04/2021"
Output: html_document
---
Back to [Mapping and Polishing](2_mapping.md) or [README](README.md)

For variant calling, I used both the R265 reference genome and parental genomes as the reference. I decided to use parental genomes as references because I want to get _de novo_ SNPs between progeny and each parent. I think it should work, and even better if I use parental genomes as the reference to call variants. One of my concerns is that using reference genome to call variants individually and then comparing them would miss some information. For example, due to the differences between them and the reference genome, we might miss some variants. Directly using the parental genome would avoid it in a way.

## Variant calling using freebayes

I used ```GATK``` to call variants which worked for my trial round before, but it would never work again. So, I used ```freebayes``` instead. 

```Freebayes``` is a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment. Also, ```freebayes``` can set ploidy levels. In my case, the two parents are haploid, while the progeny is diploid. So, I called vairants individually with different ploidy settings. Here, I called variants using both R265 and parental genomes.

```{bash}
cd ..
mkdir freebayes ; cd freebayes

# use the reference R265 to call SNPs with quality > 20

freebayes -f /scratch/youman7/reference/R265_VGII_NCBI_genomic.fna -p 2 /scratch/youman7/mapping/YM165_sorted.markdup.bam | vcffilter -f "QUAL > 20"| vcffilter -f "TYPE= snp"> YM165_freebayes.vcf

freebayes -f /scratch/youman7/reference/R265_VGII_NCBI_genomic.fna -p 1 /scratch/youman7/mapping/CDC15_sorted.markdup.bam | vcffilter -f "QUAL > 20"| vcffilter -f "TYPE= snp"> CDC15_freebayes.vcf

freebayes -f /scratch/youman7/reference/R265_VGII_NCBI_genomic.fna -p 1 /scratch/youman7/mapping/JF109_sorted.markdup.bam | vcffilter -f "QUAL > 20"| vcffilter -f "TYPE= snp"> JF109_freebayes.vcf

# use JF109 (Parent #1) as the reference to call SNPs with quality > 20

freebayes -f /scratch/youman7/polish/JF109.fasta -p 2 /scratch/youman7/mapping/YM165_sorted.markdup.bam | vcffilter -f "QUAL > 20"| vcffilter -f "TYPE= snp"> YM165_JF109_freebayes.vcf

freebayes -f /scratch/youman7/polish/JF109.fasta -p 1 /scratch/youman7/mapping/CDC15_sorted.markdup.bam | vcffilter -f "QUAL > 20"| vcffilter -f "TYPE= snp"> CDC15_JF109_freebayes.vcf

# use CDC15 (Parent #2) as the reference to call SNPs with quality > 20

freebayes -f /scratch/youman7/polish/CDC15.fasta -p 2 /scratch/youman7/mapping/YM165_sorted.markdup.bam | vcffilter -f "QUAL > 20"| vcffilter -f "TYPE= snp"> YM165_CDC15_freebayes.vcf

freebayes -f /scratch/youman7/polish/CDC15.fasta -p 1 /scratch/youman7/mapping/CDC15_sorted.markdup.bam | vcffilter -f "QUAL > 20"| vcffilter -f "TYPE= snp"> JF109_CDC15_freebayes.vcf
```

## The second variant caller VarScan

I did not plan to use two callers in this project. When I was looking for tools for detecting _de novo_ mutations, I found ```VarScan```. Then I noticed it could also identify SNPs and INDELs. So I thought I would give it a shot; if it worked, I could merge the outputs from two callers. So I tried and it worked.

```{bash}
cd ..
mkdir varscan_snp ; cd varscan_snp

# Additional setting for Fisher's Exact Test p-value 0.05 as the threshold. The default minimum base quality is 20.

# use R265 as the reference

samtools mpileup -B -q 1 -f /scratch/youman7/reference/R265_VGII_NCBI_genomic.fna /scratch/youman7/mapping/CDC15_sorted.markdup.bam | java -Xmx2G -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar mpileup2snp --p-value 0.05 --output-vcf 1 > CDC15.varscan_snp.vcf

samtools mpileup -B -q 1 -f /scratch/youman7/reference/R265_VGII_NCBI_genomic.fna /scratch/youman7/mapping/JF109_sorted.markdup.bam | java -Xmx2G -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar mpileup2snp --p-value 0.05 --output-vcf 1 > JF109.varscan_snp.vcf

samtools mpileup -B -q 1 -f /scratch/youman7/reference/R265_VGII_NCBI_genomic.fna /scratch/youman7/mapping/YM165_sorted.markdup.bam | java -Xmx2G -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar mpileup2snp --p-value 0.05 --output-vcf 1 > YM165.varscan_snp.vcf

# use JF109 (Parent #1) as the reference

samtools mpileup -B -q 1 -f /scratch/youman7/polish/JF109.fasta /scratch/youman7/mapping/CDC15_sorted.markdup.bam | java -Xmx2G -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar mpileup2snp --p-value 0.05 --output-vcf 1 > CDC15_JF109.varscan_snp.vcf

samtools mpileup -B -q 1 -f /scratch/youman7/polish/JF109.fasta /scratch/youman7/mapping/YM165_sorted.markdup.bam | java -Xmx2G -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar mpileup2snp --p-value 0.05 --output-vcf 1 > YM165_JF109.varscan_snp.vcf

# use CDC15 (Parent #2) as the reference

samtools mpileup -B -q 1 -f /scratch/youman7/polish/CDC15.fasta /scratch/youman7/mapping/YM165_sorted.markdup.bam | java -Xmx2G -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar mpileup2snp --p-value 0.05 --output-vcf 1 > YM165_CDC15.varscan_snp.vcf

samtools mpileup -B -q 1 -f /scratch/youman7/polish/CDC15.fasta /scratch/youman7/mapping/JF109_sorted.markdup.bam | java -Xmx2G -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar mpileup2snp --p-value 0.05 --output-vcf 1 > CDC15_JF109.varscan_snp.vcf
```

## Compare and merge outputs from two callers

The two variant callers have their advantages. For example, ```freebayes``` has the ploidy option while ```VarScan``` does not. So, it is good to use multiple callers then merge the outputs. In this way, although we might lose some variants, we can get true and good-quality variants.

```BEDtools``` was used to compare and merge the variants from two callers, only the same variants were kept.

```SnpEff``` is a variant annotation and effect prediction tool. It annotates and predicts the effects of genetic variants (such as amino acid changes). The outputs also include a summary. So I used it for annotating, predicting and summarizing.

```{bash}
cd ..
mkdir snp ; cd snp

# get a symbolic link to freebayes outputs in this directory first 

ln -s /scratch/youman7/freebayes/*.vcf .

ln -s /scratch/youman7/varscan_snp/*.vcf .

# merge 

bedtools intersect -a YM165_freebayes.vcf -b YM165.varscan_snp.vcf -header > YM165_snp.vcf

bedtools intersect -a CDC15_freebayes.vcf -b CDC15.varscan_snp.vcf -header > CDC15_snp.vcf

bedtools intersect -a JF109_freebayes.vcf -b JF109.varscan_snp.vcf -header > JF109_snp.vcf

bedtools intersect -a YM165_JF109_freebayes.vcf -b YM165_JF109.varscan_snp.vcf -header > YM165_JF109_unique_snp.vcf

bedtools intersect -a CDC15_JF109_freebayes.vcf -b CDC15_JF109.varscan_snp.vcf -header > CDC15_JF109_snp.vcf

bedtools intersect -a YM165_CDC15_freebayes.vcf -b YM165_CDC15.varscan_snp.vcf -header > YM165_CDC15_unique_snp.vcf

bedtools intersect -a JF109_CDC15_freebayes.vcf -b JF109_CDC15.varscan_snp.vcf -header > JF109_CDC15_snp.vcf

# merge two parental snps files

bedtools intersect -a JF109_CDC15_snp.vcf -b CDC15_JF109_snp.vcf -header > JF109_CDC15_all_snp.vcf

# annotate and summarize

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 YM165_JF109_unique_snp.vcf > YM165_JF109_snp_ann.vcf

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 YM165_CDC15_unique_snp.vcf > YM165_CDC15_snp_ann.vcf

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 CDC15_JF109_unique_snp.vcf > CDC15_JF109_snp_ann.vcf
```

#### Annotation outputs
- [YM165_JF109_unique_snp_snpEff_summary.html](https://rpubs.com/manyou7/YM165_JF109_unique_snp)
- [YM165_CDC15_unique_snp_snpEff_summary.html](https://rpubs.com/manyou7/YM165_CDC15_unique_snp)
- [CDC15_JF109_unique_snp_snpEff_summary.html](https://rpubs.com/manyou7/CDC15_JF109_unique_snp_SnpEff)

##### NOTE: I downloaded the output HTML files from the cluster to my local laptop, however, the HTML files are local files in a ```files:///``` format. So I copied all codes of each local HTML file to a R HTML file and them get them published.

## SNP density

Now let's get the SNP density using ```vcftools```. Here, I calculated the SNP density per 10K bases

```{bash}
vcftools --vcf YM165_JF109_snp.vcf --SNPdensity 10000 --out YM165_JF109_SNPS_densityPer10Kb

vcftools --vcf YM165_CDC15_snp.vcf --SNPdensity 10000 --out YM165_CDC15_SNPS_densityPer10Kb

vcftools --vcf JF109_CDC15_snp.vcf --SNPdensity 10000 --out JF109_CDC15_SNPS_densityPer10Kb
```

The scripts for generating plots are in 6_plots. But let's go see [denovo mutations](4_denovo_mutations.md) first. 






