---
Title: "Structural Variants"
Author: "Man You"
Date: "20/04/2021"
Output: html_document
---
Back to [_de novo_ mutations](4_denovo_mutations.md) or [README](README.md)

Many tools could be used to detect structural variants, such as ```LUMPY```, ```BreakDancer``` and ```DELLY```. 

For this project, I chose ```LUMPY``` for trial. It did not work for a while. So I decided to identify SNPs (in 3_variant_calling) and _de novo_ mutations (4_denovo_mutations), as well as copy number variants (not included in this project). 

However, just recently, after searching for solutions for troubleshooting over and over again, I finally figured it out. The main reason was due to the version of python. ```LUMPY``` only works under python 2 not python 3, as well as some required dependencies. And some other minor problems. Then it finally worked! But due to time, I was not able to look into other tools, but for my research, I will use at least two tools to detect the structural variants and combine the results for analysis.

To detect _de novo_ structural variations between progeny and each parent, I decided to use the parental genome as the reference to call variants.

```LUMPY``` recommended two tools, ```SAMtools``` and ```SpeedSeq```, to get splitters and discordants. By comparing the two tools, ```SpeedSeq``` is more straightforward, easier, and faster. So I decided to use ```SpeedSeq``` for this project. 

## Structural variants between YM165 and JF109 (Parent #1)

```{bash}
cd ..
mkdir sv ; cd sv

# Use SpeedSeq to align reads and get splitter and discordant bam files

speedseq align -t 8 -R "@RG\tID:id\tSM:YM165\tLB:lib" "/scratch/youman7/polish/JF109.fasta" /scratch/youman7/reads/YM165_1.fastq /scratch/youman7/reads/YM165_2.fastq

# Run LUMPY

/scratch/youman7/lumpy-sv/scripts/lumpyexpress -B YM165_1.fastq.bam -S YM165_1.fastq.splitters.bam -D YM165_1.fastq.discordants.bam -o YM165_JF109_lumpy.vcf

# Use SVTyper to call genotypes on LUMPY output VCF files. SVTyper is a recommended tool by LUMPY.

svtyper -B YM165_1.fastq.bam -S YM165_1.fastq.splitters.bam -i YM165_JF109_lumpy.vcf > YM165_JF109_lumpy.gt.vcf

# Because LUMPY doesn't have a quality filter option, so additionally, I used SnpSift to do the quality control

cat YM165_JF109_lumpy.gt.vcf | java -Xmx2G -jar "/scratch/youman7/snpEff/SnpSift.jar" filter "(QUAL >= 10 )" > YM165_JF109_lumpy.gt.filtered.vcf

# relabeled the chromosomes

sed -i 's/CP025759.1_pilon/CP025759.1/g;s/CP025760.1_pilon/CP025760.1/g;s/CP025761.1_pilon/CP025761.1/g;s/CP025762.1_pilon/CP025762.1/g;s/CP025763.1_pilon/CP025763.1/g;s/CP025764.1_pilon/CP025764.1/g;s/CP025765.1_pilon/CP025765.1/g;s/CP025766.1_pilon/CP025766.1/g;s/CP025767.1_pilon/CP025767.1/g;s/CP025768.1_pilon/CP025768.1/g;s/CP025769.1_pilon/CP025769.1/g;s/CP025770.1_pilon/CP025770.1/g;s/CP025771.1_pilon/CP025771.1/g;s/CP025772.1_pilon/CP025772.1/g;s/CP025773.1_pilon/CP025773.1/g' YM165_JF109_lumpy.gt.filtered.vcf

# Survivor_ant is a tool to annotate SVs (VCF file) with previous known SVs (VCF file) and or with genomic features (gff and or bed files). The default settings of maximum distance to assign annotation to SVs and maximum distance to group SV together are both 1000.

"/scratch/youman7/SURVIVOR_ant/bin/survivor_ant-core-0.1.0/survivor_ant" -o YM165_JF109_lumpy_ann.vcf -g /scratch/youman7/reference/R265_chr_genomic.gff -i YM165_JF109_lumpy.gt.filtered.vcf

# Annotate and summarize the output using SnpEff

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 YM165_JF109_lumpy_ann.vcf > YM165_JF109_sv.vcf
```

## Structural variants between YM165 and CDC15 (Parent #2)

```{bash}
# align reads and get splitter and discordant bam files

speedseq align -t 8 -R "@RG\tID:id\tSM:YM165\tLB:lib" "/scratch/youman7/polish/CDC15.fasta" /scratch/youman7/reads/YM165_1.fastq /scratch/youman7/reads/YM165_2.fastq

# run LUMPY

/scratch/youman7/lumpy-sv/scripts/lumpyexpress -B YM165_1.fastq.bam -S YM165_1.fastq.splitters.bam -D YM165_1.fastq.discordants.bam -o YM165_CDC15_lumpy.vcf

# call genotypes 

svtyper -B YM165_1.fastq.bam -S YM165_1.fastq.splitters.bam -i YM165_CDC15_lumpy.vcf > YM165_CDC15_lumpy.gt.vcf

# quality control

cat YM165_CDC15_lumpy.gt.vcf | java -Xmx2G -jar "/scratch/youman7/snpEff/SnpSift.jar" filter "(QUAL >= 10 )" > YM165_CDC15_lumpy.gt.filtered.vcf

# relabeled the chromosomes

sed -i 's/CP025759.1_pilon/CP025759.1/g;s/CP025760.1_pilon/CP025760.1/g;s/CP025761.1_pilon/CP025761.1/g;s/CP025762.1_pilon/CP025762.1/g;s/CP025763.1_pilon/CP025763.1/g;s/CP025764.1_pilon/CP025764.1/g;s/CP025765.1_pilon/CP025765.1/g;s/CP025766.1_pilon/CP025766.1/g;s/CP025767.1_pilon/CP025767.1/g;s/CP025768.1_pilon/CP025768.1/g;s/CP025769.1_pilon/CP025769.1/g;s/CP025770.1_pilon/CP025770.1/g;s/CP025771.1_pilon/CP025771.1/g;s/CP025772.1_pilon/CP025772.1/g;s/CP025773.1_pilon/CP025773.1/g' YM165_CDC15_lumpy.gt.filtered.vcf

# annotate SV files

"/scratch/youman7/SURVIVOR_ant/bin/survivor_ant-core-0.1.0/survivor_ant" -o YM165_CDC15_lumpy_ann.vcf -g /scratch/youman7/reference/R265_chr_genomic.gff -i YM165_CDC15_lumpy.gt.filtered.vcf

# annotate and summarize the output

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 YM165_CDC15_lumpy_ann.vcf > YM165_CDC15_sv.vcf
```

## Structural variants between two parents using CDC15 (Parent #2) as the reference

```{bash}
# align reads and get splitter and discordant bam files

speedseq align -t 8 -R "@RG\tID:id\tSM:JF109\tLB:lib" "/scratch/youman7/polish/CDC15.fasta" /scratch/youman7/reads/JF109_1.fastq /scratch/youman7/reads/JF109_2.fastq

# run LUMPY

/scratch/youman7/lumpy-sv/scripts/lumpyexpress -B JF109_1.fastq.bam -S JF109_1.fastq.splitters.bam -D JF109_1.fastq.discordants.bam -o JF109_CDC15_lumpy.vcf

# call genotypes 

svtyper -B JF109_1.fastq.bam -S JF109_1.fastq.splitters.bam -i JF109_CDC15_lumpy.vcf > JF109_CDC15_lumpy.gt.vcf

# quality control

cat JF109_CDC15_lumpy.gt.vcf | java -Xmx2G -jar "/scratch/youman7/snpEff/SnpSift.jar" filter "(QUAL >= 10 )" > JF109_CDC15_lumpy.gt.filtered.vcf

# relabeled the chromosomes

sed -i 's/CP025759.1_pilon/CP025759.1/g;s/CP025760.1_pilon/CP025760.1/g;s/CP025761.1_pilon/CP025761.1/g;s/CP025762.1_pilon/CP025762.1/g;s/CP025763.1_pilon/CP025763.1/g;s/CP025764.1_pilon/CP025764.1/g;s/CP025765.1_pilon/CP025765.1/g;s/CP025766.1_pilon/CP025766.1/g;s/CP025767.1_pilon/CP025767.1/g;s/CP025768.1_pilon/CP025768.1/g;s/CP025769.1_pilon/CP025769.1/g;s/CP025770.1_pilon/CP025770.1/g;s/CP025771.1_pilon/CP025771.1/g;s/CP025772.1_pilon/CP025772.1/g;s/CP025773.1_pilon/CP025773.1/g' JF109_CDC15_lumpy.gt.filtered.vcf

# annotate SV files

"/scratch/youman7/SURVIVOR_ant/bin/survivor_ant-core-0.1.0/survivor_ant" -o JF109_CDC15_lumpy_ann.vcf -g /scratch/youman7/reference/R265_chr_genomic.gff -i JF109_CDC15_lumpy.gt.filtered.vcf

# annotate and summarize the output

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 JF109_CDC15_lumpy_ann.vcf > JF109_CDC15_sv.vcf
```

## Find the commonly occurred SVs between two parents and progeny

```{bash}
# find common SVs

bedtools intersect -a YM165_JF109_sv.vcf -b YM165_CDC15_sv.vcf -header > YM165_common_sv.vcf

# relabel

sed -i 's/CP025759.1_pilon/CP025759.1/g;s/CP025760.1_pilon/CP025759.1/g;s/CP025761.1_pilon/CP025759.1/g;s/CP025762.1_pilon/CP025759.1/g;s/CP025763.1_pilon/CP025759.1/g;s/CP025764.1_pilon/CP025759.1/g;s/CP025765.1_pilon/CP025759.1/g;s/CP025766.1_pilon/CP025766.1/g;s/CP025767.1_pilon/CP025767.1/g;s/CP025768.1_pilon/CP025768.1/g;s/CP025769.1_pilon/CP025769.1/g;s/CP025770.1_pilon/CP025770.1/g;s/CP025771.1_pilon/CP025771.1/g;s/CP025772.1_pilon/CP025772.1/g;s/CP025773.1_pilon/CP025773.1/g' YM165_common_sv.vcf

# annotate and summarize

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 YM165_common_sv.vcf > YM165_common_sv_ann.vcf

# get each SV type in an individual file

egrep '(^#|SVTYPE=DEL)' YM165_common_sv.vcf > YM165_common_sv_DEL.vcf #deletion

egrep '(^#|SVTYPE=DUP)' YM165_common_sv.vcf > YM165_common_sv_DUP.vcf #duplication

egrep '(^#|SVTYPE=INV)' YM165_common_sv.vcf  > YM165_common_sv_INV.vcf #inversion

egrep '(^#|SVTYPE=BND)' YM165_common_sv.vcf > YM165_common_sv_BND.vcf #breakend notations

egrep '(^#|SVTYPE=INS)' YM165_common_sv.vcf > YM165_common_sv_INS.vcf #insertion of novel sequence #none detected

egrep '(^#|SVTYPE=CNV)' YM165_common_sv.vcf > YM165_common_sv_CNV.vcf #copy number variable region #none detected

# annotate and summarize

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 YM165_common_sv_DEL.vcf > YM165_common_sv_DEL_ann.vcf

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 YM165_common_sv_DUP.vcf > YM165_common_sv_DUP_ann.vcf

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 YM165_common_sv_INV.vcf > YM165_common_sv_INV_ann.vcf

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 YM165_common_sv_BND.vcf > YM165_common_sv_BND_ann.vcf
```

## Find _de novo_ SVs between YM165 and JF109 (Parent #1)

```{bash}
# find _de novo_ SVs

bedtools intersect -a YM165_JF109_sv.vcf -b YM165_CDC15_sv.vcf -v -header > YM165_to_JF109_unique_sv.vcf

# relabel

sed -i 's/CP025759.1_pilon/CP025759.1/g;s/CP025760.1_pilon/CP025759.1/g;s/CP025761.1_pilon/CP025759.1/g;s/CP025762.1_pilon/CP025759.1/g;s/CP025763.1_pilon/CP025759.1/g;s/CP025764.1_pilon/CP025759.1/g;s/CP025765.1_pilon/CP025759.1/g;s/CP025766.1_pilon/CP025766.1/g;s/CP025767.1_pilon/CP025767.1/g;s/CP025768.1_pilon/CP025768.1/g;s/CP025769.1_pilon/CP025769.1/g;s/CP025770.1_pilon/CP025770.1/g;s/CP025771.1_pilon/CP025771.1/g;s/CP025772.1_pilon/CP025772.1/g;s/CP025773.1_pilon/CP025773.1/g' YM165_to_JF109_unique_sv.vcf

# annotate and summarize

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 YM165_to_JF109_unique_sv.vcf > YM165_to_JF109_unique_sv.vcf

# get each SV type in an individual file

egrep '(^#|SVTYPE=DEL)' YM165_to_JF109_unique_sv.vcf > YM165_to_JF109_unique_sv_DEL.vcf

egrep '(^#|SVTYPE=DUP)' YM165_to_JF109_unique_sv.vcf > YM165_to_JF109_unique_sv_DUP.vcf 

egrep '(^#|SVTYPE=INV)' YM165_to_JF109_unique_sv.vcf > YM165_to_JF109_unique_sv_INV.vcf

egrep '(^#|SVTYPE=BND)' YM165_to_JF109_unique_sv.vcf > YM165_to_JF109_unique_sv_BND.vcf

egrep '(^#|SVTYPE=INS)' YM165_to_JF109_unique_sv.vcf > YM165_to_JF109_unique_sv_INS.vcf #none detected

egrep '(^#|SVTYPE=CNV)' YM165_to_JF109_unique_sv.vcf > YM165_to_JF109_unique_sv_CNV.vcf #none detected

# annotate and summarize

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 YM165_to_JF109_unique_sv_DEL.vcf > YM165_to_JF109_unique_DEL_ann.vcf

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 YM165_to_JF109_unique_sv_DUP.vcf > YM165_to_JF109_unique_DUP_ann.vcf

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 YM165_to_JF109_unique_sv_INV.vcf > YM165_to_JF109_unique_INV_ann.vcf

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 YM165_to_JF109_unique_sv_BND.vcf > YM165_to_JF109_unique_BND_ann.vcf
```

## Find _de novo_ SVs between YM165 and CDC15 (Parent #2)

```{bash}
# find _de novo_ SVs

bedtools intersect -a YM165_CDC15_sv.vcf -b YM165_JF109_sv.vcf -v -header > YM165_to_CDC15_unique_sv.vcf

# relabel

sed -i 's/CP025759.1_pilon/CP025759.1/g;s/CP025760.1_pilon/CP025759.1/g;s/CP025761.1_pilon/CP025759.1/g;s/CP025762.1_pilon/CP025759.1/g;s/CP025763.1_pilon/CP025759.1/g;s/CP025764.1_pilon/CP025759.1/g;s/CP025765.1_pilon/CP025759.1/g;s/CP025766.1_pilon/CP025766.1/g;s/CP025767.1_pilon/CP025767.1/g;s/CP025768.1_pilon/CP025768.1/g;s/CP025769.1_pilon/CP025769.1/g;s/CP025770.1_pilon/CP025770.1/g;s/CP025771.1_pilon/CP025771.1/g;s/CP025772.1_pilon/CP025772.1/g;s/CP025773.1_pilon/CP025773.1/g' YM165_to_CDC15_unique_sv.vcf

# annotate and summarize

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 YM165_to_CDC15_unique_sv.vcf > YM165_to_CDC15_unique_sv.vcf

# get each SV type in an individual file

egrep '(^#|SVTYPE=DEL)' YM165_to_CDC15_unique_sv.vcf > YM165_to_CDC15_unique_sv_DEL.vcf

egrep '(^#|SVTYPE=DUP)' YM165_to_CDC15_unique_sv.vcf > YM165_to_CDC15_unique_sv_DUP.vcf 

egrep '(^#|SVTYPE=INV)' YM165_to_CDC15_unique_sv.vcf > YM165_to_CDC15_unique_sv_INV.vcf

egrep '(^#|SVTYPE=BND)' YM165_to_CDC15_unique_sv.vcf > YM165_to_CDC15_unique_sv_BND.vcf

egrep '(^#|SVTYPE=INS)' YM165_to_CDC15_unique_sv.vcf > YM165_to_CDC15_unique_sv_INS.vcf #none detected

egrep '(^#|SVTYPE=CNV)' YM165_to_CDC15_unique_sv.vcf > YM165_to_CDC15_unique_sv_CNV.vcf #none detected

# annotate and summarize

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 YM165_to_CDC15_unique_sv_DEL.vcf > YM165_to_CDC15_unique_DEL_ann.vcf

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 YM165_to_CDC15_unique_sv_DUP.vcf > YM165_to_CDC15_unique_DUP_ann.vcf

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 YM165_to_CDC15_unique_sv_INV.vcf > YM165_to_CDC15_unique_INV_ann.vcf

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 YM165_to_CDC15_unique_sv_BND.vcf > YM165_to_CDC15_unique_BND_ann.vcf
```

## Find SVs between two parents

```{bash}
# find SVs between parents

bedtools intersect -a CDC15_JF109_lumpy.vcf -b JF109_CDC15_lumpy.vcf -header > parents_common_sv.vcf

# relabel

sed -i 's/CP025759.1_pilon/CP025759.1/g;s/CP025760.1_pilon/CP025759.1/g;s/CP025761.1_pilon/CP025759.1/g;s/CP025762.1_pilon/CP025759.1/g;s/CP025763.1_pilon/CP025759.1/g;s/CP025764.1_pilon/CP025759.1/g;s/CP025765.1_pilon/CP025759.1/g;s/CP025766.1_pilon/CP025766.1/g;s/CP025767.1_pilon/CP025767.1/g;s/CP025768.1_pilon/CP025768.1/g;s/CP025769.1_pilon/CP025769.1/g;s/CP025770.1_pilon/CP025770.1/g;s/CP025771.1_pilon/CP025771.1/g;s/CP025772.1_pilon/CP025772.1/g;s/CP025773.1_pilon/CP025
773.1/g' parents_common_sv.vcf

# annotate and summarize

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 parents_common_sv.vcf > parents_sv_merged.vcf

# get each SV type in an individual file

egrep '(^#|SVTYPE=DEL)' parents_common_sv.vcf> parents_sv_DEL.vcf

egrep '(^#|SVTYPE=DUP)' parents_common_sv.vcf> parents_sv_DUP.vcf

egrep '(^#|SVTYPE=INV)' parents_common_sv.vcf > parents_sv_INV.vcf

egrep '(^#|SVTYPE=BND)' parents_common_sv.vcf> parents_sv_BND.vcf

egrep '(^#|SVTYPE=INS)' parents_common_sv.vcf> parents_sv_INS.vcf #none detected

egrep '(^#|SVTYPE=CNV)' parents_common_sv.vcf> parents_sv_CNV.vcf #none detected

# annotate and summarize

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 parents_sv_DEL.vcf > parents_sv_DEL_ann.vcf

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 parents_sv_DUP.vcf > parents_sv_DUP_ann.vcf

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 parents_sv_INV.vcf > parents_sv_INV_ann.vcf

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" R265 parents_sv_BND.vcf > parents_sv_BND_ann.vcf
```

Summaries obtained from the HTML outputs by ```SnpEff```.

##### Annotation outputs
- [YM165_to_CDC15_unique_sv_SnpEff_summary](https://rpubs.com/manyou7/YM165_to_CDC15_unique_sv)
- [YM165_to_JF109_unique_sv_SnpEff_summary](https://rpubs.com/manyou7/YM165_to_JF109_unique_sv)
- [YM165_common_sv_SnpEff_summary](https://rpubs.com/manyou7/YM165_common_sv)
- [Parents_sv_SnpEff_summary](https://rpubs.com/manyou7/parents_common_sv)

Please click [HERE](6_plots.md) to visualize the data.


