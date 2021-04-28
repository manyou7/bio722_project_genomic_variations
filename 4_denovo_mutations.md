---
Title: "_de novo_ mutations"
Author: "Man You"
Date: "10/04/2021"
Output: html_document
---
Back to [variant calling](https://github.com/manyou7/Bio722_project_genomic_variation/blob/346ee8842927b17bb02dea023d0c9673bbc3098f/3_variant_calling.md) or [README](https://github.com/manyou7/Bio722_project_genomic_variation/blob/73e10cbdcfe81b789a9d6b22b8bbdb2e336b3dda/README.md)

I found ```VarScan trio``` is a useful tool that can detect de novo mutations, as well as identify transmitted alleles. It can leverage the family relationship to improve variant calling accuracy, identify apparent Mendelian Inheritance Errors (MIEs), and detect high-confidence _de novo_ mutations. Although it was designed for human genome, I think it should work for fungi as well. The command requires "mpileup" for the father, mother and child in this order. So, I treated the two parents, JF109 as the father and CDC15 as the mother, based on their mating types, and of course, the progeny as the child. 

Then I ran:

```{bash}
cd..
mkdir trio_snp ; cd trio_snp

# Generate a three-sample mpileup

samtools mpileup -B -q 1 -f /scratch/youman7/reference/R265_chr_genomic.fna JF109.bam CDC15.bam YM165.bam >trio.mpileup

# Run VarScan trio with the recommended settings

java -Xmx2G -jar $EBROOTVARSCAN/VarScan.v2.4.2.jar trio trio.mpileup trio.mpileup.output --min-coverage 20 --min-var-freq 0.20 --p-value 0.05 -adj-var-freq 0.05 -adj-p-value 0.15
```

This command produced two VCF output files: one for SNPs and one for indels. In the outputs, each SNP/INDEL was assigned a STATUS, 1=untransmitted, 2=transmitted, 3=denovo, 4=MIE. 

Here, I only annotated and summarized SNPs. INDELs will be included in my research.

```{bash}
# Annotate the output as well as creating a section counting each effect by adding -csvStats

java -Xmx2G -jar "/scratch/youman7/snpEff/snpEff.jar" eff -csvStats snpEffoutput.csv R265 trio.mpileup.output.snp.vcf  > snpEffoutput.vcf

# Get summary of counts of each effect

grep -A24 "# Count by effects" snpEffoutput.csv
```

Nonsynonymous SNPs can change the amino acid sequence of protein, while the synonymous SNPs don't. So, I am more interested in looking into the nonsynonymous ones.

Then, I used  ```SnpSift filter``` to extract nonsynonymous SNPs from the overall SNPs. 

```{bash}
cat snpEffoutput.vcf | java -Xmx2G -jar "/scratch/youman7/snpEff/SnpSift.jar" filter "( EFF[*].EFFECT = 'missense_variant' )"  > trio.nosynonymous.snp.vcf
```

Then I separated files based on the STATUS:

```{bash}
#split different types 

egrep '(^#|STATUS=3)' trio.nosynonymous.snp.vcf > denovo_mutation_snp_nonsynonymous.vcf

egrep '(^#|STATUS=1)' trio.nosynonymous.snp.vcf > untransmitted_snp_nonsynonymous.vcf

egrep '(^#|STATUS=2)' trio.nosynonymous.snp.vcf > transmitted_snp_nonsynonymous.vcf

egrep '(^#|STATUS=4)' trio.nosynonymous.snp.vcf > MIE_snp_nonsynonymous.vcf
```

Next, I counted numbers of SNPs per chromosome

```{bash}
# total counts of nonsynonymous SNPs

grep -v "^#" trio.nosynonymous.snp.vcf | cut -f 1 | sort | uniq -c > snp_nosynonymous_counts.csv

# counts of transmitted nonsynonymous SNPs

grep -v "^#" transmitted_snp_nonsynonymous.vcf | cut -f 1 | sort | uniq -c > transmitted_snp_nonsynonymous_counts.csv

# counts of untransmitted nonsynonymous SNPs

grep -v "^#" untransmitted_snp_nonsynonymous.vcf | cut -f 1 | sort | uniq -c > untransmitted_snp_nonsynonymous_counts.csv

# counts of MIE nonsynonymous SNPs

grep -v "^#" MIE_snp_nonsynonymous.vcf | cut -f 1 | sort | uniq -c > MIE_snp_nonsynonymous_counts.csv

# counts of de novo nonsynonymous SNPs

grep -v "^#" denovo_mutation_snp_nonsynonymous.vcf | cut -f 1 | sort | uniq -c > denovo_mutation_snp_nonsynonymous_counts.csv
```

To make things easier, I relabeled all chromosomes.

```{bash}
sed -i 's/CP025759.1/Chr1/g;s/CP025760.1/Chr2/g;s/CP025761.1/Chr3/g;s/CP025762.1/Chr4/g;s/CP025763.1/Chr5/g;s/CP025764.1/Chr6/g;s/CP025765.1/Chr7/g;s/CP025766.1/Chr8/g;s/CP025767.1/Chr9/g;s/CP025768.1/Chr10/g;s/CP025769.1/Chr11/g;s/CP025770.1/Chr12/g;s/CP025771.1/Chr13/g;s/CP025772.1/Chr14/g;s/CP025773.1/MT/g' snp_nosynonymous_counts.csv

sed -i 's/CP025759.1/Chr1/g;s/CP025760.1/Chr2/g;s/CP025761.1/Chr3/g;s/CP025762.1/Chr4/g;s/CP025763.1/Chr5/g;s/CP025764.1/Chr6/g;s/CP025765.1/Chr7/g;s/CP025766.1/Chr8/g;s/CP025767.1/Chr9/g;s/CP025768.1/Chr10/g;s/CP025769.1/Chr11/g;s/CP025770.1/Chr12/g;s/CP025771.1/Chr13/g;s/CP025772.1/Chr14/g;s/CP025773.1/MT/g' transmitted_snp_nonsynonymous_counts.csv

sed -i 's/CP025759.1/Chr1/g;s/CP025760.1/Chr2/g;s/CP025761.1/Chr3/g;s/CP025762.1/Chr4/g;s/CP025763.1/Chr5/g;s/CP025764.1/Chr6/g;s/CP025765.1/Chr7/g;s/CP025766.1/Chr8/g;s/CP025767.1/Chr9/g;s/CP025768.1/Chr10/g;s/CP025769.1/Chr11/g;s/CP025770.1/Chr12/g;s/CP025771.1/Chr13/g;s/CP025772.1/Chr14/g;s/CP025773.1/MT/g' untransmitted_snp_nonsynonymous_counts.csv

sed -i 's/CP025759.1/Chr1/g;s/CP025760.1/Chr2/g;s/CP025761.1/Chr3/g;s/CP025762.1/Chr4/g;s/CP025763.1/Chr5/g;s/CP025764.1/Chr6/g;s/CP025765.1/Chr7/g;s/CP025766.1/Chr8/g;s/CP025767.1/Chr9/g;s/CP025768.1/Chr10/g;s/CP025769.1/Chr11/g;s/CP025770.1/Chr12/g;s/CP025771.1/Chr13/g;s/CP025772.1/Chr14/g;s/CP025773.1/MT/g' MIE_snp_nonsynonymous_counts.csv

sed -i 's/CP025759.1/Chr1/g;s/CP025760.1/Chr2/g;s/CP025761.1/Chr3/g;s/CP025762.1/Chr4/g;s/CP025763.1/Chr5/g;s/CP025764.1/Chr6/g;s/CP025765.1/Chr7/g;s/CP025766.1/Chr8/g;s/CP025767.1/Chr9/g;s/CP025768.1/Chr10/g;s/CP025769.1/Chr11/g;s/CP025770.1/Chr12/g;s/CP025771.1/Chr13/g;s/CP025772.1/Chr14/g;s/CP025773.1/MT/g' denovo_mutation_snp_nonsynonymous_counts.csv
```

Scripts of Plots are in 6_plots.Rmd. But now let's go check [structural variants detection](https://github.com/manyou7/Bio722_project_genomic_variation/blob/346ee8842927b17bb02dea023d0c9673bbc3098f/5_structural_variants.md).



