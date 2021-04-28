---
Title: "Plots"
Author: "Man You"
Date: "12/04/2021"
Output: html_document
---

Back to [structural variants](https://github.com/manyou7/Bio722_project_genomic_variation/blob/429305d62aeffb9c0a8c3a63a39cf29a1ccd0707/5_structural_variants.md) or [README](https://github.com/manyou7/Bio722_project_genomic_variation/blob/73e10cbdcfe81b789a9d6b22b8bbdb2e336b3dda/README.md)

# Data Visualization

For SNPs and Structural variations, the HTML outputs of SnpEff already summarized and visualized the data. Although the SnpEff HTML outputs already showed the SNP density per 10KB individually, here, I showed the SNP density for all by chromosome. 

First, we need to library the packages that will be used.

```{r}
library(ggplot2)
library(multcompView)
library(RColorBrewer)
```

## Plot SNP density of all by chromosomes on one graph

```{r all, echo=FALSE}
# Set the working directory

setwd("~/R/Bio722_project_genomic_variation/snps/twocallers")

# Read data

YM165_JF109.snp <- read.table("~/R/Bio722_project_genomic_variation/snps/twocallers/YM165_JF109_unique_SNPs_densityPer10Kb.snpden", header=T)

YM165_CDC15.snp <- read.table("~/R/Bio722_project_genomic_variation/snps/twocallers/YM165_CDC15_unique_SNPs_densityPer10Kb.snpden", header=T)

CDC15_JF109.snp <- read.table("~/R/Bio722_project_genomic_variation/snps/twocallers/CDC15_JF109_unique_SNPs_densityPer10Kb.snpden", header=T)

# Add a new column called ID before merging them

YM165_JF109.snp$ID <- "YM165_JF109"

YM165_CDC15.snp$ID <- "YM165_CDC15"

CDC15_JF109.snp$ID <- "CDC15_JF109"

# Now, these are ready to be merged three datasets together

all <- rbind(YM165_JF109.snp, YM165_CDC15.snp,CDC15_JF109.snp)

ggplot(all, aes(x=BIN_START, y=SNP_COUNT, color=ID))+
  geom_jitter(size=1.5, stat = "identity", position = "jitter", alpha=0.5)+
  xlab('Positions')+
  ylab("SNP density per 10KB")+
  scale_color_brewer(palette="Dark2")+
  facet_wrap(~CHROM, scales = "free")+
  ggtitle("SNP density of unqiue SNPs between samples")+
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  theme(axis.title.x=element_text(face = "bold", size = 10))+
  theme(axis.title.y = element_text(face = "bold", size = 10))+
  theme(axis.text.x=element_text(size=7, color='black', face='bold', vjust=0.7))+
  theme(axis.text.y=element_text(size=7,face='bold'))
```

## Plot the nonsynonymous SNPs called by ```VarScan trio```

```{r trio_snp, echo=FALSE}
# Set the working directory

setwd("~/R/Bio722_project_genomic_variation/snps/trio")

# Read data

snp <- read.csv("snp_counts.csv")
nonsy.snp <- read.csv("snp_nosynonymous_counts.csv")
denovo.snp <- read.csv("denovo_mutation_snp_nonsynonymous_counts.csv")
MIE.snp <- read.csv("MIE_snp_nonsynonymous_counts.csv")
trans.snp <- read.csv("transmitted_snp_nonsynonymous_counts.csv")
untrans.snp <- read.csv("untransmitted_snp_nonsynonymous_counts.csv")

# Add a new column called TYPE

snp$TYPE <- "Total"
nonsy.snp$TYPE <- "Total_nonsyn"
denovo.snp$TYPE <- "DENOVO"
MIE.snp$TYPE <- "MIE"
trans.snp$TYPE <- "Transmitted"
untrans.snp$TYPE <- "Untransmitted"

# Merge them

trio_snp <- rbind(snp, nonsy.snp, denovo.snp, MIE.snp, trans.snp, untrans.snp)

# Now plot it. However, bar chart is not good for this

mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 7))

ggplot(data=trio_snp, aes(x=CHROM, y=COUNTS, fill=TYPE)) +
  geom_bar(stat="identity", position="stack",color="black")+
  xlab('Chromosomes')+
  ylab("Counts")+
  scale_color_manual(values = mycolors)+
  ggtitle("Different types of SNPs between progeny and parents)")

# dots are better than bar chart

ggplot(trio_snp, aes(x=CHROM, y=COUNTS, color=TYPE))+
  geom_jitter(size=1.5, stat = "identity", position = "jitter", alpha=0.5)+
  xlab('Chromosomes')+
  ylab("SNP counts")+
  scale_color_brewer(palette="Dark2")+
  ggtitle("Different types of SNPs between progeny and parents")+
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  theme(axis.title.x=element_text(face = "bold", size = 10))+
  theme(axis.title.y = element_text(face = "bold", size = 10))+
  theme(axis.text.x=element_text(size=7, color='black', face='bold', vjust=0.7))+
  theme(axis.text.y=element_text(size=7,face='bold'))

# Summarize the data in a table would make more sense.
```

## Plot _de novo_ SVs between progeny and each parent

```{r svs, echo=FALSE}
# Set the working directory

setwd("~/R/Bio722_project_genomic_variation/svs")

# read the data (the data were collected from the output files)

ym165_sv <- file.choose("YM165_denovo_sv.csv") #somehow read.csv can't read the csv file directly, so I did this additional step, and it worked

svs <- read.csv(ym165_sv)

# due to the length differences, so I convert the data to log10

svs$LENGTH <- log10(svs$LENGTH)

# check 

svs$LENGTH

# plot it without BND data (there was no length for BND)

ggplot(svs, aes(x=factor(CHROM, levels=unique(CHROM)), y=LENGTH, color=GROUP, shape=TYPE))+
  geom_jitter(size=3, stat = "identity", position = "jitter", alpha=0.7)+
  xlab('Chromosomes')+
  ylab("LOG10 SV length (bases)")+
  scale_color_brewer(palette="Dark2")+
  ggtitle("De novo SVs between progeny and each parent")+
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  theme(axis.title.x=element_text(face = "bold", size = 10))+
  theme(axis.title.y = element_text(face = "bold", size = 10))+
  theme(axis.text.x=element_text(size=7, color='black', face='bold', vjust=0.7))+
  theme(axis.text.y=element_text(size=7,face='bold'))
```

That's all.

Hope you enjoyed it and please feel free to comment on any steps you think need to be improved or include additional steps! 

Thank you!

