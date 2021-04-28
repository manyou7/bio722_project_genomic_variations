# Bio722 project: Genomic variation between parents and progeny in the human pathogenic _Cryptococcus_ speices complex

A variety of genetic changes can affect the biology of species, including single-nucleotide polymorphisms (SNPs), small insertion-deletion events (INDELs), transposon insertions and large structural variations (SVs). For example, SVs can contribute to the evolution of an organism through disruption of an existing gene, creation of a new gene, or a chimeric gene product through gene fusions. 

The human pathogenic _Cryptococcus_ species complex is the causative agent of cryptococcosis. It consists of two divergent lineages, _Cryptococcus neoformans_ and _Cryptococcus gattii_. Each lineage has multiple molecular types/subspecies with a range of nucleotide sequence divergence from 2% to 15%. Previously, chromosomal rearrangements, segmental duplications and whole-genome copy number variations have been described for _C. neoformans_. Recently, genome variation between some molecular types of _C. gattii_ has been reported. However, the genomic variations between parental strains and their hybrid progeny remain unknown. 

I have a total of seven crosses constructed from my previous study, and each of them includes two parents and one progeny. All samples were whole-genome sequenced on Illumina HiSeq X with PE 150 bp. My research aims to investigate the genomic variations between parental strains and their progeny in the human pathogenic _Cryptococcus_ species complex. 

For this Bio722 project, I chose one cross (two parents and one progeny), as an example to figure out the pipeline. 

Here show all steps and scripts for the pipeline, including pre-processing steps (QC and trimming), read mapping and polishing, SNPs calling and classcifying, structural variation detection, functional annotation, and some plots.

### Pipeline:
1. [Pre-processing](https://github.com/manyou7/bio722_project_genomic_variations/blob/3361a74fb74382bfb4f92567513722da75cf1b3d/1_pre_processing.md)
2. [Mapping and Polishing](https://github.com/manyou7/bio722_project_genomic_variations/blob/d7c37105a7ac40b83e0206fc979dce817ba4a9db/2_mapping.md)
3. [Variant calling](https://github.com/manyou7/bio722_project_genomic_variations/blob/01e4810d1d0725e7a8bde9d12f16d0368fbc6666/3_variant_calling.md)
4. [_De novo_ mutations calling](https://github.com/manyou7/bio722_project_genomic_variations/blob/605838e20dc4f866779096ea812baea25fa72b95/4_denovo_mutations.md)
5. [Structural variants detection](https://github.com/manyou7/bio722_project_genomic_variations/blob/01e4810d1d0725e7a8bde9d12f16d0368fbc6666/5_structural_variants.md)
6. [Data visualization](https://github.com/manyou7/bio722_project_genomic_variations/blob/01e4810d1d0725e7a8bde9d12f16d0368fbc6666/6_plots.md)

#### NOTE: Links might not work if the branch is not showing as main. I was not able to figure out why the branch name sometimes changes to a random name and the links won't work. But manually changing it back to main branch, all links will work again in each .md file.
