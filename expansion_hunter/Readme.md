# Expansion Hunter

## What does this app do?

This app implements Expansion Hunter tool (v 4.0.2) for estimating repeat sizes. Expansion Hunter aims to estimate sizes of such repeats by performing a targeted search through a BAM file for reads that span, flank, and are fully contained in each repeat. For more information about Expansion Hunter, please visit https://github.com/Illumina/ExpansionHunter. 

## What are typical use cases for this app?

Use the app to detect repeat expansions for previously defined loci for *one* sample in exome or genome sequencing data (WES/WGS).

## What data are required for this app to run?

Required:

- sample BAM file (.bam)
- sample BAM file index (.bai)
- FASTA file with reference genome
- FASTA file index
- Variant catalog defining repeat expansion loci

Optional:

- Sample sex (only affects loci on sex chromosomes)
- Region extension length (length of region around repeat to examine, default: 1000)
- Output JSON (default: False)


## What does this app output?

- BAMlet containing alignments of reads that overlap or located in close proximity to each variant
- VCF file with variant genotypes
- JSON file with variant genotypes (optional)

## How does this app work?

This app performs the following steps:
  - Runs ExpansionHunter tool (https://github.com/Illumina/REViewer) for one sample
  - Runs samtools to sort ExpansionHunter BAMlet

