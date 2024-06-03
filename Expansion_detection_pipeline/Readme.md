# Expansion Detection Pipeline

## Overview

The Expansion Detection Pipeline is designed to identify and analyze short tandem repeat (STR) expansions in NGS data. The pipeline integrates three main tools: **ExpansionHunter**, **Stranger**, and **REViewer**.

## Tools Included

1. **ExpansionHunter**: Used for detecting and genotyping short tandem repeats (STRs) from sequence data (WES or WGS).
2. **Stranger**: A tool to interpret and annotate STR expansions.
3. **REViewer**: Provides a graphical review of the detected STRs, aiding in validation and visualization.

## Pipeline Stages

### 1. EH_RAW Stage

This stage involves the initial processing of genomic data to prepare it for STR expansion detection using **ExpansionHunter**. Key outputs include a variant catalog and a BAM file with corresponding index.

### 2. Stranger Stage

This stage processes the variant catalog and VCF file using **Stranger** to detect STR expansions.

### 3. Reviewer Stage

This stage involves reviewing the processed data for final verification and validation of STR expansions using **REViewer**.

More information about input files and output files of each stage can be found in README docs of each app.
