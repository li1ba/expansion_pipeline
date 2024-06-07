# Repeat Expansion Detection Workflow for DNAnexus

This repository contains tools and a workflow for detecting short tandem repeat (STR) expansions using DNAnexus. The workflow has been tested using exome sequencing (ES) data. The workflow integrates the following tools:

- **ExpansionHunter**
- **Stranger**
- **Reviewer**

Additionally, a TSV viewer tool is included to assist in viewing and sorting workflow output.

## Table of Contents
1. [Introduction](#introduction)
2. [Prerequisites](#prerequisites)
3. [Setup](#setup)
    - [Running on DNAnexus](#running-on-dnanexus)
    - [Running Locally](#running-locally)
4. [TSV file Viewer](#tsv-viewer)
5. [Variant Catalog](#variant-catalog)
6. [Additional information](#additional-info)

## Introduction

This repository provides a workflow to detect repeat expansions (RE) in ES data. It uses several tools to analyze sequencing data and can be executed on the DNAnexus platform. The workflow pulls necessary Docker images, creates applet objects, and assembles them into a workflow object on DNAnexus.

## Prerequisites

Before using this workflow, ensure you have the following:
- DNAnexus account and API access
- DNAnexus CLI (`dx-toolkit`) installed
- Git installed

## Setup

Clone this repository and navigate to the directory:

```sh
git clone https://github.com/li1ba/expansion_pipeline.git
cd expansion_pipeline
```

### Building on DNAnexus
Make sure you are logged in DNAnexus:
```sh
dx login
```
Run the build script to set up the workflow:
```sh
bash build_workflow.sh
```
This script performs several functions:
- Pulls required Docker images.
- Creates Docker image files for each application.
- Generates applet objects in DNAnexus.
- Assembles the workflow into a workflow object in DNAnexus.

This workflow is configured for the **AWS region `eu-central-1`**. If your DNAnexus region is different, you must update the `regionalOptions` in each applet's `dxapp.json` file to match your designated region.

### Local Usage

To use the workflow locally as a bash script, run:
```sh
bash repeat_expansion_detection.sh -b <bamname> -f <fastaname> -v <varcat> -s <sex>
```
Replace `<bamname>`, `<fastaname>`, `<varcat>`, and `<sex>` with your specific file names and parameters.

-b: BAM file name (BAM index .bai should also be in the same directory as FASTA file)

-f: FASTA file name

-v: Variant catalog file

-s: Sample sex (e.g., male, female)


Output will be created in the same directory as the script. REViewer output will be created in a directory reviewer_output.

## TSV Viewer

Included in the repository is the `TSV file viewer` script for handling TSV file output in DNAnexus platform. This script allows for searching and sorting through detected STR loci in the output files.

To use `TSV file viewer` in DNAnexus, upload it to your project:
```sh
dx upload --type FileViewer --details='{"patterns":["*.tsv"]}' "TSV file viewer"
```

## Variant Catalog

A variant catalog in .json format is provided for human reference hg38 that covers REs captured in ES dataset. 

## Additional Information

For more detailed information about the tools, visit their respective GitHub pages:

- [ExpansionHunter GitHub](https://github.com/Illumina/ExpansionHunter)
- [Stranger GitHub](https://github.com/moonso/stranger)
- [Reviewer GitHub](https://github.com/Illumina/REViewer) 
