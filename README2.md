# PROJECT NAME

:collision: Brief description

# Contents

- [Summary](#summary)
- [Experimental Information](#experimental-information)
- [Contrasts](#contrasts)
- [Workflow](#)
    - [Step 00. Sample information file](#step-00-sample-information-file)
    - [Step 0. Fastq Quality Control](#step-0-fastq-quality-control)
    - [Step 1. Trimming](#step-1-trimming)
    - [Step 2. Fastq Quality Control](#step-2-fastq-quality-control)
    - [Step 3. Mapping](#step-3-mapping)
    - [Setp 4. Differential expressed genes](#step-4-differential-expressed-genes)
- [Issues](#issues)
- [Scripts](#scripts)
    - [Running scripts](#running-scripts)
    - [Verification Scripts](#verification-scripts)
    
    # Summary 

*Note: For a more detailed information about the parameters selection consult the scripts in the :open_file_folder: Scripts*

# Experimental information
- Cell line/ Animal model / Clinical sample
- Specie: Human (*Homo sapiens*)/ Mouse (*Mus musculus*)
- Sample number: 
- Type of read: Paired-end
- Strandness: Read 2 is sense read
- Read length:  bp
- Library average size:  bp
- Library preparation: Illumina TruSeq Stranded mRNA (e.g.)
- Type of sequencer: NovaSeq 600 (e.g.)
- Treatment levels: 

# Contrasts

Contrasts perform in this analysis
- [A] vs [B]


# Workflow

## Step 00. Sample information file

## Step 01. Fastq Quality Control

## Step 02. Trimming

## Step 03. Fastq Quality Control

## Step 04. Mapping

Summary of the mapping results

## Step 05. Gene count file

## Step 06. Differential expressed genes

Results summary, problems found, ...

# Issues

## STAR genome index

In BigData, we can find the corresponding genome index for human species with a read length of 101 bp. This version is not compatible with the STAR version 2.7.10a, so we have to generate a new version. This version is saved in "W:/mponce/AC58/05_DEGs/03_STAR", and will be moved to another folder to make it more accessible.

# Scripts

## Running scripts

- 000_Sample_Info
- 00_MULTIFASTQC
- 01_VAST_TOOLS_Align
- 02_VAST_TOOLS_Combine
- 03_Matt_get_vast

## Verification Scripts

- 00_fastq_file


# References

# Software versions 

FastQC v0.12.1

STAR 2.7.10a

R version 4.2.1 in Rocky

R version 4.3.0 (local)