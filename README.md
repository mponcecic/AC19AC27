# Differentially expressed genes pipeline

:collision: Differentially expressed genes analysis workflow from scratch for the AC lab.

# Contents

- [Introduction](#introduction)
- [How this project should be run](#how-this-project-should-be-run)

    0. [First steps](#first-steps)
    1. [Fastq quality control](#fastq-quality-control)
    2. [Trimming](#trimming)
    3. [Trimmed fastq quality control](#trimmed-fastq-quality-control)
    4. [Mapping reads](#mapping-reads)
    5. [Differentially expressed genes analysis](#differentially-expressed-genes-analysis)
- [Verification scripts](#verification-scripts)
- [Additional information](#additional-information)
    - [Big Data folders](#big-data-folders)
    - [Check corrupt files](#check-corrupt-files)
    - [Project_information.md](#project_informationmd)


# Introduction

This repository contains the differential expressed genes analysis pipeline for the AC lab. This file gives general guidelines to run the pipeline but more detailed information can be found in the scripts. Before going into detail in the analysis workflow, we will comment on the folders found in the repository which are Schematics, Scripts and utils.

- Scripts Folder with all the scripts of each step, verification scripts and help scripts for different programs used. 
- utils Files to load libraries, functions created for the pipeline, how to generate an interactive session and how to generate a cutadapt environment in Rocky.
- Schematics Folder with the figures shown in the README.

The pipeline presents the following key steps, which can be summarised in the following list:

0. Run `00_Sample_info.R`
1. Fastq quality control (FastQC)
2. Trimming the reads (*optional*)
3. Fastq quality control of the trimmed reads (FastQC)
4. Mapping the reads and extracting the gene count matrix
5. Differentially expressed genes analysis using DESeq2, EdgeR, limma-voom and/or Wilcoxon test

![Fig3](/Schematics/Pipeline_flow.png)
*Complete differentially expressed genes pipeline workflow for AC lab.*

The scripts are sorted numerically, as well as, the folders in which their outputs are saved. For example, script `01_FASTQC` corresponds to the folder 01_FASTQC. There is only one script that does not follow this rule, `00_Sample_info.R`. More detailed information of each script can be found in it, through this document we will give some guidelines to run the project and understand the structure but reading the scripts can't be skipped.

All the scripts must be modified to adjust the pipeline to your project. The parts of the scripts that must be modified are found between the following characters: 
"# -------------------------------------------"

Furthermore, the following scripts can be directly run in Rocky login node without the need to ask for an interactive session: `01_FASTQC.R`, `01_MULTIFASTQC.R`, `02_TRIMMING.R`, `03_FASTQC.R`, `03_MULTIFASTQC.R` and 

In addition, each step is followed by a verification file. These verification files can check if any job failed but also can generate a summary of the results that gives a general overview of how the step went. The verification files are described in more detail in a [later section](#verification-scripts)

There are only two verification files that must be run one is the `00_fastq_file_check.R` and `04_STAR_result_check.R`. The first file gives information about the size of the files that can be used to estimate data requirements when running jobs in the cluster. The second one shows the results of the alignment, which is a key step for the analysis. If something goes wrong this file should be an alarm for you. 

After running the whole pipeline, you should have a folder in BigData with the folders represented in the following figure. These steps will be discussed in more detail in the [following section](#how-this-project-should-be-run).

![Fig2](/Schematics/Folder_structure_project.png)
*Folders are generated in BigData personal folder after running the pipeline.*

# How this project should be run

## First steps

First, clone this GitHub repository and create a local copy for modifying the scripts to run your project, as well as, creating a R project. For cloning the GitHub repository, click on the following [link](https://docs.github.com/es/repositories/creating-and-managing-repositories/cloning-a-repository)

Second, read the Library_Preparation_Report.pdf which can be found in the project folder with the fastq files and will be used in the first script and fill the `00_Sample_info.R`. This script will create a folder in BigData with the project names, cloning the folders utils and generating the `Sample_info.csv` which will be used in several parts of the pipeline. This script involves different manual adjustments, read the script for detailed information.

## Fastq quality control

The quality control of the raw reads is performed using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses that can be used to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis. This should be done with `01_FASTQC.R` and the `01_MULTIFASTQC.R` generates a summary report with all the information from FastQC using [fastqcr](https://github.com/kassambara/fastqcr). 

Common problems that can appear in the raw counts
- Alterations in base sequence content
- Alterations in the GC content
- Sequence duplications levels
- Overrepresented sequences
- Adapter sequence presence

The overrepresented and adapter sequences can be trimmed off in the following steps.

## Trimming

The trimming step is optional, in most cases, we use [cutadapt](https://cutadapt.readthedocs.io/en/stable/), which finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequences from your high-throughput sequencing reads. Indeed, we do not filter adapters only, but low-quality bps and small reads which will increase the percentage of aligned reads.

In some cases, the trimming is not optional as it happens with SMARTer Stranded Total RNA-Seq kit v2-Pico input Mammalian when the insert size is smaller than 150 bps. This trimming is described in the script, `02_TRIMMING.R`.

## Trimmed fastq quality control

The quality control of the trimmed reads must be performed to check the changes. This step is exactly equal to the [previous section](#fastq-quality-control) but running the scripts `03_FASTQC.R` and `03_MULTIFASTQC.R`.

## Mapping reads

The read alignment is performed using Spliced Transcripts Alignment to a Reference ([STAR](https://github.com/alexdobin/STAR)). STAR is a fast RNA-seq read mapper, with support for splice-junction and fusion read detection. STAR aligns reads by finding the maximal mappable prefix hits between reads (or read pairs) and the genome, using a Suffix Array index (also known as genome indexes). The script performing this step is `04_STAR_align.R`. 

In addition, the script `04_STAR_genome_index.R` aims to generate genome indexes that are not available which human and mouse genomes with a read length of 51, 101 and 151. If you have another read length, you should generate the genome indexes before running the alignment of the reads. 

The alignment result must be transformed into a gene count matrix with the rows being the Ensembl identifier of the gene and the columns, the sample names. To accomplish this be careful with the data strandedness defined by the library preparation method. The **strandedness** is discussed in more detail in the script, `04_STAR_GeneCounts.R`.

## Differentially expressed genes analysis 

> **Note**: *Consult the [Differential gene expression workshop](https://github.com/hbctraining/DGE_workshop) for an overview of the DEGs analysis*  

`05_DEGs_v1_qc.R`
`05_DEGs_v2_DESeq2.R`
`05_DEGs_v2_EdgeR.R`
`05_DEGs_v2_limma-voom.R`
`05_DEGs_v2_Wilcoxon.R`
`05_DEGs_v3_Comparison.R`






you can find a workshop given by Harvard 

![Fig4](/Schematics/DEG_flow.png)
*Differentially expressed genes analysis summary*

# Verification scripts

In the script folder, you can find the following scripts: `00_fastq_file.R`, `01_fastqc_job_check.R`, `02_trimmed_fastq_job_check.R`, `02_trimmed_fastq_result_check.R`, `03_fastqc_job_check.R`, `04_STAR_job_check.R` and `04_STAR_result_check.R`. As you can see, we can differentiate between three types of files which include *_file*, *_job_check* and *_result_check*. 
- *_file* : There is only one script with this pattern which is 00_fastq_file.R. As a result, you can see the size of each file and should be run while you run 00_Sample_info script. The result will give you valuable information to change parameters such as memory and time while running programs such as FASTQC, cutadapt or STAR in Rocky, among other applications. 
- *_job_check* : After running jobs in Rocky, run the corresponding script to verify each job was **completed**. If one job fails, you just have to consult the .err or the .out file to verify the error and then adjust the .sh to run it again.
- *_result_check* : These scripts produce the following outputs, one file with how the trimmed fastqs were affected, another with the results of the alignment (how good the mapping to the genome was) and one compeling the trimming results with the alignment results. 
**04_STAR_result_check.R MUST BE RUN TO FOLLOW THE ANALYSIS**


# Additional information

## Big Data folders

We have created different folders and files available in DATA_shared to run the whole pipeline. They can be seen in the following figure. 

![Fig1](/Schematics/Files_in_DATAshared.png)
*Folders generated in BigData with all the information concerning this pipeline.*

The information in the following folders are: 
- Genome_Rocky: The reference genomes for Human and Mouse used to align the reads.
- R_Rocky: R libraries to run the scripts in Rocky
- Fastq files corresponding to the project, which can be found in other folders but will be usually found in DATA_shared. 

## Check corrupt files

First things first, verify that the fastqs you will be using are not corrupt. You can run the following code

```
cd /drives/w/DATA_shared/PROJECT_NAME/FASTQs/

md5sum *gz > check.txt
md5sum -c *gz > check.txt 
cksum *gz > check.txt
```

## Project_information.md

In the folder Schematics, you can find a markdown file named `Project_information` in which you can include a summary of every step of the analysis, as well as, the project information. This file can be set as the README in the new GitHub project, it could be used as a quick view of all project information, comparisons made, problems, etc.

# References

* Harvard Chan Bioinformatics Core (HBC). 2022. [Differential gene expression workshop](https://github.com/hbctraining/DGE_workshop)
 
 * Love, M.I., Anders, S. and Huber, W. 2023. [Analyzing RNA-seq data with DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#log-fold-change-shrinkage-for-visualization-and-ranking)


