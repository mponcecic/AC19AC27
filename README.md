# Differentially expressed genes pipeline

:dna::computer: Differentially expressed genes analysis workflow from scratch for the AC lab.

# Contents

- [Introduction](#introduction)
- [How this project should be run](#how-this-project-should-be-run)
    - [First steps](#star-first-steps)
    - [Fastq quality control](#heavy_check_mark-fastq-quality-control)
    - [Trimming](#scissors-trimming)
    - [Trimmed fastq quality control](#heavy_check_mark-trimmed-fastqs-quality-control)
    - [Mapping reads](#mapping-reads)
    - [Differentially expressed genes analysis](#dnatrophy-differentially-expressed-genes-analysis)
- [Verification scripts](#open_file_folder-verification-scripts)
- [Log files](#open_file_folder-log-files)
- [Additional information](#mag_right-additional-information)
    - [Big Data folders](#open_file_folder-big-data-folders)
    - [Check corrupt files](#biohazard-check-corrupt-files)
    - [Project_information.md](#page_with_curl-project_informationmd)
- [References](#book-references)
- [Software version](#gear-software-version)


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

The scripts are sorted numerically, as well as, the folders in which their outputs are saved. For example, script `01_FASTQC` corresponds to the folder 01_FASTQC. There is only one script that does not follow this rule, `00_Sample_info.R`. More detailed information on each script can be found in it, through this document we will give some guidelines to run the project and understand the structure but reading the scripts can't be skipped.

All the scripts must be modified to adjust the pipeline to your project. The parts of the scripts that must be modified are found between the following characters: 

"# -------------------------------------------"

Furthermore, the following scripts can be directly run in Rocky login node without the need to ask for an interactive session: `01_FASTQC.R`, `01_MULTIFASTQC.R`, `02_TRIMMING.R`, `03_FASTQC.R`, `03_MULTIFASTQC.R` and 

In addition, each step is followed by a verification file. These verification files can check if any job failed but also can generate a summary of the results that gives a general overview of how the step went. The verification files are described in more detail in a [later section](#verification-scripts)

There are only two verification files that must be run one is the `00_fastq_file_check.R` and `04_STAR_result_check.R`. The first file gives information about the size of the files that can be used to estimate data requirements when running jobs in the cluster. The second one shows the results of the alignment, which is a key step for the analysis. If something goes wrong this file should be an alarm for you.  

After running the whole pipeline, you should have a folder in BigData with the folders represented in the following figure. These steps will be discussed in more detail in the [following section](#how-this-project-should-be-run).

![Fig2](/Schematics/Folder_structure_project.png)
*Folders are generated in BigData personal folder after running the pipeline.*

# How this project should be run

In this section, we will discuss the different pipeline steps, previously mentioned. The following figure shows a general view of differentially expressed genes pipeline. 

![Fig3](/Schematics/Pipeline_flow.png)
*Complete differentially expressed genes pipeline workflow for AC lab.*

## :star: First steps

First, clone this GitHub repository and create a local copy for modifying the scripts to run your project, as well as, creating a R project. For cloning the GitHub repository, click on the following [link](https://docs.github.com/es/repositories/creating-and-managing-repositories/cloning-a-repository)

Second, read the Library_Preparation_Report.pdf which can be found in the project folder with the fastq files and will be used in the first script and fill the `00_Sample_info.R`. This script will create a folder in BigData with the project names, cloning the folders utils and generating the `Sample_info.csv` which will be used in several parts of the pipeline. This script involves different manual adjustments, read the script for detailed information.

## :heavy_check_mark: Fastq quality control

The quality control of the raw reads is performed using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses that can be used to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis. This should be done with `01_FASTQC.R` and the `01_MULTIFASTQC.R` generates a summary report with all the information from FastQC using [fastqcr](https://github.com/kassambara/fastqcr). 

Common problems that can appear in the raw counts
- Alterations in base sequence content
- Alterations in the GC content
- Sequence duplications levels
- Overrepresented sequences
- Adapter sequence presence

The overrepresented and adapter sequences can be trimmed off in the following steps.

## :scissors: Trimming

The trimming step is optional, in most cases, we use [cutadapt](https://cutadapt.readthedocs.io/en/stable/), which finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequences from your high-throughput sequencing reads. Indeed, we do not filter adapters only, but low-quality bps and small reads which will increase the percentage of aligned reads.

In some cases, the trimming is not optional as it happens with SMARTer Stranded Total RNA-Seq kit v2-Pico input Mammalian when the insert size is smaller than 150 bps. This trimming is described in the script, `02_TRIMMING.R`.

## :heavy_check_mark: Trimmed fastqs quality control

The quality control of the trimmed reads must be performed to check the changes. This step is exactly equal to the [previous section](#fastq-quality-control) but running the scripts `03_FASTQC.R` and `03_MULTIFASTQC.R`.

## :dna: Mapping reads

The read alignment is performed using Spliced Transcripts Alignment to a Reference ([STAR](https://github.com/alexdobin/STAR)). STAR is a fast RNA-seq read mapper, with support for splice-junction and fusion read detection. STAR aligns reads by finding the maximal mappable prefix hits between reads (or read pairs) and the genome, using a Suffix Array index (also known as genome indexes). The script performing this step is `04_STAR_align.R`. 

In addition, the script `04_STAR_genome_index.R` aims to generate genome indexes that are not available which human and mouse genomes with a read length of 51, 101 and 151. If you have another read length, you should generate the genome indexes before running the alignment of the reads. 

The alignment result must be transformed into a gene count matrix with the rows being the Ensembl identifier of the gene and the columns, the sample names. To accomplish this be careful with the data strandedness defined by the library preparation method. The **strandedness** is discussed in more detail in the script, `04_STAR_GeneCounts.R`.

## :dna::trophy: Differentially expressed genes analysis 

> **Note**: *Consult the [Differential gene expression workshop](https://github.com/hbctraining/DGE_workshop) for an overview of the DEGs analysis*  


The differentially expressed genes analysis presents different steps, (1) quality control of the samples, (2) applying DESeq2 and/or other methods per contrast for DEGs and (3) comparing the results among different methods. 

![Fig4](/Schematics/DEG_flow.png)
*Differentially expressed genes analysis summary*

The first step is to perform quality control (`05_DEGs_v1_qc.R`) among all the samples without considering the treatment. As a result, plots and data tables. The plots, saved in :open_file_folder: QC, will highlight if there are outliers or another source variance in the data that is not explained by the experimental condition. 

In a folder named :open_file_folder: Results, the data tables are saved and all relevant information in the following analysis will be saved here. The results files are metadata, filtered data and transformed data (VST/RLOG, CPM) for raw and filtered data. The filtered data is based on `filterByExprs()` in the `EdgeR` package which is applied to all the data. This step is key for EdgeR and limma-voom.

The second step is to apply one or different methods, which are DESeq2 (`05_DEGs_v2_DESeq2.R`), EdgeR (`05_DEGs_v2_EdgeR.R`), limma-voom (`05_DEGs_v2_limma-voom.R`) and Wilcoxon test (`05_DEGs_v2_Wilcoxon.R`), to obtain the DEGs. The method used for the analysis can be chosen due to several reasons, here you can find a summary table. 

| Method | Data | Group of interest |
|:----------:|:---------:|:------------:| 
| DESeq2 | Cell lines/Model organism | Researchers |
| EdgeR | Cell lines/Model organism | Computational team |
| limma-voom | Cell lines/Model organism | Computational team |
| Wilcoxon | Clinical data | Computational team |

> **Note**: DESeq2 method is applied on raw counts (DESeq2_NoFilter) and filtered counts (DESeq2). **The results given to the researchers are the obtain based on the raw countsafter applying DESeq2**.

This step is specific per each comparison, so a subset of the matrix is made for the analysis. 
As a result, a folder with the name of the comparison is generated and inside this folder, you can find different folders with the name of the method used and Results folder.
In the folder with the method name, all the figures for that method are saved. In the Results folder, the contrast result for this comparison and method is saved.

However, all the contrast results per method are saved and compiled in three different files saved in the Result folder outside the comparison. These files are: 

- `Method_Project_name;All_VSTblindFALSE_padj_fdrcutoff_lfc_lfccutoff.xlsx` **This is the file given to the researchers**. It includes each comparison in separate sheets. Take a brief look at how the table should look like.

![Fig6](/Schematics/output_xlsx.png)

- `Method_Project_name;All_VSTblindFALSE_padj_fdrcutoff_lfc_lfccutoff.txt` This file contains the same information as before but with all the comparisons in it. This file is key for following analyses such as trajectories. 

- `Summary_tab_Method_Project_name_padj_fdrcutoff_lfc_lfccutoff.csv` Summary of the design formula, DEGs, upregulated and downregulated genes of all comparisons. Offers a quick view of the numeric results.

The third step compares the results among the different methods (`05_DEGs_v3_Comparison.R`). The figures are saved in a folder named Comparison and the data is saved in both Results folders. This step is optional. 

As a result, of running this part of the pipeline you will obtain the following folders. In purple are marked the folders in which relevant figures to the researchers are found. In addition, the only xlsx file found in Results must be given to. 

![Fig7](/Schematics/05_DEG_folder.png)

# :open_file_folder: Verification scripts

In the script folder, you can find the following scripts 
`00_fastq_file.R`, `01_fastqc_job_check.R`, `02_trimmed_fastq_job_check.R`, `02_trimmed_fastq_result_check.R`, `03_fastqc_job_check.R`, `04_STAR_job_check.R` and `04_STAR_result_check.R`. 

As you can see, we can differentiate between three types of files which include *_file*, *_job_check* and *_result_check*. 
- *_file* : There is only one script with this pattern which is 00_fastq_file.R. As a result, you can see the size of each file and should be run while you run 00_Sample_info script. The result will give you valuable information to change parameters such as memory and time while running programs such as FASTQC, cutadapt or STAR in Rocky, among other applications. 
- *_job_check* : After running jobs in Rocky, run the corresponding script to verify each job was **completed**. If one job fails, you just have to consult the .err or the .out file to verify the error and then adjust the .sh to run it again.
- *_result_check* : These scripts produce the following outputs, one file with how the trimmed fastqs were affected, another with the results of the alignment (how good the mapping to the genome was) and one compeling the trimming results with the alignment results. 
**04_STAR_result_check.R MUST BE RUN TO FOLLOW THE ANALYSIS**

# :open_file_folder: Log files

The log files are another result of the scripts and are saved in the log folder. These log files are files used to save paths, variables, contrasts, colors, etc which will be used in other steps of the analysis. It's used to save all the relevant information of each step, to avoid errors due to a variable appearing in more than one script and as a way to increase reproducibility. 

The log files generated if you run all the scripts are `0_Sample_info_DATE.log`, `4_STAR_DATE.log`, `4_STAR_GeneCounts_DATE.log`, `5_DEG_qc_DATE.log`, `5_DEG_v2_DESeq2_DATE.log`, `5_DEG_v2_DESeq2_NoFilter_DATE.log`, `5_DEG_v2_EdgeR_DATE.log`, `5_DEG_v2_limma-voom_DATE.log` and `5_DEG_v3_Methodscompare_DATE.log`. 

# :mag_right: Additional information

## :open_file_folder: Big Data folders

We have created different folders and files available in DATA_shared to run the whole pipeline. They can be seen in the following figure. 

![Fig1](/Schematics/Files_in_DATAshared.png)
*Folders generated in BigData with all the information concerning this pipeline.*

The information in the following folders are: 
- Genome_Rocky: The reference genomes for Human and Mouse used to align the reads.
- R_Rocky: R libraries to run the scripts in Rocky
- Fastq files corresponding to the project, which can be found in other folders but will be usually found in DATA_shared. 

## :biohazard: Check corrupt files

First things first, verify that the fastqs you will be using are not corrupt. You can run the following different codes to see if the files are corrupted. 

```
cd /drives/w/DATA_shared/PROJECT_NAME/FASTQs/

md5sum *gz > check.txt
md5sum -c *gz > check.txt 
cksum *gz > check.txt
```

## :page_with_curl: Project_information.md

In the folder Schematics, you can find a markdown file named `Project_information` in which you can include a summary of every step of the analysis, as well as, the project information. This file can be set as the README in the new GitHub project, it could be used as a quick view of all project information, comparisons made, problems, etc.

# :book: References

* Harvard Chan Bioinformatics Core (HBC). 2022. [Differential gene expression workshop](https://github.com/hbctraining/DGE_workshop)
 
 * Love, M.I., Anders, S. and Huber, W. 2023. [Analyzing RNA-seq data with DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#log-fold-change-shrinkage-for-visualization-and-ranking)


# :gear: Software version

FastQC v0.12.0

cutadapt version 4.3

STAR version 2.7.10a

R version 4.2.1 (Rocky)