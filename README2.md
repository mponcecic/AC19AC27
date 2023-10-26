# PROJECT NAME

:collision: Brief description

# Contents

- [Summary](#summary)
- [Experimental Information](#experimental-information)
- [Contrasts](#contrasts)
- [Workflow](#)
    - [Step 0. Fastq Quality Control](#step-0-fastq-quality-control)
    - [Step 1. Mapping](#step-1-mapping)
    - [Step 2. Hypothesis testing](#setp-2-hypothesis-testing-with-matt)
    - [Step 3. Filtering differentially spliced events](#step-3-selecting-differentially-spliced-events)
    - [Step 4. Differentially splicing analysis](#step-4-differentially-splicing-analysis)
- [Issues](#issues)
- [Scripts](#scripts)
    - [Running scripts](#running-scripts)
    - [Verification Scripts](#verification-scripts)
    
    # Summary 

*Note: For a more detailed information about the parameters selection consult the scripts in the :open_file_folder: Scripts*

# Experimental information
- Cell line DUI45
- Specie: Human (*Homo sapiens*)
- Sample number: 16
- Type of read: Paired-end
- Strandness: Read 2 is sense read
- Read length: 101 bp
- Library average size: 327 bp
- Library preparation: Illumina TruSeq Stranded total RNA with rRNA depleted
- Type of sequencer: NovaSeq 600
- Time course: 4, 24, 48 hours

# Contrasts

All posible contrast must be performed in this data set. A list of the different contrasts is shown: 
- N vs 4
- N vs 24
- N vs 48
- 4 vs 24
- 4 vs 48
- 24 vs 48

# Workflow

## Step 0. Fastq Quality Control


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

- Vast-tools 

Tapial, J., Ha, K.C.H., Sterne-Weiler, T., Gohr, A., Braunschweig, U., Hermoso-Pulido, A., Quesnel-Vallières, M., Permanyer, J., Sodaei, R., Marquez, Y., Cozzuto, L., Wang, X., Gómez-Velázquez, M., Rayón, M., Manzanares, M., Ponomarenko, J., Blencowe, B.J., Irimia, M. (2017). An atlas of alternative splicing profiles and functional associations reveals new regulatory programs and genes that simultaneously express multiple major isoforms. Genome Res, 27(10):1759-1768

- Matt

Gohr, A., Irimia, M. (2018). Matt: Unix tools for alternative splicing analysis. Bioinformatics, 35:130-132.

- Hypoxia paper which inspired this study

Han, J., Li, J., Ho, J. C., Chia, G. S., Kato, H., Jha, S., ... & Lee, K. L. (2017). Hypoxia is a key driver of alternative splicing in human breast cancer cells. Scientific reports, 7(1), 4108.

- Time course paper

Vivori, C., Papasaikas, P., Stadhouders, R., Di Stefano, B., Rubio, A. R., Balaguer, C. B., ... & Valcárcel, J. (2021). Dynamics of alternative splicing during somatic cell reprogramming reveals functions for RNA-binding proteins CPSF3, hnRNP UL1, and TIA1. Genome Biology, 22(1), 1-30.

- Inspiration of the previous paper

Cieply, B., Park, J. W., Nakauka-Ddamba, A., Bebee, T. W., Guo, Y., Shang, X., ... & Carstens, R. P. (2016). Multiphasic and dynamic changes in alternative splicing during induction of pluripotency are coordinated by numerous RNA-binding proteins. Cell reports, 15(2), 247-255.

# Software versions 

VAST-TOOLS v2.5.1

Matt 1.3.1

R version 4.3.0 