# Differentially expressed genes pipeline

:collision: Differential expressed genes analysis workflow from scratch for the AC lab.

# Contents

- [Check corrupt files](#check-corrupt-files)
- [How this scripts should be run](#how-this-scripts-should-be-run)
- [Verification scripts](#verification-scripts)
- [Folder and scripts](#folder-and-scripts)

# Check corrupt files

First things first, verify that the fastqs you will be using are not corrupt. You can run the following code

```
cd /drives/w/DATA_shared/PROJECT_NAME/FASTQs/

md5sum *gz > check.txt
md5sum -c *gz > check.txt 
cksum *gz > check.txt
```

> ***Note***: There is a file name README2 in which you can include a summary of every step of the analysis, as well as, the project information. This file can be set as the README in the new Github project, it could be used as quick view to all project information, comparisons made, problems, etc.

# How this scripts should be run 

First, you should run the script name 00_Sample_info.R, which generates Sample_info.csv file needed in different analysis steps.
This scripts involves different manual adjustments and you should have access to the library preparation pdf file. In addition, 
this file creates a folder in BigData to save the project analysis which is name as the project.

# Verification scripts

In the verification script folder, you can find the following scripts
00_fastq_file.R, 01_fastqc_output.R, 02_trimmed_fastq_output.R, 02_trimmed_fastq_result_check.R, 03_fastqc_output.R, 04_STAR_output.R and 04_STAR_result_check.R. As you can see, we can differenciate between three types of files which include *_file*, *_output* and *_result_check*. 
- *_file* : There is only one script with this pattern which is 00_fastq_file.R. As a result, you can see the size of each file and should be run while you run 00_Sample_info script. The result will give you valuable information to change parameters such as memory, time, .. while running programs such as FASTQC, cutadapt or STAR in Rocky, among other applications. 
- *_output* : After running jobs in Rocky, run the corresponding script to verify each job was **completed**. THIS IS A MUST FOR ADVANCING INTO ANOTHER STEP. If one job failed, you just have to consult the .err or the .out file to verify the fail and them adjust the .sh to run it again.
- *_result_check* : This scripts produce the following outputs, one file with how the trimmed fastqs were affected, another with the results of the alignment (how good the mapping to the genome was) and one compeling the trimming results with the alignment results.  


# Folder and scripts


![Folders generated in BigData with all the information concerning this pipeline.](/Schematics/Folder_structure2.png)

![Folders generated in BigData personal folder after running the pipeline](/Schematics/Folder_structure_project2.png)

![Differentially expressed genes analysis summary.](/Schematics/DEG_flow.png)


