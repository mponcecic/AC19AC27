# Differentially expressed genes pipeline

:collision: Differential expressed genes analysis workflow from scratch for the AC lab.

# Contents

- []()

- [How this scripts should be run](#how-this-scripts-should-be-run)
- [Verification scripts](#verification-scripts)
- [Folder and scripts](#folder-and-scripts)
- [Additional information](#additional-information)
- [Check corrupt files](#check-corrupt-files)
- [Project_information.md](#project_informationmd)



# How this scripts should be run 

First, you should run the script name 00_Sample_info.R, which generates Sample_info.csv file needed in different analysis steps.
This scripts involves different manual adjustments and you should have access to the library preparation pdf file. In addition, 
this file creates a folder in BigData to save the project analysis which is name as the project.

# Verification scripts

In the script folder, you can find the following scripts: 00_fastq_file.R, 01_fastqc_job_check.R, 02_trimmed_fastq_job_check.R, 02_trimmed_fastq_result_check.R, 03_fastqc_job_check.R, 04_STAR_job_check.R and 04_STAR_result_check.R. As you can see, we can differenciate between three types of files which include *_file*, *_job_check* and *_result_check*. 
- *_file* : There is only one script with this pattern which is 00_fastq_file.R. As a result, you can see the size of each file and should be run while you run 00_Sample_info script. The result will give you valuable information to change parameters such as memory, time, .. while running programs such as FASTQC, cutadapt or STAR in Rocky, among other applications. 
- *_job_check* : After running jobs in Rocky, run the corresponding script to verify each job was **completed**. If one job failed, you just have to consult the .err or the .out file to verify the fail and them adjust the .sh to run it again.
- *_result_check* : This scripts produce the following outputs, one file with how the trimmed fastqs were affected, another with the results of the alignment (how good the mapping to the genome was) and one compeling the trimming results with the alignment results. **MUST RUN AT LEAST 04_STAR_result_check.R TO SEE ALIGNMENT QUALITY**


# Folder and scripts

First, clone this GitHub repository and create a local copy to adjust for your project. 

Second, we created different folders and files available in DATA_shared to run the whole pipeline. In \ref{label1}, you can find the folder with the reference genome, libraries to run R in Rocky and the project fastq files which can be found in other folder but will be usually found in DATA_shared. 

![Folders generated in BigData with all the information concerning this pipeline.\label{label1}](/Schematics/Files_in_DATAshared.png)

After running the first script (`00_Sample_info.csv`), folder with the project name is created in your personal BigData folder, otherwise you must specify the path, a file name `Sample_info.csv` and the folder `utils`. Make sure you have all this folders and files to keep going with the analysis

![Folders generated in BigData personal folder after running the pipeline.\label{label2}](/Schematics/Folder_structure_project.png)

![Complete differentially expressed genes pipeline workflow for AC lab.\label{label3}](/Schematics/Pipeline_flow.png)


![Differentially expressed genes analysis summary.\label{label4}](/Schematics/DEG_flow.png)


# Additional information

## Check corrupt files

First things first, verify that the fastqs you will be using are not corrupt. You can run the following code

```
cd /drives/w/DATA_shared/PROJECT_NAME/FASTQs/

md5sum *gz > check.txt
md5sum -c *gz > check.txt 
cksum *gz > check.txt
```

## Project_information.md

There is a markdown file named `Project_information` in which you can include a summary of every step of the analysis, as well as, the project information. This file can be set as the README in the new GitHub project, it could be used as quick view to all project information, comparisons made, problems, etc.
