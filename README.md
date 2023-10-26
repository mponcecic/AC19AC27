# DEG_Reference

:collision: Differential expressed genes analysis workflow from scratch for the AC lab.

# Prior step

First things first, verify that the fastqs you will be using are not corrupt. You can run the following code

```
cd /drives/w/DATA_shared/PROJECT_NAME/FASTQs/

md5sum *gz
md5sum -c *gz 
cksum *gz
```

In addition, there is a file name README2 in which you can include a summary of every step of the analysis. This 
file can be set as the README in the new Github project you will create and will show in a quick view all the 
information from the project.

# How this scripts should be run 

First, you should run the script name 00_Sample_info.R, which generates Sample_info.csv file needed in different analysis steps.
This scripts involves different manual adjustments and you should have access to the library preparation pdf file. In addition, 
this file creates a folder in BigData to save the project analysis which is name as the project.

# Work done 

Sample information script 

# Folder
- Analysis scripts
- Verification scripts
