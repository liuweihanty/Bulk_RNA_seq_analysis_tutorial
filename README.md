# Table of contents <br>
 - [Introduction](#introduction)
 - [Demo data](#demo_data)
 - [Analysis workflow](#analysis_workflow)
 - [Step by step analysis](#Step_by_step_analysis)
 - [Other situations](#Other_situations) 
    - [Single end sequencing analysis](#single_end_sequencing_analysis)

  
## Introduction <br>
This tutorial walks step-by-step tutorial of analysis pipeline for bulk RNA seq. 

## Demo data
I will be using an example data set to illustrate this workflow. This is a paired-end bulk RNA seq human CD34+ HSPC, investigating the effect of CUX1 KO on human HSPC. There are two conditions, control gHPRT and knockout gCUX1. There are two replicates for each condition, each replicate is done in different days

gHPRT <br>
rep1(day1): <br>
MMc-RS5-SK-01_S131_R1_001.fastq.gz <br>
MMc-RS5-SK-01_S131_R2_001.fastq.gz <br>

rep2(day2): <br>
MMc-RS5-SK-02_S132_R1_001.fastq.gz <br>
MMc-RS5-SK-02_S132_R2_001.fastq.gz <br>

gCUX1 <br>
rep1(day1): <br>
MMc-RS5-SK-03_S133_R1_001.fastq.gz <br>
MMc-RS5-SK-03_S133_R2_001.fastq.gz <br>

rep2(day2): <br>
MMc-RS5-SK-05_S135_R1_001.fastq.gz <br>
MMc-RS5-SK-05_S135_R2_001.fastq.gz <br>


## Analysis workflow
![GitHub Logo](https://github.com/liuweihanty/Bulk_RNA_seq_analysis_tutorial/blob/main/figures/bulk_RNA_seq_analysis_workflow.png)

The steps in grey are run on the cluster, in blue are run locally on your computer

## Step by step analysis
* ### Run fastqc 
  You can do this in either linux or R. [fastqcr](http://www.sthda.com/english/wiki/fastqcr-an-r-package-facilitating-quality-controls-of-sequencing-data-for-large-numbers-of-samples) package provided an easy implementation in R language. You can run this in your local desktop. The fastqcr code is attached in this folder named fastqc.Rmd. Please see the [document](https://github.com/liuweihanty/ChIP_analysis_tutorial/blob/f4982c5fd9c9e25d493fb50f1813dc429562869b/fastqc.Rmd) for details.

* ### Remove illumina library adaptor.
  * You will obtain adaptor sequence from the sequencing facility. If not, from fastqc results, in the last section "adaptor sequence" you will see it. Typical Illumina library sequencing adaptor sequences could be seen [here](https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314) <br>
  * Use [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) to trim the adaptors. See [here](https://cutadapt.readthedocs.io/en/stable/installation.html) for how to install Cutadapt <br>
  * Run the code below, swap the AACCGGTT for your actual adaptor sequence <br>
   ```cutadapt -a AACCGGTT -o output.fastq input.fastq```

* ### (Optional) Run fastqcr again to ensure the adaptors are successfully removed
  
* ### Set up your working directory. (the demo example folder names are written in parenthesis)
  * create your project folder. **/CD34_CUX1_KO/**
  * create four sub-folders underneath your project folder
     * **/CD34_CUX1_KO/input** $~~~$ trimmed fastqs
     * **/CD34_CUX1_KO/output** $~~~$ the analysis output
     * **/CD34_CUX1_KO/logs** $~~~$ the error and output records files for debugging
     * **/CD34_CUX1_KO/scripts** $~~~$ the analysis scripts <br>
     
  Your working directory should look like this by now: <br>
     <img src="https://github.com/liuweihanty/Bulk_RNA_seq_analysis_tutorial/blob/main/figures/working_directory_demo.png" alt="repo_demo" width="350" height="200">

           
* ### Set up the job running scripts
     Now that we have the adaptor trimmed fastqs, it's time to proceeed to next steps. Genome alignment using STAR, samtools filtering and featureCounts are automated using the two scripts below <br>
    * **job_submission.sh**: this script specify and select the two fastq files(forward and reverse reads) for each sample, and send the fastqs to the "run_job.sh" script below
    * **run_job.sh**:  this scripts takes in the forward and reverse fastqs for each sample from the job_submission.sh file above, and performs steps stated above(STAR, samtools, featurecounts) on each sample (except IDR analysis), in a paralelled fashion( all samples will be simutaneously analyzed), so no matter how many samples you have, you just need to run once. <br>

    Now let's take a look inside of an example of each file and I will explain what each code chunk do: <br>
    
    **job_submission.sh** For each sample, this script below find the forward read(R1) fastq file, and subsequently locate the reverse read(R2) file for that same sample(it can do so because the fastq file names you got from the sequencing core differ only in "R1" and "R2" part for the file name). This script essently locate the forward and reverse reads fastq files parallelly for each sample, and feed them into the "run_jobs.sh" file to run all the analysis steps. **What you need to do**: change the directory path "project_dir" to your project directory.
    ```bash
   
    #!/bin/bash
    #specify your project directory
    project_dir=/gpfs/data/mcnerney-lab/NGS_analysis_tutorials/RNA_seq/CD34_CUX1_KO
    
    #change directory to where the fastq files are
    cd $project_dir/input
    
    #this for loop will take the input fastq files and run the scripts for all of them one pair after another
      
    for i in $(ls *R1*.gz)
    do
    otherfilename="${i/R1/R2}"
    echo $i
    echo $otherfilename
    
    sbatch --job-name=run_RNA_seq_wrapper --time=12:00:00 \
       -o $project_dir/logs/run_RNA_seq.out \
       -e $project_dir/logs/run_RNA_seq.err \
       --partition=tier2q \
       --nodes=1 \
       --ntasks-per-node=8 \
       --mem-per-cpu=10000 \
       --wrap="sh $project_dir/scripts/run_job.sh $i $otherfilename $project_dir"
          
    done
   
    ```
          
    **run_job.sh** This is the script that performs the actual analysis for each sample. The input are fastq files, and it will output:<br>
    *the aligned and q30 filtered bam file <br>
    *individual bigwig files for each replicate, each sample <br>
    *Raw read count matrix (to be fed into DESeq2) <br>

    **You don't need to modify anything for this script** <br>

    **Important:** Do no remove duplicated reads. Dups are expected for RNA seq.
  
* ### Run the job
    *change directory to the scripts folder that contains your job_submission.sh and run_job.sh file, type in ``` chmod +x * ```, this give execution rights to your scripts <br>
    * type in ```./job_submission.sh``` and the job should start to run


* ### Differential expression analysis 
  Please see the DESeq.Rmd scripts in this folder for details. Run this step locally on your computer using RStudio. The input is the featurecounts.txt result raw count matrix from last step. This current step will generate the DEG result table and graphs(eg PCA,heatmap, volcano plot, MA plot etc.)


