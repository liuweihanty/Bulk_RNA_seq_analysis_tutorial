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
  * create your project folder. **/CD34_CUX1_CnR/**
  * create four sub-folders underneath your project folder
     * **/CD34_CUX1_CnR/input** $~~~$ trimmed fastqs
     * **/CD34_CUX1_CnR/output** $~~~$ the analysis output
     * **/CD34_CUX1_CnR/logs** $~~~$ the error and output records files for debugging
     * **/CD34_CUX1_CnR/scripts** $~~~$ the analysis scripts <br>
     
  Your working directory should look like this by now: <br>
     <img src="https://github.com/liuweihanty/ChIP_analysis_tutorial/blob/7064d2d19c974fa2adacc081f888f956b49ce070/figures/working_directory_before_run.png" alt="repo_demo" width="350" height="200">

           
* ### Set up the job running scripts
     Now that we have the adaptor trimmed fastqs, it's time to proceeed to next steps. In the flow chart above, we finished steps 1 and 2 so far. Step 3 to 6 will be implemented in an automated workflow, which is organized into two bash scripts: <br>
    * **job_submission.sh**: this script specify and select the two fastq files(forward and reverse reads) for each sample, and send the fastqs to the "run_job.sh" script below
    * **run_job.sh**:  this scripts takes in the forward and reverse fastqs for each sample from the job_submission.sh file above, and performs steps 3-6 in the flow chart on each sample (except IDR analysis), in a paralelled fashion( all samples will be simutaneously analyzed), so no matter how many samples you have, you just need to run once. <br>

   
    Now let's take a look inside of an example of each file and I will explain what each code chunk do: <br>
    
    **job_submission.sh** For each sample, this script below find the forward read(R1) fastq file, and subsequently locate the reverse read(R2) file for that same sample(it can do so because the fastq file names you got from the sequencing core differ only in "R1" and "R2" part for the file name). This script essently locate the forward and reverse reads fastq files parallelly for each sample, and feed them into the "run_jobs.sh" file to run all the analysis steps. **What you need to do**: change all the directory path to your project directory.
    ```bash
    
     #!/bin/bash
     
     #PBS -N CD34_CUX1_CnR
     #PBS -S /bin/bash
     #PBS -l walltime=24:00:00
     #PBS -l nodes=1:ppn=8
     #PBS -l mem=32gb
     #PBS -o /gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/logs/run_CnR_wrapper.out
     #PBS -e /gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/logs/run_CnR_wrapper.err
     
     date
     module load gcc/6.2.0
     
     
     #specify your project directory
     project_dir=/gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR
     
     #change directory to where the fastq files are
     cd $project_dir/input
     
     #this for loop will take the input fastq files and run the scripts for all of them one pair after another
     
     for i in $(ls *R1*.gz)
     do
     otherfilename="${i/R1/R2}"
     echo $i
     echo $otherfilename
     
     
     #here you need to specify whether to perform macs2 peak calling by include the -macs2 flag or not. If you include, you need to specify either -p or -q significance threshold followed by a number. Do not specify both p and q values
     
     qsub -v project_path=$project_dir,fq_F=$i,fq_R=$otherfilename,-macs2,-p=0.1 /gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/scripts/run_job.sh 
           
     done

   
    ```
    * **Important**: Within the job_submission.sh file in the last line qsub,  you have the choice of specifying whether to run macs2 peak calling step. There are three mandatory flags in the qsub command:
        * -macs2: whether to run macs2, if you include this flag, the problem will run macs2 peak caller, if not, the program will skip macs2.
        * -p: p value for macs2. Use p value as the significant thrshold and specify it to be 0.1 if you are running IDR.
        * -q: q value for macs2. If you are not running IDR (eg you just have one rep), use this adjusted p value for macs2 instead. Only specify either p or q, not both!
          
    **run_job.sh** This is the script that performs the actual analysis for each sample. The input are fastq files, and it will output:<br>
    *the aligned and filtered bam file <br>
    *bigwig files for each individual replicate <br>
    *individual bigwig files for each replicate, each sample <br>
    *MACS2 called peaks in narrowpeak format <br>

    **You don't need to modify anything for this script** <br>

    **Note:** <br>
    This tutorial uses p=0.1 for MACS2 peak calling, this is because the downstream [IDR](https://github.com/nboley/idr) workflow requires a loose significance threshold. IDR find consensus peaks across two biological replicates. It's best to use IDR if you have replicates. If you just have one rep (eg for a pilot study), since you are not using IDR in this case, you can set q=0.1 etc for MACS2 for an actual robust significance threshold.<br>
    
  
* ### Run the job
    *change directory to the scripts folder that contains your job_submission.sh and run_job.sh file, type in ``` chmod +x * ```, this give execution rights to your scripts <br>
    * type in ```./job_submission.sh``` and the job should start to run


* ### IDR analysis (If you chose to run macs2)
    * After the program finishes run, go to the macs2 result folder <project_folder/output/macs2>, and download the narrowPeak files to your computer. These files are extended bed files that contain each called peak as a row, with additional columns(etc. significance, binding intensity). At this stage, you can change the narrowPeak file name to something that makes sense to you.
    * If you don't have IDR installed, install it on your local computer [Instructions here](https://github.com/nboley/idr)
    * Now let's run IDR
      ```
      ##Sort your narrowPeak files by the -log10(p-value) column
      sort -k8,8nr CD34_CUX1_CnR_rep1.narrowPeak > CD34_CUX1_CnR_rep1.sorted.narrowPeak
      sort -k8,8nr CD34_CUX1_CnR_rep2.narrowPeak > CD34_CUX1_CnR_rep2.sorted.narrowPeak
      ##run idr
      idr --samples CD34_CUX1_CnR_rep1.sorted.narrowPeak CD34_CUX1_CnR_rep2.sorted.narrowPeak \
          --input-file-type narrowPeak \
          --rank p.value \
          --output-file CD34_CUX1_CnR_idr \
          --plot \
          --log-output-file CD34_CUX1_CnR.idr.log
     * The output file CD34_CUX1_CnR.idr contains the consensus peaks across the two replicates. For a detailed explanantion of what columns are inside, please see [HBC training](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html). What we care about here is column 5, which contains the scaled IDR value = -125*log2(IDR) For example, peaks with an IDR of 0.1 have a score of 415, peaks with an IDR of 0.05 have a score of int(-125log2(0.05)) = 540.
     * Use this command to filter out peaks with IDR < 0.05 and compile to a bed file. This is the final file that contains your consensus peaks.
       ``` awk '{if($5 >= 540) print $0}' CD34_CUX1_CnR_idr > CD34_CUX1_CnR_idr_005.bed ``` 




## Other situations
   * ### Single end sequencing analysis
   If you are running single end sequencing (eg.ChIP-seq), please modidy the job submission code as follows:
   
    ```
     #!/bin/bash
     
     #PBS -N CD34_CUX1_CnR
     #PBS -S /bin/bash
     #PBS -l walltime=24:00:00
     #PBS -l nodes=1:ppn=8
     #PBS -l mem=32gb
     #PBS -o /gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/logs/run_CnR_wrapper.out
     #PBS -e /gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR/logs/run_CnR_wrapper.err
     
     date
     module load gcc/6.2.0
     
     
     #specify your project directory
     project_dir=/gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ChIP_seq/CD34_CUX1_CnR
     
     #change directory to where the fastq files are
     cd $project_dir/input
     
     for i in $(ls *.gz)
     do
     echo $i
  
     qsub -v project_path=$project_dir,fq=$i,-macs2,-p 0.1 /gpfs/data/mcnerney-lab/.../CD34_CUX1_CnR/logs/scripts/run_job.sh
           
     done
    
    ```

   And modify the following places in the run_job.sh script <br>
    1.At the start, instead of ``` echo $fq_F and $fq_R ```, just change to ``` echo $fq ``` <br>
    2.Do the same thing in the bwa mem code <br>
    3.Change the ``` -f BAMPE ``` flag in macs2 to ``` -f BAM ``` <br>
