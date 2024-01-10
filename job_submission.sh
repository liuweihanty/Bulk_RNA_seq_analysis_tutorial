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
