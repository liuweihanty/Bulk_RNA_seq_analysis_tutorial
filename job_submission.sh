
#!/bin/bash

#PBS -N run_RNA_seq_wrapper
#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=32gb
#PBS -o /gpfs/data/mcnerney-lab/NGS_analysis_tutorials/RNA_seq/CD34_CUX1_KO/logs/run_RNA_seq_wrapper.out
#PBS -e /gpfs/data/mcnerney-lab/NGS_analysis_tutorials/RNA_seq/CD34_CUX1_KO/logs/run_RNA_seq_wrapper.err

date
module load gcc/6.2.0

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

qsub -v project_folder=$project_dir,fq_F=$i,fq_R=$otherfilename $project_dir/scripts/run_job.sh
      
done
