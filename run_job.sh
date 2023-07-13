#!/bin/bash
#This scrip takes the forward and reverse read fastq files of RNA seq, runs FastQC, STAR alignment, samtool filtering and featureCounts
#the final output is a count matrix
#USAGE: sh run_RNA_seq.sh <path to forward read fastq> <path to reverse read fastq> 

#PBS -N run_RNA_seq
#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=32gb
#PBS -o /gpfs/data/mcnerney-lab/NGS_analysis_tutorials/RNA_seq/CD34_CUX1_KO/logs/run_RNA_seq.out
#PBS -e /gpfs/data/mcnerney-lab/NGS_analysis_tutorials/RNA_seq/CD34_CUX1_KO/logs/run_RNA_seq.err


#specify the input fastq files
echo $fq_F
echo $fq_R
echo $project_folder

#specify the location of reference genome directory
GENOME_DIR=/gpfs/data/mcnerney-lab/liuweihan/CUX1_CASP_diff_ref_transcriptome_bulk/hg19_refgenome_CUX1_CASP_diff

#specify the input GTF file
GTF=/gpfs/data/mcnerney-lab/liuweihan/CUX1_CASP_diff_ref_transcriptome_bulk/hg19.refGene_CUX1_CASP_diff.gtf


#grab base name of the fastq files
base=`basename $fq_F .fastq.gz`
echo "Sample name is $base"


#specify the number of cores to use for various downstream analysis
cores=8


#create all the output directories
# The -p option means mkdir will create the whole path if it does not exist and refrain from complaining if it does exist

mkdir -p $project_folder/output/STAR
mkdir -p $project_folder/output/bigwigs
mkdir -p $project_folder/output/featureCounts

# set up output filenames and locations

align_out=$project_folder/output/STAR/${base}_
samtools_q30_in=$project_folder/output/STAR/${base}_Aligned.sortedByCoord.out.bam
samtools_q30_out=$project_folder/output/STAR/${base}_Aligned.sortedByCoord.out.q30.bam

bamCoverage_in=$project_folder/output/STAR/${base}_Aligned.sortedByCoord.out.q30.bam
bamCoverage_out=$project_folder/output/STAR/${base}_Aligned.sortedByCoord.out.q30.bw



#run the jobs

#genome alignment
echo "Run STAR aligner"

module load gcc/6.2.0
module load STAR/2.6.1d

cd $project_folder/input

STAR --runThreadN $cores \
    --genomeDir $GENOME_DIR \
    --readFilesIn $fq_F $fq_R \
    --outFileNamePrefix $align_out \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate 



#bam filtering
echo "Run samtools filter and index"

module load gcc/6.2.0
module load samtools/1.10

samtools view -@ 8 -bh -q 30 $samtools_q30_in -o $samtools_q30_out
samtools index $samtools_q30_out

#generate bigwigs
#you need deepTools for this, and deepTools could be called on as long as you loaded the correct python version, so you don't need to load deepTools separately.
echo "generating bigwig files"

module load gcc/6.2.0
module load python/3.8.1

bamCoverage -b $bamCoverage_in -o $bamCoverage_out

#count RNA reads
echo "Run featureCounts"

module load gcc/6.2.0
module load intel/2017
module load subread/1.5.3

featureCounts -T 8 \
	-a $GTF \
	-p \
	-t exon \
	-o $project_folder/output/featureCounts/featurecounts.txt \
	$project_folder/output/STAR/*.q30.bam


done

