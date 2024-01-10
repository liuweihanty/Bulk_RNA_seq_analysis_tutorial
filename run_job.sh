#!/bin/bash

#specify the input fastq files
fq_F=$1
fq_R=$2
project_dir=$3

#specify the location of reference genome directory
GENOME_DIR=/gpfs/data/mcnerney-lab/reference_genomes_CUX1_CASP_diff/human/hg19_refgenome

#specify the input GTF file
GTF=/gpfs/data/mcnerney-lab/reference_genomes_CUX1_CASP_diff/GTF_files/human/hg19.refGene_CUX1_CASP_diff.gtf


#grab base name of the fastq files
base=`basename $fq_F .fastq.gz`
echo "Sample name is $base"


#specify the number of cores to use for various downstream analysis
cores=8


#create all the output directories
# The -p option means mkdir will create the whole path if it does not exist and refrain from complaining if it does exist

mkdir -p $project_dir/output/STAR
mkdir -p $project_dir/output/bigwigs
mkdir -p $project_dir/output/featureCounts

# set up output filenames and locations

align_out=$project_dir/output/STAR/${base}_
samtools_q30_in=$project_dir/output/STAR/${base}_Aligned.sortedByCoord.out.bam
samtools_q30_out=$project_dir/output/STAR/${base}_Aligned.sortedByCoord.out.q30.bam

bamCoverage_in=$project_dir/output/STAR/${base}_Aligned.sortedByCoord.out.q30.bam
bamCoverage_out=$project_dir/output/bigwigs/${base}_Aligned.sortedByCoord.out.q30.bw



#run the jobs

#genome alignment
echo "Run STAR aligner"

module load gcc/11.3.0
module load gcc/12.1.0
module load intel/2022.2
module load llvm/14.0.5
module load star/2.7.10b

cd $project_dir/input

STAR --runThreadN $cores \
    --genomeDir $GENOME_DIR \
    --readFilesIn $fq_F $fq_R \
    --outFileNamePrefix $align_out \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate 



#bam filtering
echo "Run samtools filter and index"

module load gcc/12.1.0
module load intel/2022.2
module load llvm/14.0.5
module load samtools/1.17

samtools view -@ 8 -bh -q 30 $samtools_q30_in -o $samtools_q30_out
samtools index $samtools_q30_out

#generate bigwigs
#you need deepTools for this, and deepTools could be called on as long as you loaded the correct python version, so you don't need to load deepTools separately.
echo "generating bigwig files"

module load gcc/12.1.0
module load intel/2022.2
module load llvm/14.0.5
module load bedtools/2.30.0
module load python/3.10.5

bamCoverage -b $bamCoverage_in -o $bamCoverage_out

#count RNA reads
echo "Run featureCounts"

module load gcc/12.1.0
module load intel/2022.2
module load llvm/14.0.5
module load subread/2.0.5

featureCounts -T 8 \
	-a $GTF \
	-p \
	-t exon \
	-o $project_dir/output/featureCounts/featurecounts.txt \
	$project_dir/output/STAR/*.q30.bam


done

