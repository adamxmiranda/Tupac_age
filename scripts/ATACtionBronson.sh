#!/bin/bash

################################################################################
###This is the script for the wRapper function for ATAC-seq processing##########
#############Better known as ATACtionBronson###################################
################################################################################

usage(){
echo "
ATACtionBronson
Description:
Personalized wRapper function for automating ATAC-seq Processing

How to use this function:
ATACtionBronson.sh path=/path/to/fastq/files

Requirements:
A build of the bbmap package in your .bashrc

Parameters:
path=/path/to/fastq/files         Specify the path to the fastq files. This wRapper
                                  will create multiple subdirectories within this
                                  directory as part of the processing.

"
}
if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi
#ARGUMENT parsing function from Stackexchange
#shout out JRichardsz
for ARGUMENT in "$@"
do
   KEY=$(echo $ARGUMENT | cut -f1 -d=)

   KEY_LENGTH=${#KEY}
   VALUE="${ARGUMENT:$KEY_LENGTH+1}"

   export "$KEY"="$VALUE"
done
# Print variables
echo "path = $path"
if [ -z "$path" ]; then
   echo "Path does not exist"
   exit;
 else
   echo "Path exists"

fi

#create subdirectories for staged outputs
if [ ! -d "$path/raw_reads" ]; then
  mkdir $path/raw_reads
else
  echo "$path/raw_reads already exists"
fi
mv $path/*.fastq $path/raw_reads
if [ ! -d "$path/trimmed_reads" ]; then
  mkdir $path/trimmed_reads
else
  echo "$path/trimmed_reads already exists"
fi
TRIM_DIR="$path/trimmed_reads"
if [ ! -d "$path/trimmed_reads/fastQC" ]; then
  mkdir $path/trimmed_reads/fastQC
else
  echo "$path/trimmed_reads/fastQC already exists"
fi
FASTQC_DIR="$path/trimmed_reads/fastQC"
if [ ! -d "$path/mapped_reads" ]; then
  mkdir $path/mapped_reads
else
  echo "$path/mapped_reads already exists"
fi
MAP_DIR="$path/mapped_reads"
if [ ! -d "$path/filtered_reads" ]; then
  mkdir $path/filtered_reads
else
  echo "$path/filtered_reads already exists"
fi
FILT_DIR="$path/filtered_reads"

#Loop to perform trimming and quality control
#Loop only works for files named with the same convention as what VANTAGE uses
#change to the raw reads directory
cd $path/raw_reads
#load the modules for this step
module load GCC/5.4.0-2.26
module load cutadapt/1.9.1-Python-3.5.2
module load FastQC/0.11.9
SUFF_1="R1_001.fastq.gz"
SUFF_2="R2_001.fastq.gz"
for filename in *R1_001.fastq.gz
do
  echo "begin trimming ${filename}"
  base=$(basename $filename R1_001.fastq.gz)
  echo ${base}
  trim_galore --fastqc --fastqc_args "--outdir ${FASTQC_DIR}" \
  --paired --retain_unpaired \
  --output_dir ${TRIM_DIR} \
  ${base}${SUFF_1} ${base}${SUFF_2}
  echo "end trimming ${base}"
done

#loop for bbmap repair, removes singletons and sorts .fq files
#satisfies whiny mapping alignment algorithms
cd ${TRIM_DIR}
for filename in *R1_001_val_1.fq.gz
do
  base=`basename $filename _S1_L005_R1_001_trimmed.fq`
  echo "Can we fix it?"
  echo "repair $base"
  #use bbmap_repair to filter mismatched reads
  repair.sh in1=${base}R1_001_val_1.fq.gz in2=${base}R2_001_val_2.fq.gz \
  out1=${base}_Mate1_repaired.fq out2=${base}_Mate2_repaired.fq \
  outs=${base}_singletons.fq repair
  echo "Yes we can!"
  echo "repaired $base"
done

#load BWA modules
module load GCC/6.4.0-2.28
module load BWA/0.7.17
module load SAMtools/1.6
#In the case I need to make another index
# cd /home/mirandax/bwa_index/hg38
# bwa index -a bwtsw -p $PREFIX $REF
#BWA mapping loop
REF="/data/park_lab/adam/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
for filename in *_Mate1_repaired.fq
do
  base=`basename $filename _Mate1_repaired.fq`
  echo $filename
  echo "let's map dis tang"
  bwa mem -t 12 $REF ${TRIM_DIR}/${base}_Mate1_repaired.fq ${TRIM_DIR}/${base}_Mate2_repaired.fq | \
  samtools view -@ 12 -Sb - > ${MAP_DIR}/${base}_bwa.bam
  echo "end map"
done

#Filtering loop
#Filtering for mapQ>40 (bwa max is 60)
#Removes mitochondrial DNA
#Removes blacklisted sites
module load GCC/8.2.0 SAMtools/1.9
cd ${MAP_DIR}
for filename in *_bwa.bam
do
  base=`basename $filename _bwa.bam`
  echo "convert sam to bam & filter for MAPQ > 40 & sort bam"
  samtools view -@ 12 -S -b -q 40 ${base}_bwa.bam > ${base}_bwa_mapq40.bam
  samtools sort -@ 12 -o ${base}_bwa_mapq40_sorted.bam ${base}_bwa_mapq40.bam
  echo " index sorted bam files"
  samtools index -@ 12 -b ${base}_bwa_mapq40_sorted.bam ${base}_bwa_mapq40_sorted.bam.bai
  echo "remove mtDNA reads"
  samtools view -@ 12 -b ${base}_bwa_mapq40_sorted.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
  chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > ${base}_bwa_no_chrM.bam
  echo "Index no_ChrM files"
  samtools index -@ 12 -b ${base}_bwa_no_chrM.bam  ${base}_bwa_no_chrM.bam.bai
  echo "remove blacklisted regions"
  samtools view -@ 12 -b -L /data/hodges_lab/hg38_genome/hg38.blacklist.bed \
  -U ${base}_bwa_filtered.unsorted.bam ${base}_bwa_no_chrM.bam > ${base}_bwa_blacklisted.bam
  samtools sort -@ 12 ${base}_bwa_filtered.unsorted.bam > ${FILT_DIR}/${base}_bwa_filtered.bam
  echo "Index filtered files"
  samtools index -@ 12 -b ${FILT_DIR}/${base}_bwa_filtered.bam ${FILT_DIR}/${base}.filtered.bam.bai
done

#remove sequencing duplicates
module load picard/2.18.27
cd ${FILT_DIR}
PICARD="$EBROOTPICARD/picard.jar"
for filename in *_filtered.bam
do
  echo "removing sequencing duplicates"
  echo ${filename}
  base=`basename $filename _filtered.bam`
  java -jar $PICARD MarkDuplicates I=${base}_filtered.bam O=${base}.unique.bam \
  M=${base}_marked_dup_metrics-all.txt REMOVE_DUPLICATES=TRUE
  echo "dupes removed"
done
