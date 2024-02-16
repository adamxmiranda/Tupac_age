#!/bin/bash
#####################################################################
###This is the script for the wRapper function of RNA-seq analysis###
#############Better known as RNA_DMC#################################
#####################################################################

usage(){
echo "
RNA_DMC
Description:
Personalized wRapper function for automating RNA-seq Processing

How to use this function:
RNA_DMC.sh path=/path/to/fastq/files

Requirements:
A build of the Subread package in your .bashrc
I use
/subread-2.0.0-Linux-x86_64/bin

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
# #create subdirectories for staged outputs
if [ ! -d "$path/raw_reads" ]; then
  mkdir $path/raw_reads
else
  echo "$path/raw_reads already exists"
fi
mv $path/*.fastq.gz $path/raw_reads
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
if [ ! -d "$path/counts_file" ]; then
  mkdir $path/counts_file
else
  echo "$path/counts_file already exists"
fi
COUNT_DIR="$path/counts_file"
Loop to perform trimming and quality control
Loop only works for files named with the same convention as what VANTAGE uses
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

# Loop to perform mapping
#load modules
module load GCC/6.4.0-2.28
module load STAR/2.5.4b
#hg38 index on park_lab
INDEX="/data/park_lab/STAR_index"
TRIM_DIR="$path/trimmed_reads"
MAP_DIR="$path/mapped_reads"
FILT_DIR="$path/filtered_reads"
COUNT_DIR="$path/counts_file"
#suffixes for trimmed files
TRIM_SUFF_1="R1_001_val_1.fq.gz"
TRIM_SUFF_2="R2_001_val_2.fq.gz"
cd ${TRIM_DIR}
#mapping loop for trimmed reads created in the first loop
for filename in ${TRIM_DIR}/*${TRIM_SUFF_1}
do
  echo "mapping ${filename}"
  base=$(basename $filename R1_001_val_1.fq.gz)
  echo ${base}
  STAR --runMode alignReads \
  --runThreadN 12 \
  --genomeDir ${INDEX} \
  --readFilesCommand zcat \
  --readFilesIn ${base}${TRIM_SUFF_1} ${base}${TRIM_SUFF_2} \
  --outFileNamePrefix ${MAP_DIR}/${base}.sam
  echo "finished mapping ${base}"
done

#Loop to perform Filtering and sorting of mapped reads
#load modules
module load GCC/5.4.0-2.26
module load SAMtools/1.5
cd ${MAP_DIR}
for filename in ${MAP_DIR}/*.samAligned.out.sam
do
  echo ${filename}
  base=$(basename $filename .samAligned.out.sam)
  echo $base
  echo "${filename}: convert sam to bam & filter for MAPQ > 30 & sort bam"
  samtools view -@ 8 -S -b -q 30 ${filename} | samtools sort -@ 8 -n -o ${FILT_DIR}/${base}.bam
  echo "filtered and sorted ${base}"
done

#Feature Counts
cd ${FILT_DIR}
featureCounts -F GTF -a /data/park_lab/STAR_index/gtf/genes.gtf \
              -G /data/park_lab/STAR_index/fasta/genome.fa \
              -o ${COUNT_DIR}/featureCounts_all_samples.txt \
              -T 12 -p \
              -O -t exon ${FILT_DIR}/*.bam
