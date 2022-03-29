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
A build of the Subread package in your .bashrc
I use


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
