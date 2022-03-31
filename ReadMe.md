# Tupac_age
![Logo](https://github.com/adamxmiranda/Tupac_age/blob/'progress'/logo/Tupac_age-logo.svg?raw=true)
This is the Kendrick La_markdown of the Tu-pac-age of wRapper functions to automate my work for me. Herein will be many rap-centric (Big) puns to make the tedium of writing this much more tolerable.

These functions are personalized for use on Vanderbilt's ACCRE cluster. Efforts to generalize this for wider use likely won't happen.

## RNA-seq Processing
```bash
RNA_DMC.sh
```
This function contains all of the necessary steps to go from .fastq to featureCounts.
Make your analysis a little less Tricky!
Tricky, Tricky, Tricky

This wrapper makes use of functions from:
[TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
[SAMtools](http://www.htslib.org/)
[STAR](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)
[Subread](http://subread.sourceforge.net/)

```bash
#Usage Example
sh RNA_DMC.sh /path/to/fastq/files
```
### Further RNA-seq Analysis
Under construction is a R shiny web app wRapper of DESeq2 and clusterProfiler for further analysis of RNA-seq

```bash
Pdiffyxpression.R
```

## ATAC-seq Processing
```bash
ATACtionBronson.sh
```
With this wRapper you don't have to be "Blue". This wrapper processes ATAC-seq data from .fastq through mapping and filtering. This function wraps functions from:
[TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
[bbmap](https://sourceforge.net/projects/bbmap/)
[BWA](http://bio-bwa.sourceforge.net/)
[SAMtools](http://www.htslib.org/)
[Picard](https://broadinstitute.github.io/picard/)

Functionality for peak calling is under construction

```bash
#Usage Example
sh ATACtionBronson.sh /path/to/fastq/files
```
