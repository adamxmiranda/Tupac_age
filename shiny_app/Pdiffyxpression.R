# Load packages
library(shiny)
library(shinythemes)
library(tidyverse)
library(janitor)
library(scales)
library(DESeq2)
library(biomaRt)
library(ggrepel)
library(apeglm)
library(genefilter)
library(pheatmap)
library(cluster)
library(factoextra)
library(RColorBrewer)

# Data Input
## This part will need to be customized by each person using the app to add their own data
## What you will need is a counts file (The output from RNA_DMC will work just fine) and a sample_data file
### The sample_data file should include each row as a sample name !IN THE SAME ORDER AS IN THE COUNTS FILE!
### Each column should include the corresponding details of an independent variable
#### ex. Genotype column and each sample with WT or Mutant in the cell
#If using the output from featureCounts/
counts_file <- read
sample_data <- read_tsv(file = "/path/to/sample_data.txt")
