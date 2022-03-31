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
counts_file <- #Read in feature counts files using read.delim
counts_file <- read.delim(file = "", skip = 1,
                              col.names = c("geneID", "Chr", "Start", "End", "Strand", "Length", "R702R_KO_6F_rep2", "R702R_KO_6G_rep1", "R702R_KO_6G_rep2", "K700E_TWT1_rep1", "K700E_TWT1_rep2", "K700E_TWT2_rep1", "K700E_TWT2_rep2", "K700E_KO_1H_rep1", "K700E_KO_1H_rep2", "K700E_KO_3A_rep1", "R702R_TWT1_rep1", "K700E_KO_3A_rep2", "K700E_KO_3B_rep1", "K700E_KO_3B_rep2", "K700E_KO_6B_rep1", "K700E_KO_6B_rep2", "R702R_TWT1_rep2", "R702R_TWT2_rep1", "R702R_TWT2_rep2", "R702R_KO_2E_rep1", "R702R_KO_2E_rep2", "R702R_KO_5F_rep1", "R702R_KO_5F_rep2", "R702R_KO_6F_rep1"),
                              row.names = "geneID") %>%
                              dplyr::select(!c("geneID", "Chr", "Start", "End", "Strand", "Length"))
sample_data <- read_tsv(file = "/path/to/sample_data.txt")
