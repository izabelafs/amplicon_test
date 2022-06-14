#####Load required libraries #####
library(renv)
library(BiocManager)
library(VariantAnnotation)
library(valr)
library(SEQprocess)
library(tidyverse)
library(GenomicAlignments)
library(ggplot2)
library(ggseqlogo)
library(drake)
library(here)

#Binary paths
binary <- "/opt/homebrew/bin"
cat <- "/bin/cat"
bwa <- file.path(binary,"bwa")
samtools <- file.path(binary,"samtools")
bcftools <- file.path(binary,"bcftools")
bedtools <- file.path(binary, "bedtools")
tabix <- file.path(binary, "tabix")

#Paste fasta files related to the same gene
paste <- c("amplicon_1","amplicon_2","amplicon_12")

#Align each amplicon read to the reference genome
amplicon <- c("amplicon_12","amplicon_3","amplicon_3.test_2","amplicon_3.test_3")

#Compute the .bed files to access the coordinates
list_genes <- c("gene1","gene2")