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

#Binary paths
binary <- "/opt/homebrew/bin"
cat <- "/bin/cat"
bwa <- file.path(binary,"bwa")
samtools <- file.path(binary,"samtools")
bcftools <- file.path(binary,"bcftools")
bedtools <- file.path(binary, "bedtools")
tabix <- file.path(binary, "tabix")

###### System calls ######
#Index genomeFile using bwa 
#system2(bwa, args = "index data/genome.fa" )

#Align each amplicon read to the reference genome
system2(cat, args = "data/amplicon_1.fa data/amplicon_2.fa > data/amplicon_12.fa")
system2(bwa, args = "mem data/genome.fa data/amplicon_12.fa > data/amplicon_12.sam")
system2(bwa, args = "mem data/genome.fa data/amplicon_3.fa > data/amplicon_3.sam")
system2(bwa, args = "mem data/genome.fa data/amplicon_3.test_2.fa > data/amplicon_3.test_2.sam")
system2(bwa, args = "mem data/genome.fa data/amplicon_3.test_3.fa > data/amplicon_3.test_3.sam")


#Clean up read pairing information and flags:
system2(samtools, args = "fixmate -O bam data/amplicon_12.sam data/amplicon_12_fixmate.bam")
system2(samtools, args = "fixmate -O bam data/amplicon_3.sam data/amplicon_3_fixmate.bam")
system2(samtools, args = "fixmate -O bam data/amplicon_3.test_2.sam data/amplicon_3_test2.fixmate.bam")
system2(samtools, args = "fixmate -O bam data/amplicon_3.test_3.sam data/amplicon_3_test3.fixmate.bam")

#Sort each aligned .sam  
system2(samtools, args = "sort -O bam -o data/amplicon_12_sorted.bam -T /tmp/lane_temp data/amplicon_12_fixmate.bam")
system2(samtools, args = "sort -O bam -o data/amplicon_3_sorted.bam -T /tmp/lane_temp data/amplicon_3_fixmate.bam")
system2(samtools, args = "sort -O bam -o data/amplicon_3_test2_sorted.bam -T /tmp/lane_temp data/amplicon_3_test2.fixmate.bam")
system2(samtools, args = "sort -O bam -o data/amplicon_3_test3_sorted.bam -T /tmp/lane_temp data/amplicon_3_test3.fixmate.bam")

#Index each BAM using samtools:
system2(samtools, args = "index data/amplicon_12_sorted.bam")
system2(samtools, args = "index data/amplicon_3_sorted.bam")
system2(samtools, args = "index data/amplicon_3_test2_sorted.bam")
system2(samtools, args = "index data/amplicon_3_test3_sorted.bam")

##Compute the .bed files to access the coordinates
gene1 <- read_bed("data/gene1.bed")
colnames(gene1) <- c("ensembl","chr","start","end","score")
gene2 <- read_bed("data/gene2.bed")
colnames(gene2) <- c("ensembl","chr","start","end","score")
genes <- bind_rows(gene1,gene2)
amp.coord <- read_bed("data/amplicon_coordinates.bed")

#Generate GRanges for each to make the track comparison
genes <- makeGRangesFromDataFrame(genes,
                                  keep.extra.columns=T,
                                  seqnames.field=c("seqnames", "seqname",
                                                   "chromosome", "chrom",
                                                   "chr", "chromosome_name",
                                                   "seqid"),
                                  start.field="start",
                                  end.field="end")
names(genes) <- genes@elementMetadata$ensembl
amp.coord <- makeGRangesFromDataFrame(amp.coord,
                                      keep.extra.columns=T,
                                      seqnames.field=c("seqnames", "seqname",
                                                       "chromosome", "chrom",
                                                       "chr", "chromosome_name",
                                                       "seqid"),
                                      start.field="start",
                                      end.field="end")

#Make sure the seqlevels are the same
seqlevels(amp.coord) <- "chr7"
seqlevels(genes) <- "chr7"

##Find overlaps between the bed files 
overlaps <- findOverlaps(amp.coord, genes)
hits_ov <- genes[subjectHits(overlaps)]

#Read amplicons to extract the sequences
amplicon_12 <-  scanBam("data/amplicon_12_sorted.bam")
amplicon_3 <-  scanBam("data/amplicon_3_sorted.bam")
amplicon_3_test2 <-  scanBam("data/amplicon_3_test2_sorted.bam")
amplicon_3_test3 <-  scanBam("data/amplicon_3_test3_sorted.bam")

seq_12 <- amplicon_12[[1]][["seq"]]
seq_3 <- amplicon_3[[1]][["seq"]]
seq_3_test2 <- amplicon_3_test2[[1]][["seq"]]
seq_3_test3 <- amplicon_3_test3[[1]][["seq"]]

#Define function to count mutations
find_supporting <- function(seq,amplicon) {
  seq <- as.data.frame(as.character(seq))
  diff <- unique(seq)
  result <-
    data.frame(
      gene = 1:nrow(diff),
      index = 1:nrow(diff),
      ref = 1:nrow(diff),
      mut = 1:nrow(diff),
      can = 1:nrow(diff),
      mm = 1:nrow(diff),
      n_can = 1:nrow(diff),
      n_mm = 1:nrow(diff),
      fq_a = 1:nrow(diff),
      fq_c = 1:nrow(diff),
      fq_g = 1:nrow(diff),
      fq_t = 1:nrow(diff),
      a_mm_fq = 1:nrow(diff),
      c_mm_fq = 1:nrow(diff),
      g_mm_fq = 1:nrow(diff),
      t_mm_fq = 1:nrow(diff)
    )
  for (i in 1:nrow(diff)) {
    if (i %% 2 == 0) {
      next
    } else {
      result$index[i] <- amplicon
      result$ref[i] <- diff[i,]
      result$mut[i] <- diff[i + 1,]
      sub1[i] <- strsplit(diff[i, ], "")[[1]]
      sub2[i] <- strsplit(diff[i + 1, ], "")[[1]]
      result$can[i] <- sub1[sub1 != sub2]
      result$mm[i] <- sub2[sub1 != sub2]
      result$fq_a[i] <- sum(sub1 == "A")/length(sub1)
      result$fq_c[i] <- sum(sub1 == "C")/length(sub1)
      result$fq_g[i] <- sum(sub1 == "G")/length(sub1)
      result$fq_t[i] <- sum(sub1 == "T")/length(sub1)
      result$a_mm_fq[i] <- sum(sub2 == "A")/length(sub2)
      result$c_mm_fq[i] <- sum(sub2 == "C")/length(sub2)
      result$g_mm_fq[i] <- sum(sub2 == "G")/length(sub2)
      result$t_mm_fq[i] <- sum(sub2 == "T")/length(sub2)
      result$n_can[i] <- length(which(seq == diff[i, ]))
      result$n_mm[i] <- length(which(seq == diff[(i + 1), ]))
      if (amplicon == "amplicon_12"){
        result$gene <- hits_ov$ensembl[1]
      }
      else if (amplicon == "amplicon_3"){
        result$gene <- hits_ov$ensembl[3]
      }
    }
  }
  return(result)
}

result_12 <- find_supporting(seq_12,"amplicon_12")
result_3 <- find_supporting(seq_3,"amplicon_3")

test2 <- find_supporting(seq_3_test2,"amplicon_3")
test3 <- find_supporting(seq_3_test3,"amplicon_3")

#Format the results accordingly
results <-
  rbind(result_12, result_3) %>% subset(n_mm > 10) %>% mutate(frequency = n_mm/n_can,
                                                              type = paste(can, sep = ">", mm)) %>% select(-mm,-can)


test.results <-
  rbind(test2, test3) %>% subset(n_mm > 10) %>% mutate(frequency = n_mm /
                                                         n_can,
                                                       type = paste(can, sep = ">", mm)) %>% select(-starts_with("fq"),-mm,-can)

#Optional plots to summarize mismatches

bar <- ggplot(results, aes(fill=gene, y=frequency, x=index)) + 
  geom_bar(position="stack", stat="identity")

print(bar)

frequency_can <- results %>% select(starts_with("fq")) %>% as.matrix() %>% t()
frequency_mm <- results %>% select(ends_with("fq")) %>% as.matrix() %>% t()

rownames(frequency_can) <-  c('A','C','G','T')
rownames(frequency_mm) <-  c('A','C','G','T')

ggseqlogo(frequency_can, method='custom', seq_type='dna') + ylab('Canonical amplicons')
ggseqlogo(frequency_mm, method='custom', seq_type='dna') + ylab('Mutated amplicons')
                           