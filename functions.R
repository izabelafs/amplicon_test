##Define data-handle functions to process the amplicon data

#Paste fasta files related to the same gene
paste_amplicon <- function(paste){
  system2(
    cat,
    args = paste0(
      "data/",
      paste[1],
      ".fa",
      " ",
      "data/",
      paste[2],
      ".fa",
      " ",
      ">"," ",
      "data/",paste[3],".fa"
    )
  )
}

#Align each amplicon read to the reference genome
align_amplicon <- function(amplicon) {
  for (i in 1:length(amplicon)) {
    system2(
      bwa,
      args = paste0(
        "mem data/genome.fa",
        " ",
        "data/",
        amplicon[i], ".fa",
        " ",
        ">",
        " ",
        "data/",
        amplicon[i],
        ".sam"
      )
    )
  }
}

#Clean up read pairing information and flags:
clean_up <- function(amplicon) {
  for (i in 1:length(amplicon)) {
    system2(
      samtools,
      args = paste0(
        "fixmate -O bam ",
        "data/",
        amplicon[i],
        ".sam ",
        "data/",
        amplicon[i],
        "_fixmate.bam"
      )
    )
  }
}

#Sort each aligned .sam 
sort_aligned <- function(amplicon) {
  for (i in 1:length(amplicon)) {
    system2(
      samtools,
      args = paste0(
        "sort -O bam -o ",
        "data/",
        amplicon[i],
        "_sorted.bam ",
        "-T /tmp/lane_temp ",
        "data/",
        amplicon[i],
        "_fixmate.bam"
      )
    )
  }
}

#Index each BAM using samtools:
index_sorted <- function(amplicon){
  for (i in 1:length(amplicon)) {
    system2(samtools, args = paste0("index ","data/",amplicon[i],"_sorted.bam"))
  }
}
#Compute the .bed files to access the coordinates
read_genes <- function(genes){
  i = 1:2
  list = read_bed(paste0("data/",genes[i],".bed")) %>% bind_rows()
  colnames(list) <- c("ensembl","chr","start","end","score")
  return(list)
}

#Generate the GRanges objects to get the overlaps between the coordinates
find_overlaps <- function(genes, amp.coord) {
  df1 <- makeGRangesFromDataFrame(
    genes,
    keep.extra.columns = T,
    seqnames.field = c(
      "seqnames",
      "seqname",
      "chromosome",
      "chrom",
      "chr",
      "chromosome_name",
      "seqid"
    ),
    start.field = "start",
    end.field = "end"
  )
  names(df1) <- df1@elementMetadata$ensembl
  df2 <- makeGRangesFromDataFrame(
    amp.coord,
    keep.extra.columns = T,
    seqnames.field = c(
      "seqnames",
      "seqname",
      "chromosome",
      "chrom",
      "chr",
      "chromosome_name",
      "seqid"
    ),
    start.field = "start",
    end.field = "end"
  )
  names(df2) <- df2@elementMetadata$ensembl
  seqlevels(df1) <- seqlevels(df2)
  overlaps <- findOverlaps(df1, df2)
  hits <- df1[queryHits(overlaps)]
  return(hits)
}
#Read amplicons to extract the sequences
output_seq <-  function(amplicon){
  path <- file.path(paste0("data/", amplicon,"_sorted.bam"))
  amplicon <- scanBam(path)
  seq <- amplicon[[1]][["seq"]]
  return(seq)
}

#Define function to count the number of mutated sequences
find_supporting <- function(seq,amplicon,hits) {
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
      n_mm = 1:nrow(diff)
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
      result$n_can[i] <- length(which(seq == diff[i, ]))
      result$n_mm[i] <- length(which(seq == diff[(i + 1), ]))
      if (amplicon == "amplicon_12"){
        result$gene <- hits$ensembl[1]
      }
      else if (amplicon == "amplicon_3"){
        result$gene <- hits$ensembl[3]
      }
    }
  }
  return(result)
}

#Format the results accordingly
format <- function(df1, df2) {
  df <-
    rbind(df1, df2) %>% subset(n_mm > 10) %>% mutate(frequency = n_mm /
                                                                  n_can,
                                                                type = paste(can, sep = ">", mm)) %>% select(-mm, -can)
  return(df)
}
#Define function to count nucleotides
return_n_count <- function(seq) {
  seq <- as.data.frame(as.character(seq))
  diff <- unique(seq)
  df <-
    data.frame(
      fq_a = 1:length(diff),
      fq_c = 1:length(diff),
      fq_g = 1:length(diff),
      fq_t = 1:length(diff),
      a_mm_fq = 1:length(diff),
      c_mm_fq = 1:length(diff),
      g_mm_fq = 1:length(diff),
      t_mm_fq = 1:length(diff)
    )
  for (i in 1:length(diff)) {
    if (i %% 2 == 0) {
      next
    } else {
      sub1[i] <- strsplit(diff[i, ], "")[[1]]
      sub2[i] <- strsplit(diff[i + 1, ], "")[[1]]
      df$fq_a[i] <- sum(sub1 == "A") / length(sub1)
      df$fq_c[i] <- sum(sub1 == "C") / length(sub1)
      df$fq_g[i] <- sum(sub1 == "G") / length(sub1)
      df$fq_t[i] <- sum(sub1 == "T") / length(sub1)
      df$a_mm_fq[i] <- sum(sub2 == "A") / length(sub2)
      df$c_mm_fq[i] <- sum(sub2 == "C") / length(sub2)
      df$g_mm_fq[i] <- sum(sub2 == "G") / length(sub2)
      df$t_mm_fq[i] <- sum(sub2 == "T") / length(sub2)
    }
    return(df)
  }
}
