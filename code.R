source(here::here("functions.R"))
source(here::here("config.R"))

###### System calls ######
#Index genomeFile using bwa 
if (!file.exists("data/genome.fa.fai"))
  system2(bwa, args = "index data/genome.fa" )

#Generate plan to keep track of the data manipulations after indexing of reference genome 
system_plan <-
  drake_plan(
    paste_amplicon(paste),
    align_amplicon(amplicon),
    clean_up(amplicon),
    sort_aligned(amplicon),
    index_sorted(amplicon)
  )

###### R calls ######
data_plan <- drake_plan(
  genes = read_genes(list_genes),
  amp.coord =
    read_bed("data/amplicon_coordinates.bed"),
  hits = find_overlaps(genes, amp.coord),
  seq_12 = output_seq("amplicon_12"),
  seq_3 = output_seq("amplicon_3"),
  seq_3_test2 =
    output_seq("amplicon_3_test2"),
  seq_3_test3 =
    output_seq("amplicon_3_test3"),
  result_12 =
    find_supporting(seq_12, "amplicon_12",hits),
  result_3 =
    find_supporting(seq_3, "amplicon_3",hits),
  test2 =
    find_supporting(seq_3_test2, "amplicon_3",hits),
  test3 =
    find_supporting(seq_3_test3, "amplicon_3",hits),
  results = format(result_12, result_3),
  test.results = format(test2, test3),
  frequency =
    bind_rows(return_n_count(seq_12), return_n_count(seq_3))
)

make(system_plan)
make(data_plan)
loadd(data_plan$target)
### Optional plots to summarize mismatches ###
#BarPlot
bar <-
  ggplot(results, aes(fill = gene, y = frequency, x = index)) + geom_bar(position =
                                                                           "stack", stat = "identity")
print(bar)

#Matrix of nucleotide counts to generate ggseqlogo
frequency_can <- frequency %>% select(starts_with("fq")) %>% as.matrix() %>% t()
frequency_mm <- frequency %>% select(ends_with("fq")) %>% as.matrix() %>% t()
rownames(frequency_can) <-  c('A','C','G','T')
colnames(frequency_can) <- c("amplicon_12","amplicon_3")
rownames(frequency_mm) <-  c('A','C','G','T')
colnames(frequency_can) <- c("amplicon_12","amplicon_3")

#ggseqlogo
ggseqlogo(frequency_can, method='custom', seq_type='dna') + ylab('Canonical amplicons')
ggseqlogo(frequency_mm, method='custom', seq_type='dna') + ylab('Mutated amplicons')
                           