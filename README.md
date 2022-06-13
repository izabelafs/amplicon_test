# Problem Description
Several input files are provided (see detailed file explanation below).
The candidate should read the input files and process them.
Each amplicon fasta input files is a collection of 200 DNA sequences of the same length, whose coordinates are provided in the amplicon_coordinates.bed file.
Some mutations have been introduced in the sequences.
Amplicon1 and amplicon2 are partially overlapping amplicons, so some mutations might be present in both amplicons and therefore they should be combined.
Amplicon 3 is an independent one.
At a high level, the candidate should:
1)	Read sequences from the amplicon_*.fa.gz file;
2)	Identify mutations;
3)	Report frequency of each mutation.
The fasta reference file is provided as genome.fa.gz.
See next paragraph “Questions” for a detailed list of all questions and tasks

##Read amplicon files
1)	Read input files for the three amplicons files amplicon_1.fa.gz, amplicon_2.fa.gz
and amplicon_3.fa.gz and report to file a table containing:
a)	Gene name
b)  Nature of the mutation (e.g. C->T etc)
c)	frequency of the mutation
d)	number of amplicons supporting the mutation
f)	reference sequence for the amplicon
g)	mutated sequence for the amplicon

They will most likely need previous alignment to the reference genome \genome.fa.gz\

2)	Load amplicon file amplicon_3.test_2.fa, and report in a similar format as in 1)
all the mutations in the amplicon
3)	Extra Challenge: display in graphical format one or more of the previously identified mutation.
Suggestions: the candidate might want to report the frequency of the mutation with a plot or bar plot;
or use a sequence logo; or summarize the average mutation frequency per amplicon; etc. Format and solution are open

##Primer presence

Amplicon libraries are prepared by PCR amplification of a specific target,
for example the V4 hypervariable region of the bacterial 16S rRNA gene.
All reads from this type of library are expected to be nearly identical.

Expected results are:

Extremely biased per base sequence content
Extremely narrow distribution of GC content
Very high sequence duplication levels
Abundance of overrepresented sequences
In cases where the PCR target is shorter than the read length,
the sequence will read through into adapters
