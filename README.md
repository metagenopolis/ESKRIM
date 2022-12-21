# ESKRIM: EStimate with K-mers the RIchness in a Microbiome #


ESKRIM is a reference-free tool that compares microbial richness in shotgun metagenomic samples by counting k-mers

## Installation via pip ##
```
pip install eskrim
```

## Usage ##
### Basic usage ###
In this example, k-mer richness in a sample (sample1) consisting in two paired-end runs (run1 and run2) is computed.\
Forward fastq files are taken as input. Results are saved in the file *sample1.eskrim_stats.tsv*

```
eskrim -i sample1.run1_1.fastq.gz sample1.run2_1.fastq.gz -n sample1 -s sample1.eskrim_stats.tsv
```
**Quality control** (adapters removal, read trimming) and **contaminant removal** (reads from the host genome) should be performed before using ESKRIM.

Run ESKRIM similarly for each sample to be compared. All TSV output files can be merged manually.

### Advanced usage ###

#### Adjusting target read count for subsampling ####
Depending on the sequencing depth, the target number of reads to randomly draw from each sample (default = 10M) can be adjusted with the *-r* parameter.
```
eskrim -i sample1.run1_1.fastq.gz sample1.run2_1.fastq.gz -n sample1 -s sample1.eskrim_stats.tsv -r 5000000
```

#### Adjusting read length ###
All reads are trimmed to a given length (default = 80) because read length can vary between samples.\
This length can be changed with the *-l* parameter.
```
eskrim -i sample1.run1_1.fastq.gz sample1.run2_1.fastq.gz -n sample1 -s sample1.eskrim_stats.tsv -l 100
```

#### Reproducibility ###
ESKRIM ensures reproducibility when using the same random number generator seed (default = 0).\
To make read subsampling vary across executions, the parameters --seed can be used.
```
eskrim -i sample1.run1_1.fastq.gz sample1.run2_1.fastq.gz -n sample1 -s sample1.eskrim_stats.tsv --seed 1234
```

## Interpreting the output file ##
ESKRIM saves the results in a TSV file consisting in several columns (*-s* parameter).
* *sample_name* : sample name specified with *-n* parameter.
* *total_num_reads* : number of reads in the sample before subsampling.
* *num_Ns_reads_ignored* : number of reads with undetermined bases that were discarded.
* *num_too_short_reads_ignored* : number of reads with undetermined bases that were discarded.
* *target_num_reads* : target number of reads to draw during the subsampling step.
* *num_selected_reads* : number of reads actually drawn after subsampling.
* *read_length* : length at which reads were trimmed (*-l* parameter).
* *kmer_length* : length of counted k-mers (*-k* parameter).
* *num_distinct_kmers* : number of distinct kmers in subsampled reads.
* *num_solid_kmers* : number of kmers seen at least twice.
* *num_mercy_kmers* : number of non-solid kmers occuring in a read where all k-mers are not solid.\
__This value is a proxy to compare microbial richness between samples.__

**WARNING**: Do not consider results when *num_selected_reads* is strictly lower than *target_num_reads*.\
In this case, ignore the samples concerned or decrease the number of reads to be drawn randomly (*-r* parameter).

## Authors ##

* Florian Plaza Oñate: florian.plaza-onate@inrae.fr
* Emmanuelle Le Chatelier: emmanuelle.le-chatelier@inrae.fr

