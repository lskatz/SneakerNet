# SYNOPSIS

`sn_iontorrent_assembleAll.pl`

Assembles ion torrent raw read data and
predicts genes

# Software requirements

* SPAdes
* CG-Pipeline
* Prodigal

# Algorithm

Assembles raw read data with SPAdes with optimal options,
predicts genes with Prodigal,
then reports assembly and gene prediction metrics.

# Outputs

Table with columns

* sample
* genomeLength
* CDS - number of coding sequences
* N50
* longestContig
* numContigs
* avgContigLength
* assemblyScore - log(N50/numContigs * penalty) where a genome is penalized for being different than the expected genome length
* minContigLength
* expectedGenomeLength
* kmer21 - a quantifier for how much 21-mers are duplicated. A marker for overassembly.
* GC

