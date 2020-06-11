# SYNOPSIS

`assembleAll.pl`

Assembles raw reads into assembly and then
predicts genes.

# Software requirements

In version 2, Shovill was added and many requirements were added

* Shovill - has many dependencies itself
  * BWA
  * Flash
  * Java
  * Lighter
  * Mash
  * Megahit (not used but checked when shovill runs)
  * pigz
  * Pilon
  * Samclip
  * Samtools
  * Seqtk
  * Skesa
  * Trimmomatic
  * Velvet  (not used but checked when shovill runs)
* CG-Pipeline
* Prodigal

# Algorithm

Assembles genomes with Skesa, predicts genes with prodigal,
then produces assembly and gene prediction metrics.

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


