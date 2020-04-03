# SYNOPSIS

`sn_assemblyWorkflow_init.pl`

Initializes an assembly-only workflow so that other plugins can read the proper data.

# Software requirements

* CG-Pipeline
* Prodigal

# Algorithm

Copies each assembly in root SneakerNet-formatted folder into `$dir/SneakerNet/assembly/$sample`.
Runs assembly metrics.
Runs prodigal on all assemblies.
Runs gene prediction metrics.

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


