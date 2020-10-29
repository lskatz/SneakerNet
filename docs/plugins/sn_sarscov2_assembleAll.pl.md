# SYNOPSIS

`sn_sarscov2_assembleAll.pl`

Runs a reference assembly against Wuhan-1
and produces assembly metrics.

# Software requirements

* Bowtie2
* samtools
* bcftools
* cutadapt
* seqtk

# Algorithm

* Downloads the reference genome into your installation directory
* Trims adapters with cutadapt
* Maps reads against genome with bowtie2
* Calls SNPs with samtools
* Filters SNPs
* Consensus fasta from SNPs and the reference genome

# Outputs

Assemblies are found in `SneakerNet/assemblies/[samplename]/*.fasta`

Table with columns

* File - the sample name with some file extensions
* genomeLength
* N50
* longestContig
* numContigs
* avgContigLength
* assemblyScore - log(N50/numContigs * penalty) where a genome is penalized for being different than the expected genome length
* minContigLength
* expectedGenomeLength
* kmer21 - a quantifier for how much 21-mers are duplicated. A marker for overassembly.
* GC
* effectiveCoverage - the number of mapped base pairs against the assembly. Will be lower than the coverage noted in readMetrics.
* percentNs - the percentage of Ns in the assembly, from 0 to 1

