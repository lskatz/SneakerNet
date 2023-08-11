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
  * Mash (before shovill 1.1.0)
  * KMC (in shovill 1.1.0)
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
* `gfa_connector`, part of the skesa package

# Algorithm

Assembles genomes with Skesa,
creates an assembly graph file with `gfa_connector`,
predicts genes with prodigal,
then produces assembly and gene prediction metrics.
Contigs with length < 500 are not considered in the assembly metrics.

# Outputs

Assemblies are found in `SneakerNet/assemblies/[samplename]/*.fasta`
GFA graph file found in `SneakerNet/assemblies/[samplename]/*.gfa`
Annotated assembly found in `SneakerNet/assemblies/[samplename]/*.gbk`

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
* effectiveCoverage - the number of mapped base pairs against the assembly. Will be lower than the coverage noted in readMetrics.

Starting with v2.5, four more output files are found in
`SneakerNet/assemblies/[samplename]/` but not included in the final report.
These files are used to make the aforementioned combined metrics table.

* predictionMetrics.tsv - gene prediction metrics such as CDS count
* assemblyMetrics.tsv - N50, genomeLength, etc
* depth.tsv.gz - `samtools depth` output: a three column file with contig, pos, depth of coverage
* effectiveCoverage.tsv - contains effective coverage

