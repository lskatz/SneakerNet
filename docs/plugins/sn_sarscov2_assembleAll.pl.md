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

## Database building

* Uses `reference_fasta_id` to download a fasta and genbank file from NCBI
* Downloads to SneakerNet installation, into a subfolder `db/fasta`
* Indexes the fasta file with `samtools` and `bowtie2-build`
* Gene annotations are not indexed and can be read on the fly

## Assembly

* Downloads the reference genome into your installation directory
* Trims adapters with cutadapt
* Maps reads against genome with bowtie2
* Calls SNPs with samtools
* Filters SNPs
* Consensus fasta from SNPs and the reference genome

## Annotation

* copies coordinates of reference genome genes to assembly
* Counts intact coding sequences designated by the CDS tag in the genbank file
* Counts intact coding sequences in the new assembly

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
* altCdsCount - number of CDS in new assembly
* refCdsCount - number of CDS in reference assembly
* expectedCdsPercentage - altCdsCount / refCdsCount
* longestContiguous - the longest contiguous part of the genome between any two Ns



