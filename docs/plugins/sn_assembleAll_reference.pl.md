# SYNOPSIS

`sn_assembleAll_reference.pl.md`

Runs a reference assembly against a reference genome assembly
and produces assembly metrics.

# Software requirements

* Bowtie2
* samtools
* bcftools
* trimmomatic
* seqtk

# Algorithm

## Database building

* Uses `reference_fasta_id` to download a fasta and genbank file from NCBI
  * Found in taxonProperties.conf
  * Can use comma separated list, e.g., `reference_fasta_id=MH185784,MH185777,MH185772,MN367319,MN367321,MN367323,MN367326,MH430075`
* Downloads to SneakerNet installation, into a subfolder `db/fasta`
* Indexes the fasta file with `samtools` and `bowtie2-build`
* Gene annotations are not indexed and can be read on the fly

## Assembly

* Downloads the reference genome into your installation directory
* Trims adapters with trimmomatic
  * primers are found in taxonProperties.conf
  * e.g., primers_bed_url=https://raw.githubusercontent.com/artic-network/primer-schemes/master/nCoV-2019/V3/nCoV-2019.primer.bed
  * If no primers, then will not trim and will produce a warning
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

Some intermediate files are in `SneakerNet/assemblies/[samplename]/consensus
* depth.tsv - output of `samtools depth`
* out.vcf.gz.masked.vcf.gz and related `csi` file
* sorted.bam and related `bai` file

A graph of all depths of coverage in `SneakerNet/forEmail/depth.png`

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



