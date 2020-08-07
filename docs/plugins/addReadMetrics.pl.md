# SYNOPSIS

`addReadMetrics.pl`

Adds raw read metrics for each fastq.gz file

# Software requirements

* CG-Pipeline

# Algorithm

Quantifies statistics like average read length.  Only measures 1% of the reads
and extrapolates for the whole file, and so some inconsistencies might be found
such as a single sample of a set that has a minimum read length of 38 when
all the others have 35 -- in this case, a minimum read length of 35 might be 
very rare but present.

# Outputs

Table with columns describing metrics of each reads file.
Reads are not interleaved and are quantified independently of either R1 or R2.

* File
* avgReadLength
* totalBases - i.e., total number of nucleotides
* minReadLength
* maxReadLength
* avgQuality
* numReads
* PE? - usually will be `no` because each split read is analyzed separately
* coverage - genome coverage calculated by `totalBases` divided by expected genome size found in `taxonProperties.conf`
* readScore - TODO
* medianFragmentLength - TODO
