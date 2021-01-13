# SYNOPSIS

`sn_kraken.pl`

Runs Kraken1 on a sample.

# Software requirements

* Kraken1
* Database formatted for Kraken1

# Algorithm

Runs Kraken1 on a sample.
It looks for raw sequence reads (fastq) as input.
If not found, runs on genome assemblies under
runDir/SneakerNet/assemblies/sampleName/something.fasta.

# Outputs

No table output.
Creates files under runDir/SneakerNet/kraken/sampleName:

* kraken.report - a summary of counts per taxon with the following fields. More information here: https://ccb.jhu.edu/software/kraken/MANUAL.html#sample-reports
  * Percentage of reads covered by the clade rooted at this taxon
  * Number of reads covered by the clade rooted at this taxon
  * Number of reads assigned directly to this taxon
  * A rank code, indicating (U)nclassified, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. All other ranks are simply '-'.
  * NCBI taxonomy ID
  * indented scientific name
* kraken.filtered.report - the same as above, but all hits < 1% are removed
* kraken.taxonomy - a tab-delimited file
  * first column is the number of reads
  * second column, and all subsequent columns are the taxonomy starting with root and going down as far as possible, e.g., scientific name.
* report.html - graphical output in Krona format

