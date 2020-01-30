# SYNOPSIS

`sn_detectContamination-kraken.pl`

Attempts to guess a taxon for each sample and then lists
at most one major contaminant for the sample.

# Software requirements

* Kraken1
* Database formatted for Kraken1

# Algorithm

## species ID

Runs sample raw reads through Kraken1. Quantifies percentage of
reads that match each taxon.
Starting with species, if at least 25% of the reads correspond
with one species taxon, then records the majority species for
the sample.  If 25% are not captured at species, then it moves
onto genus, etc.

## contamination ID

For the rank that is identified (species, genus, etc),
if at least 5% of reads conflict with the species identified,
lists the major conflicting taxon.

# Outputs

Table with columns

* sample
* assumed taxon
* best-fitting taxon
* percent of reads that match the best-fitting taxon
* major conflicting taxon
* percent of reads supporting the major conflicting taxon

