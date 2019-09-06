# SYNOPSIS

Attempts to guess a taxon for each sample

# Software requirements

* Kraken1
* Database formatted for Kraken1

# Algorithm

Runs sample raw reads through Kraken1. Quantifies percentage of
reads that match each taxon.
Compares the assigned taxon vs the assumed taxon.

# Outputs

Table with columns

* sample
* assumed taxon
* best-fitting taxon
* percent of reads that match the best-fitting taxon

