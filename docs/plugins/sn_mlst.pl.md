# SYNOPSIS

`sn_mlst`

Runs 7-gene MLST on each sample

# Software requirements

* Torsten _mlst_

# Algorithm

For each genome assembly, runs _mlst_ and reports
an individual file in each sample subdirectory.

# Outputs

Table with columns

* sample
* mlst_scheme
* 7-gene_ST
* locus1
* locus2
* locus3
* locus4
* locus5
* locus6
* locus7

