# SYNOPSIS

`sn_genotype.pl`

A plugin to run genotyping on a sample

# Software requirements

KMA version 1.4 or above

In taxonProperties.conf, need to have a URL to a fasta formated database,
e.g., `https://bitbucket.org/genomicepidemiology/serotypefinder_db/raw/ada62c62a7fa74032448bb2273d1f7045c59fdda/H_type.fsa`.
This needs to be done per taxon, and there is no option for all taxa.

# Algorithm

This algorithm only deploys if there is an entry for the given
taxon under `taxonProperties.conf`.

1. Download to the databases directory under SneakerNet (typically, `SneakerNet/db/genotype`)
2. Index with `kma index` default options
3. Query with `kma`
4. Record contents from `*.res` files

# Outputs

Table with columns directly from `.res` files

* Template
* Score
* Expected
* Template\_length
* Template\_Identity
* Template\_Coverage
* Query\_Identity
* Query\_Coverage
* Depth
* q\_value
* p\_value

