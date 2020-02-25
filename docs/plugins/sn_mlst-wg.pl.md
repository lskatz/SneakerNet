# SYNOPSIS

`sn_mlst-wg.pl`

Runs a whole-genome MLST-style analysis. Scheme can be
a whole genome MLST scheme, core genome MLST scheme,
or anything smaller.

# Software requirements

* chewBBACA
* BLAST+
* Database

## The database

The database must be listed under taxonProperties.conf. At least
one example is provided at this time for Salmonella.

The raw database must be downloaded. One method is with https://github.com/Public-Health-Bioinformatics/pubmlst_client

For example,

    pubmlst_download --scheme_name salmonella --scheme_id 4 --outdir ./salmonella_enterobase.raw

The database must be formatted according to chewBBACA instructions
at https://github.com/B-UMMI/chewBBACA/wiki/1.-Schema-Creation#14-using-an-external-schema

    chewBBACA.py PrepExternalSchema -i salmonella_enterobase.raw -o wgMLST/salmonella.enterobase.chewBBACA --cpu 4
    rm -rf salmonella_enterobase.raw

Next, the database must be located under your SneakerNet installation
directory, under `SneakerNet/db/wgMLST/something` where _something_
is the name of the scheme you provided in taxonProperties.conf.

# Algorithm

Runs chewBBACA which is a very very smart BLAST. For more information,
please see https://github.com/B-UMMI/chewBBACA

# Outputs

Summary table

For more information, please see https://github.com/B-UMMI/chewBBACA/wiki/2.-Allele-Calling#allele-call-statistics-output-results_statisticstxt

* Genome (input fasta file)
* db (the database listed under taxonProperties.conf)
* EXC alleles which have exact matches (100% DNA identity) with previously identified alleles
* INF inferred new alleles using Prodigal CDS predictions
* LNF loci not found.
* PLOT possible loci on the tip of the query genome contigs.
* NIPH non-informative paralogous hit
* NIPHEM similar to NIPH classification (NIPH with exact match), but specifically referring to exact matches
* ALM alleles 20% larger than length mode of the distribution of the matched loci
* ASM similar to ALM but for alleles 20% smaller than length mode distribution of the matched loci 
