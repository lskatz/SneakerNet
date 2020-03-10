# SYNOPSIS

`sn_staramr.pl`

Runs antimicrobial resistance analysis

# Software requirements

* StarAMR >= 0.2.2
* StarAMR database

## StarAMR database

Build the database into the correct folder under SneakerNet
like so

    staramr db build --dir SneakerNet/db/staramr

# Algorithm

Runs StarAMR on genome assemblies.  Pointfinder is not currently enabled.

# Outputs

Table with columns

* Sample - sample name
* Assembly - which assembly was input
* Genotype 
* Predicted Phenotype

