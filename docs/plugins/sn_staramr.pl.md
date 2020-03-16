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

Under taxonProperties.conf, add the following line under each taxon for pointfinder (SNP)
resistance detection. Most taxa are not supported right now. At the time of this documentation,
only these are supported: salmonella, campylobacter, enterococcus faecalis and enterococcus faecium.

    pointfinder=salmonella

# Algorithm

Runs StarAMR on genome assemblies.  Pointfinder is not currently enabled.

# Outputs

Table with columns

* Sample - sample name
* Assembly - which assembly was input
* Genotype 
* Predicted Phenotype

