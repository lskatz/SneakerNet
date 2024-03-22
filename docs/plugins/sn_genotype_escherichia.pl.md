# SYNOPSIS

`sn_genotype_escherichia.pl`

A plugin to interpret `sn_genotype.pl` raw results for just Escherichia

# Software requirements

This depends on `sn_genotype.pl` having been run on the serotypefinder database
found at <https://bitbucket.org/genomicepidemiology/serotypefinder_db/raw/ada62c62a7fa74032448bb2273d1f7045c59fdda/H_type.fsa>.

# Algorithm

For both O and H results, reports the genotype with the highest score.

# Outputs

Table with columns directly from `.res` files

* sample
* O
* H

