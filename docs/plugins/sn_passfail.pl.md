# SYNOPSIS

`sn_passfail.pl`

Marks each sample as passing or failing

# Software requirements

# Algorithm

Reads upstream plugins for pass/fail criteria.
Currently these criteria are only coverage and
read quality, from the read metrics plugin.

# Outputs

Table with columns

* sample
* assembly - if all sub quality checks pass, then this passes. If all are unknown, then this is unknown. If any fail, then this fails.
  * number of Ns (threshold defined in `taxonProperties.conf`) or number of genes
  * longest contiguous contig - longest tract of genome between two Ns
* coverage - pass or fail
* quality  - pass or fail
* kraken - pass or fail 

1 indicates failure; 0 indicates pass; -1 indicates unknown (e.g., if Kraken was never run)

