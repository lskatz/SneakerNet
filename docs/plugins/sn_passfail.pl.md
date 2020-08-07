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
* coverage - pass or fail
* quality  - pass or fail
* kraken - pass or fail 

1 indicates failure; 0 indicates pass
