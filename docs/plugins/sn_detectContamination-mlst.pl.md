# SYNOPSIS

`sn_detectContamination-mlst.pl`

Detects contamination using 7-gene MLST.

# Software requirements

* ColorID
* _mlst_

# Algorithm

Uses Torsten _mlst_ for the database. Maps raw reads
against the MLST database using ColorID. Because the
seven genes are housekeeping genes, then a non-seven
allele answer for any isolate indicates possible 
contamination.

# Outputs

Table with columns `mlst-contamination-detection.tsv`

* sample
* Scheme - scheme that this script estimates to be the correct scheme
* NumLociFound
* questionableLoci

