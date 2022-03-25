# SYNOPSIS

`sn_saveFailedGenomes.pl`

Sends semi-failed genome reads to a QC_Fails subfolder

# Software requirements

* rsync

# Algorithm

If a genome has >= 10x coverage but less than the required coverage,
sends the genome reads to a subfolder called `QC_Fails`. 
For example, /mnt/CalculationEngineReads.test/Campy/QC_Fails/

If any genomes are "saved," then a warning will be emitted for the SneakerNet report.

# Outputs

A table `qc_fails.tsv` with two columns

* sample - sample name
* filename - file name that was sent to the QC_Fails subfolder

