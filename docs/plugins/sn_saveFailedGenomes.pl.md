# SYNOPSIS

`sn_saveFailedGenomes.pl`

Sends semi-failed genome reads to a QC_Fails subfolder

# Software requirements

* rsync

# Algorithm

If a genome has >= 10x coverage but less than the required coverage,
sends the genome reads to a subfolder, e.g., /mnt/CalculationEngineReads.test/Campy/QC_Fails/

# Outputs

none

