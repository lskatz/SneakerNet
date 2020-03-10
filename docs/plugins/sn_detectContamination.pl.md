# SYNOPSIS

`sn_detectContamination.pl`

Runs a kmer histogram to detect potential contamination.

# Software requirements

* Perl module `Bio::Kmer`

# Algorithm

Runs a kmer histogram. The null hypothesis is that there will be
a peak around 1x coverage due to errors in sequencing, and a
second peak around the expected genome coverage. If something
else were sequenced at the same time at a different coverage,
then there will be another peak at that coverage level.
This other peak will be indicative of potential contamination.

# Outputs

Table with columns

* File
* numPeaks
* finalDelta - the amplitude difference between one coverage and the next to determine what a valley is
* hist       - ascii representation
* firstPeak  - coverage for the first peak
* firstValley - coverage for the first valley
* secondPeak
* secondValley...

More columns are possible with additional valleys

