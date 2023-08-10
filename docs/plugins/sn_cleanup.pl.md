# SYNOPSIS

`sn_cleanup.pl`

Removes files from a SneakerNet subdirectory, usually to save space.

# Software requirements

# Algorithm

Files from the following directories are removed using the perl `unlink()` and `rmdir()` subroutines

* `$dir/SneakerNet/assemblies/*/shovill` 
  * Note: actual fasta assembly files are retained at `$dir/SneakerNet/assemblies/*/*.fasta`.
* `$dir/SneakerNet/assemblies/*/prodigal`
  * Note: actual gbk annotations are retained at `$dir/SneakerNet/assemblies/*/*.gbk`.

# Outputs

A message to stderr will state how many files were removed and how much space was saved.

Will put the following into properties.txt:

* numRemoved
* sizeRemoved - size in gigabytes

