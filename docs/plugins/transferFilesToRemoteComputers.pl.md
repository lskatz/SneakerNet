# SYNOPSIS

`transferFilesToRemoteComputers.pl`

Transfers files to a remote computer as specified in
config/settings.conf

# Software requirements

* rsync

# Algorithm

Will transfer each set of raw reads to the remote
computer.  Will not transfer if transfer was not
requested in samples.tsv. Also will not transfer
if did not pass in `sn_passfail.pl`.

# Outputs

none

