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
Will only transfer however, if `route` is set in
[`samples.tsv`](/docs/SneakerNetInput.md#samplestsv)
.

If `--careful`, then it will:

1. Look for a lock file. If found, the script will die.
If not found, the script will place a lock file.
2. Transfer the files with `rsync` (this step happens
regardless of `--careful`).
3. Remove the remote lock file.

If `--force`, then ignore warnings.
Removes the remote lock file if found.

If `--force-transfer`, then transfer the reads despite
any routing information in the spreadsheet.

If `--debug`, no files will actually be transferred.

# Outputs

none

