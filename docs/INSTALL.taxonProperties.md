# Taxon Properties

Each organism has an entry in the file `config/taxonProperties.conf`.
This is useful for configuring each taxon.
For example, if your taxon is _Klebsiella_, you can add add a new entry 
`Klebsiella` into `taxonProperties.conf` and describe exactly the
thresholds you expect and all the properties of the taxon, such that
SneakerNet understands exactly how to treat this taxon.

To invoke the correct taxon, simply use it in `samples.tsv`.
See the taxon key in the info field in [SneakerNetInput.md](SneakerNetInput.md) for more information.

## Example

There is an example taxon (`SAMPLE_TAXON`) in `taxonProperties.conf` that you can follow.
Here is the entry as of version 0.9.9

    [SAMPLE_TAXON]
    # What coverage level you are setting as a threshold.
    # Below this coverage, the sample fails.
    coverage=30
    # If a taxon is not specified in samples.tsv, then
    # create a regex to see if we can guess what the sample is
    # based on the filename
    regex='^\d+q\-\d+'
    # What is the genome size in bp
    genomesize=5000000
    # What is the minimum average Q score you accept
    # before failing the sample
    quality=28
    # When the reads are transferred to a remote location
    # as specified in settings.conf, what subfolder do they
    # go to?
    dest_subfolder=Example
    # What subfolder do you have a wgMLST scheme in?
    wgMLST=Example
    # Which option to use for staramr for pointfinder
    pointfinder=example

## Specification

### Format

Each entry has to start with a header in brackets. In the example, the header is `[SAMPLE_TAXON]`.
The header corresponds to the entry in `samples.tsv`.

Comments are allowed per line.

Keys and values are separated by equals sign. Whitespace is allowed around the equals sign.

Values can be comma separated and the order is retained.

Further documentation on the format can be found in `[Config::Simple](https://metacpan.org/pod/Config::Simple)`.

### Specific properties

| Key            | Description |
|:---------------|:------------|
|coverage        |The minimum coverage required for a sample to pass. Without this coverage, a genome will be marked as "failed" which will have implications such as not transferring to the destination in the plugin `transferFilesToRemoteComputers.pl`.|
|regex           | If a taxon is not specified, SneakerNet will use this regex (Regular expression pattern) on the filename(s) to try to guess which taxon it is. Must be specified but you can give some nonsense value to ensure it does not match anything. |
|genomesize      | The estimated genome size in bp. Will be used to estimate genome coverage.|
|quality         | Minimum average phred score accepted before a sample will fail. Any individual file can cause the sample to fail, even if all others pass.|
|dest_subfolder  | Subfolder that raw reads will be transferred to in the plugin `transferFilesToRemoteComputers.pl`. See `transfer_destination_string` in `settings.conf` for more information. |
|wgMLST          | A subfolder in `SneakerNet/db/wgMLST` that contains a chewBBACA-formatted MLST scheme.|
|pointfinder     | A value to give to staramr for the pointfinder option.|


