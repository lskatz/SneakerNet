# SYNOPSIS

`sn_report.pl`

Creates an HTML report for a SneakerNet run. Usually,
this is the second-to-last plugin to run, right before
the email plugin.

# Software requirements

# Algorithm

## key value pairs

Reads all properties in `SneakerNet/properties.txt`, which
has three columns: plugin, key, value. If a duplicate
plugin/key combination is found, then the latest value
is used.

## tables

If the key is `table`, then the corresponding table is 
converted to HTML and included in the report. All other
key/value combinations are included in the HTML report,
under their corresponding plugins.

# Outputs

HTML report

## emoji output

This script itself outputs a table and is therefore converted to a table in the HTML report.
The table has columns:

* sample - sample name
* emoji - reflective of the score.  Happiest emojis reflect 100%.
As of v0.14, emoticons range from &#128515; (best), &#129320;, &#128556;, and &#128561; (worst).
* score - a percentage, starting from 100.  Each item under the failure_code column subtracts an equal percentage from 100%.  These possible failures are shown as columns in the [passfail plugin](sn_passfail.pl.md).  If there are three possible items, then each penalty is 33%.  By default in SneakerNet version 0.10, there are three possible items: coverage, quality, and kraken.
* qual - quality scores of R1 and R2, separated by space.
* cov - genome coverages of R1 and R2, separated by space.
* taxon - the calculated taxon. The calculated taxon is, in priority order: the taxon listed on the sample spreadsheet, guessed from Kraken, pattern matching on the filename, and lastly "UNKNOWN". For example, if the taxon is specified on the sample spreadsheet, it will be used and not overwritten.

