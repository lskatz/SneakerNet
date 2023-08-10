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

Properties are recorded into `SneakerNet/properties.txt` from other plugins.
With the perl plugins, these are recorded using `SneakerNet::recordProperties()` which is documented via `perldoc lib/perl5/SneakerNet.pm`.

## Warnings

If the plugin key is `warnings`, then it will display any warnings at the top of the report.

## Errors

If the plugin key is `errors`, then it will display any errors at the top of the report.

## tables

If the value has a suffix `.csv` or `.tsv`, then the corresponding table is 
converted to HTML and included in the report. All other
key/value combinations are included in the HTML report,
under their corresponding plugins.

## images

_For this plugin version >= 2.7_

If the value has a suffix `.png` or `.gif`, then the corresponding
file path will be converted to base64 and embedded in the HTML report.

# Outputs

HTML report

## emoji output

This script itself outputs a table and is therefore converted to a table in the HTML report.
The table has columns:

* sample - sample name
* emoji - reflective of the score.  Happiest emojis reflect 100%.
  * As of v0.15, emoticons range from &#128515; (best), &#129320;, &#128556;, and &#128561; (worst).
* score - a percentage, starting from 100.  Each item under the failure_code column subtracts an equal percentage from 100%.  These possible failures are shown as columns in the [passfail plugin](sn_passfail.pl.md).  If there are three possible items, then each penalty is 33%.  By default in SneakerNet version 0.10, there are three possible items: coverage, quality, and kraken.
* qual - quality scores of R1 and R2, separated by space.
* cov - genome coverages of R1 and R2, separated by space.
* taxon - the calculated taxon. The calculated taxon is, in priority order: the taxon listed on the sample spreadsheet, guessed from Kraken, pattern matching on the filename, and lastly "UNKNOWN". For example, if the taxon is specified on the sample spreadsheet, it will be used and not overwritten.

