# SYNOPSIS

`sn_report.pl`

Creates an HTML report for a SneakerNet run. Usually,
this is the second-to-last plugin to run, right before
the email plugin.

# Software requirements

# Algorithm

Reads all properties in `SneakerNet/properties.txt`, which
has three columns: plugin, key, value. If a duplicate
plugin/key combination is found, then the latest value
is used.

If the key is `table`, then the corresponding table is 
converted to HTML and included in the report. All other
key/value combinations are included in the HTML report,
under their corresponding plugins.

# Outputs

HTML report

