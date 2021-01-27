# SYNOPSIS

`sn_helloWorld.py`

Python example plugin for helping developers.
Because this is an example, please see the Algorithm section below for more details.

# Software requirements

None

# Algorithm

This script is an example plugin for the basics.
The steps are:

1. Argument parsing with `argparse` to conform to the required flags and positional parameter of every plugin
2. In `main()`, check on whether there are some options where we print something and exit: `--help`, `--version`, or `--check-dependencies`.  This is done in `setup(args)`.
3. If we have not exited already, move onto `main(args)`
4. Make a temporary directory to avoid collision with other instances of this script. However, the user could have also supplied `--tempdir` and so the script has to accept that parameter if supplied.
5. Make the required subfolders just in case they are not there: `dir/SneakerNet/forEmail`
6. The example analysis
    1. `samples = readWriteSamples(args.dir)`: rewrite `SampleSheet.csv` to `dir + "/SneakerNet/forEmail/helloworld.py.tsv"`
    2. `readWriteFlags(args.dir, args)`: record all parameters to the output table at `dir + "/SneakerNet/forEmail/helloworld.py.tsv"`
7. Wrap up: write all properties to the central properties table. Because other plugins write to this, the plugin _appends_ and does not truncate.
    1. `writeProperties(args.dir, samples)`
    2. The central table for all plugins is at `dir + "/SneakerNet/properties.txt"`
    3. Two entries are placed into `properties.txt`: the version and the table path.

## Epilogue

1. `properties.txt` will be read by `sn_report.pl` which will email a formatted table in the plugin `emailWhoever.pl`.
2. `helloworld.py.tsv` will be included in the email because it is in the folder `SneakerNet/forEmail`. All files in that folder will be attached to the email when `emailWhoever.pl` runs.

# Outputs

## Table

This is a combined table of keys/values and samples/sampleCounter.
An example table is shown below, when the script was run with `--debug`.

|sample             | sampleCounter |
|-------------------|---------------|
|2010EL-1786        |1|
|FA1090             |2|
|LT2                |3|
|Philadelphia_CDC   |4|
|contaminated       |5|
|version            |False|
|check_dependencies |False|
|citation           |False|
|debug              |True|
|force              |False|
|tempdir            ||
|numcpus            |1|
|dir                |../t/M00123-18-001-test|

