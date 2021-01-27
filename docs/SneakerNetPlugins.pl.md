# SYNOPSIS

`SneakerNetPlugins.pl`

Runs all plugins for a SneakerNet-formatted run

# Software requirements

See: [dependencies](/docs/INSTALL.md#dependencies)

# Algorithm

Reads the workflow from `snok.txt` in an input folder.
If one is not provided, will run the default workflow.
Workflows are described in [PLUGINS.md](/docs/PLUGINS.md).

Each workflow's plugins and the order of plugins
are determined by `plugins.conf`, documented in [PLUGINS.md](/docs/PLUGINS.md).

# Outputs

Virtually all workflows create:

* [`report.html`](/docs/plugins/sn_report.pl.md)
* an [email](/docs/plugins/emailWhoever.pl.md) with the report

This script also produces a variety of outputs described
in their own documentation [SneakerNetOutput.md](/docs/SneakerNetOutput.md).

