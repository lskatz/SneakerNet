---

name: Plugin contribution
about: Contribute a new plugin
title: "[PLUGIN CONTRIBUTION]"
labels: ''
assignees: ''

---

DOCUMENTATION
============

How to develop a plugin - (PLUGINSDEV.md)[https://github.com/lskatz/SneakerNet/blob/master/docs/PLUGINSDEV.md]

What is a plugin - (PLUGINS.md)[https://github.com/lskatz/SneakerNet/blob/master/docs/PLUGINSDEV.md]

CHECKLIST
=======

* [ ] First positional argument is the SneakerNet-formatted directory
* [ ] All flags are accepted in PLUGINS.md
* Soft-coded variables
  * [ ] Added to a conf file in `settings/`
  * [ ] Documented in the plugin usage
* [ ] Sample information is read from `samples.tsv`
* Outputs
  * [ ] File outputs are in `$rundirectory/SneakerNet/pluginName`
  * [ ] Any files to be attached to the SneakerNet report are copied into `$rundirectory/SneakerNet/forEmail`
* properties recorded after the plugin is run each time, in `$rundirectory/SneakerNet/properties.tsv`.
  * [ ] version is added to properties.tsv
  * [ ] any table paths are added to properties.tsv
