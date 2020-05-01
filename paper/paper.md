---
title: 'SneakerNet: A modular quality assurance and quality check workflow for primary genomic and metagenomic read data'
authors:
- affiliation: 1
  name: Taylor Griswold
- affiliation: "1, 2"
  name: Curtis Kapsak
  orcid: 0000-0002-8735-1190
- affiliation: 1
  name: Jessica C. Chen
- affiliation: 3
  name: Henk C. den Bakker
  orcid: 0000-0002-4086-1580
- affiliation: 1
  name: Grant Williams
- affiliation: "2, 4"
  name: Alyssa Kelley
- affiliation: "1, 5"
  name: Eshaw Vidyaprakash
- affiliation: "1, 3"
  name: Lee S. Katz
  orcid: 0000-0002-2533-9161
date: "29 April, 2020"
bibliography: paper.bib
tags:
- QA/QC
affiliations:
- index: 1
  name: Enteric Diseases Laboratory Branch (EDLB), Centers for Disease Control and Prevention,
    Atlanta, GA, USA
- index: 2
  name: Weems Design Studio, Inc., Suwanee, GA, USA
- index: 3
  name: Center for Food Safety, University of Georgia, Griffin, GA, USA
- index: 4
  name: Waterborne Disease Prevention Branch (WDPB), Centers for Disease Control and Prevention,
    Atlanta, GA, USA
- index: 5
  name: IHRC, Atlanta, GA, USA
  
---

# Summary

Receiving a set of primary data from whole genome sequencing or metagenomics sequencing has become commonplace and perhaps ubiquitous in bioinformatics.
However, there is a need to standardize the quality assurance and quality control process (QA/QC).

There are very few published workflows for performing an analysis on primary sequence data that span the breadth of initial standardized QA/QC (e.g., sequence yields, contamination checks, and subtyping).
For example, the Pandoo pipeline can be given a set of genomes to run analyses: species inference, 7-gene multilocus sequence typing (MLST), resistance gene profile, plasmid profile, virulence profile, and raw read QC [@Pandoo].
The Nullarbor pipeline is similar to Pandoo, but focused on public health datasets [@Nullarbor].
Another example is the ASA3P pipeline that runs raw read trimming, assembly, annotation, taxonomic classification, MLST, antibiotic resistance detection, virulence factor detection, reference mapping, and single nucleotide polymorphism (SNP) detection [@Schwengers654319].
However, no existing "broad stroke" QA/QC pipelines seem to be focused on a plugins-based architecture for batches of unrelated bacterial sequences or for batches of bacteria from different species.
To that end, we have created SneakerNet.
The major design principles are centered around the ability to collaboratively design plugins.
With the plugins architecture, SneakerNet can dynamically change for current and future needs
with input from the bioinformatics and public health community.

# Implementation

## Plugin design

SneakerNet has a modular plugin design, where the main program calls each plugin in an ordered succession.
Each plugin, in turn, reads a set of genomes as input.
Each plugin accepts specific flagged parameters, such that
the main program can call each plugin in a standardized way.
Workflows are thus defined as a specified order of plugins. 
An example workflow order might be genome assembly, followed by MLST, and finalized with report generation.
At the time of this writing, 25 plugins are available.
These plugins are listed in the documentation in a summary table,
and each plugin has its own documentation page.

## Plugin development

The plugin system has drastically lowered the activation energy needed to develop a new step in a
SneakerNet workflow. Documentation has been provided on how to develop a new plugin,
and 'Hello World' plugins have been published in three different languages: Perl, Python, and Bash.
Because plugins are not tied to any specific language, SneakerNet collaborators do not have to be bound by any specific language.

## Configuration

SneakerNet is highly configurable as described in the installation documentation.
There are many configurations.
We would like to highlight some ways that SneakerNet can be configured.

For some genera, SneakerNet comes packaged with some recommended configurations (e.g., _Salmonella_ or _Legionella_),
and an example genus with all options commented.
These options include the minimum coverage needed for a sample to pass QC
and even some detailed options to help customize a taxon for a particular plugin such as the antimicrobial resistance plugin.
Therefore, a user could easily add a taxon to customize the workflow for his or her instance of SneakerNet.
In fact, SneakerNet has been recently configured to accommodate the protist _Cryptosporidium_ successfully with input from the CDC Waterborne laboratory.

Users can also customize the workflow.
SneakerNet comes packaged with a default workflow which specifies the order of plugins that are run.
However, if a certain analysis is not needed, e.g., 7-gene MLST, then it can be removed from the configuration.
Likewise, if a new plugin is needed, it can be added into the workflow.

# Acknowledgements

The authors would like to thank the following SneakerNet users in the Enteric Diseases Laboratory Branch for useful feedback: Heather Carleton, Katie Dillon, Blake Dinsmore, Yang Gao, Jessica Halpin, Jasmine Hensley, Monica Im, Justin Kim, Charlotte Lane, Ana Lauer, Rebecca Lindsey, Angela Poates, Zachary Rigney, Katie Roache, Ashley Sabol, Peyton Smith, Cheryl Tarr, Jenny Truong, Maryann Turnsek; Additionally, SneakerNet users from other branches: Shatavia Morrison. The authors would also like to thank Aaron Petkau for useful feedback for staramr.

# References

