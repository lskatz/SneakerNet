---
title: 'SneakerNet: a modular quality assurance and quality check workflow for raw genomic and metagenomic read data'
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
- affiliation: "1, 4"
  name: Alyssa Kelley
- affiliation: "1, 5"
  name: Eshaw Vidyaprakash
- affiliation: "1, 3"
  name: Lee S. Katz
  orcid: 0000-0002-2533-9161
date: "10 April, 2020"
bibliography: paper.bib
tags:
- QA/QC
affiliations:
- index: 1
  name: Enteric Diseases Laboratory Branch (EDLB), Centers for Disease Control and Prevention,
    Atlanta, GA, USA
- index: 2
  name: Weems Design Studio, Inc., Suwanee, GA, USA
- index: 5
  name: IHRC, Atlanta, GA, USA
- index: 3
  name: Center for Food Safety, University of Georgia, Griffin, GA, USA
- index: 4
  name: Waterborne Disease Prevention Branch (WDPB), Centers for Disease Control and Prevention,
    Atlanta, GA, USA

---

# Summary

Receiving a set of raw reads from whole genome sequencing or metagenomics sequencing has become commonplace and perhaps ubiquitous in bioinformatics.
However, there is a need to standardize the quality assurance and quality control process (QA/QC).
Therefore, we have created SneakerNet, a pipeline to standardize the QA/QC of a set of genomic or metagenomic reads.

There are very few standardized workflows for performing an initial QA/QC on those data.
For example, the Pandoo pipeline can be given a set of genomes to run analyses: species inference, 7-gene multilocus sequence typing (MLST), resistance gene profile, plasmid profile, virulence profile, and raw read QC [@Pandoo].
The Nullarbor pipeline is similar to Pandoo, but focused on public health datasets [@Nullarbor].
Another example is the ASA3P pipeline that runs raw read trimming, assembly, annotation, taxonomic classification, MLST, antibiotic resistance detection, virulence factor detection, reference mapping, and single nucleotide polymorphism (SNP) detection [@Schwengers654319].

# Implementation

## Plugin design

SneakerNet has a modular plugin design, where the main program calls each plugin in an ordered succession.
In turn, each plugin reads a set of genomes as input.
Each plugin accepts specific flagged parameters, such that
the main program can call each plugin in a standardized way.
Workflows are thus defined as a specified order of plugins. One example order might be genome assembly, followed by MLST, followed by report generation.
At the time of this writing, 25 plugins are available.

## Plugin development

The plugin system has drastically lowered the activation energy needed to develop a new step in a
SneakerNet workflow. Documentation has been provided on how to develop a new plugin,
and 'Hello World' plugins have been published in three different languages: Perl, Python, and Bash.
Because plugins are not tied to any specific language, SneakerNet has the ability to be a very collaborative project.

# Acknowledgements

SneakerNet users for useful feedback in the Enteric Diseases Laboratory Branch: Heather Carleton, Katie Dillon, Blake Dinsmore, Yang Gao, Jasmine Hensley, Monica Im, Justin Kim, Charlotte Lane, Rebecca Lindsey, Angela Poates, Zachary Rigney, Katie Roache, Ashley Sabol, Peyton Smith, Cheryl Tarr, Jenny Truong, Maryann Turnsek; SneakerNet users from other branches: Shatavia Morrison; Useful feedback for staramr: Aaron Petkau

# References

