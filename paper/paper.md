---

title: 'SneakerNet: a modular quality assurance and quality check workflow for raw genomic and metagenomic read data'
authors:
- affiliation: 1
  name: Taylor Griswold
- affiliation: 1
  name: Curtis Kapsak
- affiliation: 1
  name: Jessica C. Chen
- affiliation: 1
  name: Grant Williams
- affiliation: 1
  name: Eshaw Vidyaprakash
- affiliation: "1, 2"
  name: Lee S. Katz
  orcid: 0000-0002-2533-9161
date: "28 Februrary, 2020"
bibliography: paper.bib
tags:
- QA/QC
affiliations:
- index: 1
  name: Enteric Diseases Laboratory Branch, Centers for Disease Control and Prevention,
    Atlanta, GA, USA
- index: 2
  name: Center for Food Safety, University of Georgia, Griffin, GA, USA

---

# Summary

Receiving a set of raw reads from whole genome sequencing or metagenomics sequencing has become commonplace and perhaps ubiquitous in bioinformatics.
However, there is a need to standardize the quality assurance and quality control process (QA/QC).
Therefore, we have created SneakerNet, a pipeline to standardize the QA/QC of a set of genomic or metagenomic reads.

There are very few standardized workflows for performing an initial quality control (QC) on those data.
For example, the Pandoo pipeline can be given a set of genomes to run analyses: species inference, 7-gene multilocus sequence typing (MLST), resistance gene profile, plasmid profile, virulence profile, and raw read QC [@Pandoo].
The Nullarbor pipeline is similar to Pandoo, but focused on public health datasets [@Nullarbor].
Another example is the ASA3P pipeline that runs raw read trimming, assembly, annotation, taxonomic classification, MLST, antibiotic resistance detection, virulence factor detection, reference mapping, and single nucleotide polymorphism (SNP) detection [@Schwengers654319].

# Implementation

## Plugin design

SneakerNet has a modular plugin design, where the main program calls each plugin in an ordered succession.
Each plugin is defined as being able to accept specific flagged and positional parameters, such that
the main program can call each plugin in a standardized way.
Workflows are thus defined as a specified order of plugins. One example order might be genome assembly, followed by MLST, followed by report generation.
At the time of this writing, 25 plugins were available.

## Plugin development

The plugin system has drastically lowered the activation energy needed to develop a new step in a
SneakerNet workflow. Documentation has been provided on how to develop a new plugin,
and 'Hello World' plugins have been published in three different languages: Perl, Python, and Bash.
Therefore SneakerNet has the ability to be a very collaborative project.

# References

