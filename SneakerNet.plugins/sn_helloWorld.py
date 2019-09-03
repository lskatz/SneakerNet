#!/usr/bin/env python3
# Hello World example

import tempfile
import click
import sys
import os
import csv

VERSION='1.0'
CITATION='Hello World example by Lee Katz'

@click.command()
@click.option("-v", "--version",  is_flag=True, default=False, help="show the version and exit")
@click.option("-c", "--citation", is_flag=True, default=False, help="show the citation and exit")
@click.option("-d", "--debug",    is_flag=True, default=False, help="show debugging information")
@click.option("-f", "--force",    is_flag=True, default=False, help="Force output")
@click.option("-t", "--tempdir", default="", type=str, help="Define where the temporary directory is")
@click.option("-n", "--numcpus", nargs=1, default=1,  type=int, help="Number of cpus")
@click.argument("dir", nargs=1, default="")

def main(dir,version, citation, force, debug, tempdir, numcpus):
  # Take care of version, citation
  setup(version, citation)

  if not os.path.isdir(dir):
    print("ERROR: could not find directory "+dir)
    sys.exit(1)

  with tempfile.TemporaryDirectory() as tempdirname:
    # If not defined already, create a new tempdir
    if not tempdir:
      tempdir = tempdirname

    # If the user defined the new temporary directory, then
    # just accept it.
    else:
      if os.path.isdir(tempdir):
        print("ERROR: "+tempdir+" already exists!")
        sys.exit(1)
      else:
        os.mkdir(tempdir)
    

    # make the directory structure just in case it isn't there yet
    for subdir in dir + "/SneakerNet", dir + "/SneakerNet/forEmail":
      if not os.path.isdir(subdir):
        os.mkdir(subdir)

    # analysis
    samples = readWriteSamples(dir)
    writeProperties(dir, samples)

    return 0

  return 1

def writeProperties(dir, samples):
  propertiesFile = dir+"/SneakerNet/properties.txt"
  with open(propertiesFile, 'ta+') as fileout:
    fileout.write("{}\t{}\t{}\n".format(
      "sn_helloWorld.py", "version", VERSION
    ))
    fileout.write("{}\t{}\t{}\n".format(
      "sn_helloWorld.py", "table", samples
    ))

def readWriteSamples(dir):
  samplesout = dir + "/SneakerNet/forEmail/helloworld.py.tsv"
  samples = dir + "/samples.tsv"
  with open(samples, 'rt') as tsvin, open(samplesout,'wt') as tsvout:
    tsvin = csv.reader(tsvin, delimiter='\t')
    #tsvout= csv.writer(tsvout,delimiter='\t')

    tsvout.write("{}\t{}\n".format("sample","sampleCounter"))
    counter = 0
    for row in tsvin:
      counter += 1
      tsvout.write("{}\t{}\n".format(row[0],counter))
      #tsvout.writerow(row[0])

  return samplesout

def setup(version, citation):
  if citation:
    print(CITATION)
    sys.exit(0)
  if version:
    print(VERSION)
    sys.exit(0)

if __name__ == "__main__":
  # Main program
  sys.exit(main())

