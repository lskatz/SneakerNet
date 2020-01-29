#!/usr/bin/env python3
# Hello World example

import tempfile
#import click
import argparse
import sys
import os
import csv

VERSION='1.1'
CITATION='Hello World example by Lee Katz'

def main(args):
  # Take care of version, citation
  setup(args.version, args.citation)

  if not os.path.isdir(args.dir):
    print("ERROR: could not find directory "+args.dir)
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
    for subdir in args.dir + "/SneakerNet", args.dir + "/SneakerNet/forEmail":
      if not os.path.isdir(subdir):
        os.mkdir(subdir)

    # analysis
    samples = readWriteSamples(args.dir)
    writeProperties(args.dir, samples)

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
  parser = argparse.ArgumentParser()
  parser.add_argument("-v", "--version",  default=False, action='store_true', help="show the version and exit")
  parser.add_argument("-c", "--citation", default=False, action='store_true', help="show the citation and exit")
  parser.add_argument("-d", "--debug",    default=False, action='store_true', help="show debugging information")
  parser.add_argument("-f", "--force",    default=False, action='store_true', help="Force output")
  parser.add_argument("-t", "--tempdir",  default="", type=str, help="Define where the temporary directory is")
  parser.add_argument("-n", "--numcpus", nargs=1, default=1,  type=int, help="Number of cpus")
  parser.add_argument("dir", nargs=1, default="")

  args = parser.parse_args()

  # Main program
  sys.exit(main(args))

