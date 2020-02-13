#!/usr/bin/env python3
# Hello World example

import tempfile
#import click
import argparse
import sys
import os
import csv
import subprocess
import shutil

VERSION='1.2'
CITATION='Hello World example by Lee Katz'

def main(args):

  if not os.path.isdir(args.dir):
    print("ERROR: could not find directory "+args.dir)
    sys.exit(1)

  with tempfile.TemporaryDirectory() as tempdirname:
    # If not defined already, create a new tempdir
    if not args.tempdir:
      tempdir = tempdirname

    # If the user defined the new temporary directory, then
    # just accept it.
    else:
      if os.path.isdir(args.tempdir):
        print("ERROR: "+args.tempdir+" already exists!")
        sys.exit(1)
      else:
        os.mkdir(args.tempdir)
    

    # make the directory structure just in case it isn't there yet
    for subdir in args.dir + "/SneakerNet", args.dir + "/SneakerNet/forEmail":
      if not os.path.isdir(subdir):
        os.mkdir(subdir)

    # analysis goes here.
    # In the hello world analysis, this generates a table
    # with samples and with attributes that were used in
    # argparse.
    samples = readWriteSamples(args.dir) # truncate and write to the table
    readWriteFlags(args.dir, args)       # append to the table

    # Record that hello world was run under the SneakerNet
    # properties file under dir/SneakerNet/properties.txt
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

    tsvout.write("{}\t{}\n".format("sample","sampleCounter"))
    counter = 0
    for row in tsvin:
      counter += 1
      tsvout.write("{}\t{}\n".format(row[0],counter))

  return samplesout

def readWriteFlags(dir, args):
  tableout = dir + "/SneakerNet/forEmail/helloworld.py.tsv"
  with open(tableout,'ta') as tsvout:
    for arg in vars(args):
      tsvout.write("{}\t{}\n".format(arg, getattr(args,arg)))
  return ""

def setup(args):
  if args.citation:
    print(CITATION)
    sys.exit(0)
  if args.version:
    print(VERSION)
    sys.exit(0)
  if args.check_dependencies:
    for exe in ["python3", "cat"]:
      # See if the executable is in the path
      path = shutil.which(exe)
      if not path:
        print("Could not find {} in PATH".format(exe))
        exit(1)

      # See if the executable gives us versioning
      try:
        out = subprocess.check_output([exe, "--version"], stderr=subprocess.STDOUT)
        # Just look at the first line
        for line in out.splitlines():
          print(line.decode("utf-8").rstrip())
          break
      except subprocess.CalledProcessError as e:
        print("Could not find dependency {}: {}".format(exe,e))
        exit(1)
    sys.exit(0)

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-v", "--version",  default=False, action='store_true', help="show the version and exit")
  parser.add_argument("--check-dependencies",  default=False, action='store_true', help="check dependencies")
  parser.add_argument("-c", "--citation", default=False, action='store_true', help="show the citation and exit")
  parser.add_argument("-d", "--debug",    default=False, action='store_true', help="show debugging information")
  parser.add_argument("-f", "--force",    default=False, action='store_true', help="Force output")
  parser.add_argument("-t", "--tempdir",  default="", type=str, help="Define where the temporary directory is")
  parser.add_argument("-n", "--numcpus", nargs=1, default=1,  type=int, help="Number of cpus")
  parser.add_argument("dir", nargs='?', default="")

  args = parser.parse_args()

  # Take care of version, citation
  setup(args)

  if args.dir == "":
    parser.print_help()

  # Main program
  sys.exit(main(args))

