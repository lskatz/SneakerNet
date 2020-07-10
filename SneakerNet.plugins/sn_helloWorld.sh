#!/bin/bash

set -e
set -u

# Set defaults
SCRIPT=$(basename $0)
MY_VERSION="2.1"
MY_CITATION="Hello World example in Bash by Lee Katz"

function usage(){
  echo "Usage: $0 [options] run/"
  echo "Prints all options to run/SneakerNet/properties.txt"
  echo "OPTIONS
  --help
  --numcpus 1
  --debug
  --tempdir ''
  --force
  --version
  --citation
  --check-dependencies
  "
}

OPTIONSUSED=""

# Set expected variables and defaults to avoid "-u errors"
# Booleans are ints in this example
HELP=0
NUMCPUS=1
DEBUG=0
TEMPDIR="" # TODO add mktemp
FORCE=0
VERSION=0
CITATION=0
CHECK_DEPENDENCIES=0
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --help|-h)
            HELP=1
            shift 1
            ;;

        --numcpus|-n)
            NUMCPUS=$2
            shift 2
            ;;

        --debug)
            DEBUG=1
            shift 1
            ;;

        --tempdir)
            TEMPDIR=$2
            shift 2
            ;;

        --force)
            FORCE=1
            shift 1
            ;;

        --version)
            VERSION=1
            shift 1
            ;;

        --citation)
            CITATION=1
            shift 1
            ;;

        --check-dependencies)
            CHECK_DEPENDENCIES=1
            shift 1
            ;;

        *)
            POSITIONAL+=("$1")
            shift
            ;;
    esac
done

# restore positional parameters if there are any
if [[ "${#POSITIONAL[@]}" -gt 0 ]]; then 
  set -- "${POSITIONAL[@]}"
fi

# print and exit for certain options
if [[ "$CITATION" == 1 ]]; then
  # Print the citation and exit
  echo "$MY_CITATION"
  exit 0;
fi
if [[ "$VERSION" == 1 ]]; then
  # Print the version and exit
  echo "$MY_VERSION"
  exit 0;
fi
if [[ "$CHECK_DEPENDENCIES" == 1 ]]; then
  # Check for dependencies: print results to stderr and execs checked to stdout

  # Check version of bash
  echo "bash"
  bash --version 2>&1 | grep -m 1 version 1>&2

  # Check version of basename
  echo "basename"
  basename --version 2>&1 | grep -m 1 basename 1>&2

  exit 0;
fi

# exit if there are no arguments
if [[ "$HELP" == 1 ]] || [[ $# -lt 1 ]]; then
  usage;
  exit 0;
fi

# The run number is the first positional argument
RUN=$1

echo "Run was given as $RUN"

# Write all options to a table
TABLE="$RUN/SneakerNet/forEmail/helloworld.sh.txt"
echo -e "foo\tbar"                                 > $TABLE
echo -e "run\t$RUN"                               >> $TABLE
echo -e "help\t$HELP"                             >> $TABLE
echo -e "numcpus\t$NUMCPUS"                       >> $TABLE
echo -e "debug\t$DEBUG"                           >> $TABLE
echo -e "tempdir\t$TEMPDIR"                       >> $TABLE
echo -e "force\t$FORCE"                           >> $TABLE
echo -e "VERSION\t$VERSION"                       >> $TABLE
echo -e "CITATION\t$CITATION"                     >> $TABLE
echo -e "check-dependencies\t$CHECK_DEPENDENCIES" >> $TABLE
echo "Wrote to table $TABLE"

# Record properties from this plugin into properties.txt
echo -e "$SCRIPT\ttable\t$TABLE"        >> $RUN/SneakerNet/properties.txt
echo -e "$SCRIPT\tversion\t$MY_VERSION" >> $RUN/SneakerNet/properties.txt

