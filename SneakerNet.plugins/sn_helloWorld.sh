#!/bin/bash

function usage(){
  echo "Usage: $0 [options] run/"
  echo "options"
  echo "  -h, -d, -t, -f, -v, -c"
  echo "  -n 1"
}

# Set defaults
NUMCPUS=1
SCRIPT=$(basename $0)
VERSION=1

OPTIONSUSED=""

while getopts "hdtfvcn:" opt; do
  case ${opt} in
    h )
      usage;
      exit 0;
      ;;
    d )
      echo "You provided '-d'"
      OPTIONSUSED+=" d"
      ;;
    t )
      echo "You provided '-t'"
      OPTIONSUSED+=" t"
      ;;
    f )
      echo "You provided '-f'"
      OPTIONSUSED+=" f"
      ;;
    v )
      echo "You provided '-v'"
      OPTIONSUSED+=" v"
      ;;
    c )
      echo "You provided '-c'"
      OPTIONSUSED+=" c"
      ;;
    n )
      NUMCPUS=$OPTARG
      OPTIONSUSED+=" n=$NUMCPUS"
      ;;
    \? )
      echo "Invalid option specified!"
      exit 1
    ;;
  esac
done
shift $((OPTIND -1))

RUN=$1

echo $RUN

if [ -z $RUN ]; then
  usage
  exit 1
fi

TABLE="$RUN/SneakerNet/forEmail/helloworld.sh.tsv"
echo -e "foo\tbar" > $TABLE

echo -e "$0\ttable\t$TABLE"     >> $RUN/SneakerNet/properties.txt
echo -e "$0\tversion\t$VERSION" >> $RUN/SneakerNet/properties.txt
echo -e "$0\toptions\t$OPTIONSUSED" >> $RUN/SneakerNet/properties.txt

