# Author: Lee Katz <lkatz@cdc.gov>
# SneakerNet
 
PREFIX := $(PWD)
NUMCPUS := 1
SHELL   := /bin/bash

###################################

default: install 

help:
	@echo "Please see README.md for additional help"

all: install env

install: install-prerequisites

install-prerequisites: install-perlModules

install-perlModules:
	@echo "Installing Perl modules using cpanminus"
	for package in Config::Simple File::Slurp Email::Stuffer; do \
	  perl scripts/cpanm --self-contained -L . $$package; \
		if [ $$? -gt 0 ]; then exit 1; fi; \
	done;
	@echo "Done with Perl modules"

