#!/usr/bin/env perl

use strict;
use warnings;
use Test::More tests => 1;

use FindBin qw/$RealBin/;

use lib "$RealBin/../lib/perl5";
use_ok 'SneakerNet';

1;
