# SneakerNet

A pipeline for processing reads from a sequencing run. Currently only for Illumina-based runs.

# Synopsis

What do you do with a MiSeq run after it finishes? Are there basic analyses that you run 
such as counting how many reads you obtained? Checking for contamination? **SneakerNet performs
all these initial analyses with a nice report that is emailed to you at the end.** Additionally,
there is a plugins directory such that each installation is expandable.

# Installation

## As root

For these steps, log in as root.

    $ sudo su

### Create a dedicated user

There should be a special user that has permissions for sequencing runs.  Ours is 
'sequencingmaster'. All files made by sequencingmaster should have privileges for 
both group and user. In this way, if you want to give special privileges to other
users for sequencing runs, you will add them to the sequencingmaster group.
    
    $ useradd -m -s bash sequencermaster
    $ passwd sequencermaster

### Create a log file

    $ touch /var/log/SneakerNet.log
    $ chown sequencermaster.sequencermaster /var/log/SneakerNet.log
    $ chmod 664 /var/log/SneakerNet.log

### Create an 'inbox'

The inbox is some directory where sequencing runs will be deposited. In our lab, this
is a samba folder, but it could simply be a directory that user `cp` files to. The
inbox must have permissions like so:

    $ chmod +s /path/to/inbox/
    $ chown sequencermaster.sequencermaster /path/to/inbox/

## As 'sequencermaster'

For these steps, log in as sequencermaster.

    $ ssh sequencermaster@localhost

### Download the software

    $ mkdir ~/bin
    $ cd bin
    $ git clone git@git.biotech.cdc.gov:gzu2/sneakernet.git

### Set up the cron job (optional)

Log in again as sequencermaster and set up a cronjob to run once an hour. SneakerNet
will check for run directories once an hour. It will do nothing if there is no 
run directory present.

    $ crontab -e

Add this line to the crontab file, and then save it.

      1  *  *  *   *     ((date; /usr/bin/perl /home/sequencermaster/bin/SneakerNet.pl; echo;) >> /var/log/SneakerNet.log 2>&1)

### Configuration

You will need to edit some files for configuration before using SneakerNet.

    $ cd config

#### emails

List any emails here, one per line. These emails will be sent reports by default for each
SneakerNet run.

#### genomeSizes

This is a tab-delimited file with the following fields: Perl-style regular expression, genome size in bp, taxon.
If you want the read metrics script to be able to understand what taxon each file represents
by its name alone, list it here. Otherwise, leave this file blank.

#### settings

This file has certain key/values and should be left alone, unless you are a developer.

### Test it

If you have altered the configuration properly, then go ahead and test the software!

    $ ~/bin/SneakerNet.pl --test --now --numcpus 4

