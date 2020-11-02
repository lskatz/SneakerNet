#!/usr/bin/env perl
use strict;

#Input is a blast output that contains:
#
#qseqid sseqid pident length qcovhsp mismatch gapopen 
#qstart qend sstart send qlen slen qseq
#
#in that order
#The sseqid needs to be in the format of: Genus_Species_SubtypeFamily_string
# Ex: C_hominis_IfA12_string

my $input = @ARGV[0]; #tab delim blast output

my $combine;

#For the additional character at the end of
#subtypes IIc, chipmunk, and skunk
my $endChar = "";

#For now, we are only going to get the best match so we only want
#to grab the first line but later we might want to include all
#matches so keep the loop in and just limit the iterations
my $loopNum = 0;
 

open FILE, $input or die $!;
while (my $line = <FILE>) {

	$loopNum++;
	$line =~ s/\r\n/\n/; #remove any Windows characters
	chomp $line;
	

	#qseqid sseqid pident length qcovhsp mismatch gapopen qstart qend sstart send qlen slen qseq
	my ($qseqid, $sseqid, $pident, $length, $qcovhsp, $mismatch, $gapopen, $qstart, $qend, $sstart, $send, $qlen, $slen, $qseq) = split("\t", $line);

	#Check if sequence needs to be reverse complemented
	#otherwise will not properly count repeats
	if ($sstart > $send) {
		my $revcomp = reverse $qseq;
		$revcomp =~ tr/ATGCatgc/TACGtacg/;
		$qseq = $revcomp;
	}
	

	#parse subtype family from sseqid
	#sometimes the subtype family will be by itself and sometimes it
	# will be followed by repeat details always starting with "A"
	my ($genus, $species, $subFamily, $rest) = split("_", $sseqid, 4);
	
	my $subtype;
	if (index($subFamily, "A") != -1) {

		#parse based on the first occurence of A
		my ($sub, $rest) = split("A", $subFamily, 2);
		$subtype = $sub;
	} else {

		$subtype = $subFamily;
	}


	#Find percentage of Hit matched so we know how much of the
	#gp60 was found
	my $hPercent;
	if ($send > $sstart) {

		$hPercent = (($send-$sstart+1)/$slen)*100;
	} else {

		$hPercent = (($sstart-$send+1)/$slen)*100;
	}

	#Ubiquitum reference says no repeats
	if ($species eq "ubiquitum") {

		print "$input\t$qseqid\t$sseqid\t$subtype\t$length\t$pident\t$hPercent\n";
	} 

	#IIc, chipmunk, and skunk genotype have an additional
	#char at the end of their subtypes
	#>C_parvum_IIcA5G3a_AF164491.1
	#>Skunk_genotype_XVIa14a_KP099095.1
	#>Chipmunk_genotypeI_XIVaA18G2T1a_KP099082.1
	#>C_viatorum_XVaA3a_KP115936.1
	elsif ((index($sseqid, "_IIc") != -1) | (index(lc($sseqid), "skunk") != -1) | (index(lc($sseqid), "chipmunk") != -1) | (index(lc($sseqid), "viatorum") != -1)) {

		my ($triRepeat, $motif, $repeatSeq) = FindRepeatNumber($qseq);
		my ($rRepeats) = FindRRepeats($subtype, $qseq);

		$endChar = substr($subFamily, -1);

		if ($endChar =~ /^[a-zA-Z]/) {

			#print "TEST\n";##### Remove ###
			print "$input\t$qseqid\t$sseqid\t$subtype$triRepeat$rRepeats$endChar\t$length\t$pident\t$hPercent\n";
		} else {

			print "$input\t$qseqid\t$sseqid\t$subtype$triRepeat$rRepeats\t$length\t$pident\t$hPercent\n";
		}

	} else {

		my ($triRepeat, $motif, $repeatSeq) = FindRepeatNumber($qseq);
		my ($rRepeats) = FindRRepeats($subtype, $qseq);
		print "$input\t$qseqid\t$sseqid\t$subtype$triRepeat$rRepeats\t$length\t$pident\t$hPercent\n";
	}

	#Only grab the first (best) hit in the blast file for now
	#remove for testing new databases or for multiple samples
	#in one file
#########################################
	#if ($loopNum >= 1) {

	#	last;
	#}
}
close FILE;


## Subroutines ##

# Searches GP60 motifs and then finds
# the number of TCA, TCG, TCT repeats
# in the input sequence.
sub FindRepeatNumber {


	#@_[0] =~ m/(C+TGTTG[AG][GT][GA]G[CAT]|GA[GT][GA]G[CAT]|[CAT])((TC[AGT]){2,})(A)/;
	#@_[0] =~ m/(C+TGTTG[AG][GT][GA]G[CAT]|GA[GT][GA]G[CAT]|[CAT])((TC[AGT])(TC[AGT])+)(A)/;
	#@_[0] =~ m/(TC[AGT](TC[AGT])+)(A)/;
	@_[0] =~ m/([CT]C[AGT]([CT]C[AGT])+)(A)/;

	#my $repeats = $2;
	#my $motif = $1;
	my $repeats = $1;
	my $motif = "";

	#print "$input\t$repeats\n";

	my $ACount = ($repeats =~ s/TCA/TCA/g);
	my $ACount2 = ($repeats =~ s/CCA/CCA/g); # new
	$ACount += $ACount2;

	my $GCount = ($repeats =~ s/TCG/TCG/g);
	my $GCount2 = ($repeats =~ s/CCG/CCG/g); # new
	$GCount += $GCount2;

	my $TCount = ($repeats =~ s/TCT/TCT/g);
	my $TCount2 = ($repeats =~ s/CCT/CCT/g); # new
	$TCount += $TCount2;

	my $results = "";
	if ($ACount != "") { $results .= "A$ACount"; }
	if ($GCount != "") { $results .= "G$GCount"; }
	if ($TCount != "") { $results .= "T$TCount"; }

	return ($results, $motif, $repeats);
}

#Find the subtype specific R repeats
sub FindRRepeats {

#IIa = ACATCA = TS
#Ia = AAA/G ACG GTG GTA AGG
#Ia = AAA/G ACG GTG (GTA)/A AGG (from sample)
#If orignal =   C/AAG AA/GG GCA
#If dawn = C/AAG AA/GG GCA A/GAG AAG

#IXa Lihua = ATT CTG GTA CTG AAG ATA
#	     ATT CTG GTG CTG GAG GTA
#	     GTT CTG GTA CTA GCG ATT


#IXa sample= A/GTT CTG GTA/G CTG/A A/GA/CG A/GTA/T

	@_[1] =~ s/-//g; #remove gaps from alignments
	my $rCount; # is empty if there are no matched patterns
	
	if (@_[0] eq "Ia") {

		$rCount = (@_[1] =~ s/A[AG][AG]ACGGTG(GT)?AAGG//g);

	} elsif (@_[0] eq "If") {

		$rCount = (@_[1] =~  s/[CA]AGA[AG]GGCA//g); #simplest, could false positive
		#$rCount = (@_[1] =~ s/[CA]AGA[AG]GGCAA?(GTG)?[AG]AGA?AAG//g); #more stringent/more false negatives

	} elsif (@_[0] eq "IIa") {

		$rCount = (@_[1] =~ s/ACATCA//g);

	} elsif (@_[0] eq "IXa") {

		##IXa sample=         A/GTT  CTG GTA/G  CTG/A  A/GA/CG   A/GTA/T
		$rCount = (@_[1] =~ s/[AG]TTCTGGT[AG]CT[GA][AG][AC]G[AG]T[AT]//g);
	} elsif (@_[0] eq "XIIIa") { #C. erinaci

		#ACATCA same as IIa
		$rCount = (@_[1] =~ s/ACATCA//g);
	}

	if ($rCount ne "") {

		return("R$rCount");
	}
}

