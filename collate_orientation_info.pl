#!/usr/bin/perl/
use strict; use warnings;

die "parameter needed: 1 if trinity data only, 0 if trinity and 454 data\n" if (@ARGV != 1);

my ($trin_only) = @ARGV;

my $trinity_file = 'bom_trinity_orientation.txt';
my $trin454_file = 'bom454_orientation.txt';
my $para_names = 'Bb454_paralogues_names.txt';
my $in;
my $names;
my $line;
my $nameline;
my @fields;
my @namefields;
my $contig;
my $found;
my $outfile = 'bom_orientation_oldName.txt';
my $out;
my $orient;
my $comp;
my $subc;

open ($in,"< $trinity_file") or die "Can't open $trinity_file\n";
open ($out,"> $outfile") or die "Can't create $outfile\n";

while ($line = <$in>) {
	chomp($line);
	@fields = split ('\t',$line);
	if (!($fields[2] eq 'n')) {
		if ($fields[0] =~ /c(\d+)_g(\d+)/) {    # new contig names
			$comp = $1;
			$subc = $2;
			$fields[1] =~ /i(\d+)/;
			$contig = 'comp' . $comp . '_c' . $subc . '_seq' . $1;
		}
		else {
			$contig = $fields[0] . '_seq' . $fields[1];    # old contig names
		}
		print $out "$contig\t$fields[2]\n";
	}
}

close ($in);

if ($trin_only == 0) {
	open ($in,"< $trin454_file") or die "Can't open $trin454_file\n";
	# sample lines:
	# comp1_c0        0       Contig1_rev     -	12
	# comp1_c0        1       comp35024_c0_seq1       -	12

	my $count = 0;

	while ($line = <$in>) {
		chomp($line);
		@fields = split ('\t',$line);
		if (!($fields[2] =~ /n/)) {        # seqno = 0 and orientation not equal 'na'
			if ($fields[1] =~ /rev/) {			# if sequence reverse_complemented, invert orientation
				if ($fields[2] eq '+')  { $orient = '-' }
				else 			{ $orient = '+' }		
			}
			else 				{ $orient = $fields[2] }
			$found = 0;
	
		# look for the contig in the paralogue name file and use new 'comp...' name if found
		open ($names,"< $para_names") or die "Can't open $para_names\n";
		# sample lines:
		# comp31270_c0_seq100     Contig10018
		# comp24219_c0_seq100     Contig10025
		# comp39500_c0_seq100     Contig10059

			while ($nameline = <$names>) {
				chomp($nameline);
				@namefields = split ('\t',$nameline);
				if ($namefields[1] eq $fields[0]) {
					$contig = $namefields[0];
					$found = 1;
					last;
				}
			}
			close ($names);
			unless ($found) { $contig = $fields[0] }
			print $out "$contig\t$orient\n";
			$count++;
		}
	}
}

close ($in);
close ($out);
