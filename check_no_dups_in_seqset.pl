#!/usr/bin/perl/
use strict; use warnings;
# the infile name must contain the pattern 'set' so that the outfile name is correctly named

my %seqcount;
my ($infilename) = @ARGV;
my $outfilename = $infilename;
my $dup_count = 0;
my $name;
my @more_than_one;
my %found;
my $flag;
my $dup;
my $removed = 0;

open (my $in,"< $infilename") or die "Can't open $infilename\n";

while (my $line = <$in>) {
	chomp($line);
	if ($line =~ />/) {
		$line =~ s/>//;
		$seqcount{$line}++;
	}
}

foreach $name (keys %seqcount) {
	if ($seqcount{$name} > 1) {
		push @more_than_one, $name;
		$dup_count++;
	}
}

print "$dup_count sequences have more than one copy\n";
foreach $name (@more_than_one) {
	print "$name\t$seqcount{$name}\n";
}

close ($in);

if ($dup_count > 0 ) {
open ($in,"< $infilename") or die "Can't open $infilename\n";
$outfilename =~ s/set/set_no_dups/;
open (my $out, "> $outfilename") or die "Can't create $outfilename\n";


%seqcount = ();

while (my $line = <$in>) {
	chomp($line);
	if ($line =~ />/) {
		$flag = 0;
		foreach $dup (@more_than_one) {
			if ($line =~ /$dup$/) {      # $line includes >, so this is part of the hash key, but that's ok.
				$seqcount{$line}++;
				if ($seqcount{$line} > 1) {
					$flag = 1;
					$removed++;
				}
				last;
			}	
		}
	}
	unless ($flag) { print $out "$line\n" }
}

print "$removed duplicates were removed\n";

close ($in);
close ($out);
}
