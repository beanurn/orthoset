#!/usr/bin/perl
# robust_good_groups.pl
# finds all groups of Bb, Bvv and Bo that are consistently identified across IVs
use strict; use warnings;

my $groups_in = "groups_1.5_GoodGroups.txt";
my $test1_in = "groups_0.5.txt";
my $test2_in = "groups_5.0.txt";
my $groups_out = "RobustGroups_4tax.txt";
my $ingroups;
my $test;
my $outgroups;

my $contig;
my @array;
my @bombina_core;
my @outputarray;

my @testgroup;
my $target;
my $i;
my $j;

open ($ingroups, "< $groups_in") or die "Can't find $groups_in.\n";
open ($outgroups, "> $groups_out") or die "Can't create $groups_out.\n";

$i = 0;

while (my $line = <$ingroups>) {
	chomp ($line);
	@array = split(" ",$line);
	while (@array) {
		$contig = shift (@array);
		# the list of regexes in the following line must match the analogous list further down in test_if_robust() 
		push (@bombina_core, $contig) unless ($contig =~ m/Xtr/ || $contig =~ m/G_/);
		push (@outputarray, $contig) unless ($contig =~ m/G_/);
	}
	@bombina_core = sort {$a cmp $b} @bombina_core;
	if (@bombina_core == 0) {
		print "for group $i the core array is empty\n";
		die;
	}
	my $drop = 0;

	$drop = test_if_robust ($test1_in);
	$drop = test_if_robust ($test2_in);

	if (!($drop)) {
		foreach $contig (@outputarray) {
			print $outgroups "$contig ";
		}
		print $outgroups "\n";
	}
	$i++;
	if ($i % 10 == 0) {
		print "$i\t";
	}
	@bombina_core = ();
	@outputarray = ();
}


close ($ingroups);
close ($outgroups);


sub test_if_robust {

	my @newarray = ();
	my @testarray;
	my $nomatch = 0;

	my $infilename = $_[0];

	
	open ($test, "< $infilename") or die "Can't open $infilename.\n";
	while (my $testline = <$test>) {
		chomp ($testline);
		if ($testline =~ m/\Q$bombina_core[0]\E/) {
			@newarray = split (" ",$testline);
			last;
		}
	}
	close ($test);	
	
	if (@newarray == 0) { $nomatch = 1 }
	else {
		while (@newarray) {
			my $ct = shift (@newarray);
			push (@testarray, $ct) unless ($ct =~ m/Xtr/ || $ct =~ m/G_/);
		}
		if (@bombina_core != @testarray) { $nomatch = 1 }
		else {
			@testarray = sort {$a cmp $b} @testarray;
			for (my $j=0;$j < @bombina_core; $j++) {
				if (!($bombina_core[$j] eq $testarray[$j])) {
					$nomatch = 1;
					last;
				}
			} 
		}
	}

	@testarray = ();
	return $nomatch;
}


