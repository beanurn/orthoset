#!/usr/bin/perl/
use strict; use warnings;


die "Two arguments needed: pep filename and orthoMCL taxon prefix\n" unless (@ARGV == 2);
my ($infile,$prefix) = @ARGV;

my %keep;
my @startbuffer;
my @endbuffer;
my @orfbuffer;
my @orientbuffer;
my $seqbuffer;
my @linebuffer;
my $in;
my $line;
my $out;
my %known_orientation;
my @fields;
my $sequence;
my $i;
my $orf;

my $orientfile = 'orient_trinity_orientation_ed.txt';

open ($in,"< $orientfile") or die "Can't open $orientfile\n";

while ($line = <$in>) {
	chomp($line);
	@fields = split ('\t',$line);
	$known_orientation{$fields[0]} = $fields[1];
}
close ($in);

my $counter = 0;
open ($in,"< $infile") or die "Can't open $infile\n";
my $first = 1;
$seqbuffer = '';

while ($line = <$in>) {
	if ($line =~ />/) {
		chomp($line);
		if ($line =~ /(Contig\d+_cap3):/) { $sequence = $1 }
		else {
			$line =~ /(\w+\d+):/;
			$sequence = $1;
		}
		if (!($sequence eq $seqbuffer)) {
			if ($first) { $first = 0 }
			else {
				if (@startbuffer > 1) {
					select_orfs ();
				}
				else {
					$keep{$orfbuffer[0]} = 1;
				}	
				@startbuffer = ();
				@endbuffer = ();
				@orientbuffer = ();
				@orfbuffer = ();
				$seqbuffer = $sequence;
			}
		}
		$line =~ /(\d+)-(\d+)/;
		push @startbuffer, $1;
		push @endbuffer, $2;
		if ($line =~ /\+/) 	{ push @orientbuffer, '+' }
		else 			{ push @orientbuffer, '-' }
		$line =~ /\) /;
		push @orfbuffer, $';
		$counter++;
#		die if ($counter == 10); 
	}
}


if (@startbuffer > 1) {
	select_orfs ();
}
else {
	$keep{$orfbuffer[0]} = 1;
}	


close ($in);

open ($in,"< $infile") or die "Can't open $infile\n";
my $outfile = $infile . '_selected';
open ($out,"> $outfile") or die "Can't create $outfile\n";

my $flag = 0;

while ($line = <$in>) {

	chomp($line);
	if ($line =~ />/) {
		if ($flag) {
			for ($i=0;$i<@linebuffer;$i++) {
				print $out "$linebuffer[$i]\n";
			}
			@linebuffer = ();			
		}
		$line =~ /\) /;
		$orf = $';
		if ($keep{$orf} == 1)   { $flag = 1 }
		else 			{ $flag = 0 }
		if ($flag) {
			$orf = '>' . $prefix . '|' . $orf; 
			push @linebuffer, $orf;
		}
	}
	elsif ($flag) {
		$line =~ s/\*//g; 
		push @linebuffer, $line;
	 }
}

if ($flag) {
	for ($i=0;$i<@linebuffer;$i++) {
		print $out "$linebuffer[$i]\n";
	}
}


close ($in);
close ($out);

sub select_orfs {

my $i;
my $pluscount = 0;
my $minuscount = 0;
my $true_orientation;

	if (exists $known_orientation{$seqbuffer}) {
		for ($i=0;$i<@startbuffer;$i++) {
			if ($orientbuffer[$i] eq $known_orientation{$seqbuffer}) {
				$keep{$orfbuffer[$i]} = 1;
			}
			else { $keep{$orfbuffer[$i]} = 0 }
		}
	}
	else {
		for ($i=0;$i<@startbuffer;$i++) {
			if ($orientbuffer[$i] eq '+') {
				$pluscount += $endbuffer[$i] - $startbuffer[$i];
			}
			else {	$minuscount += $endbuffer[$i] - $startbuffer[$i] }
		}
		if ($pluscount >= $minuscount) 	{ $true_orientation = '+' }
		else				{ $true_orientation = '-' }

		for ($i=0;$i<@startbuffer;$i++) {
			if ($orientbuffer[$i] eq $true_orientation) {
				$keep{$orfbuffer[$i]} = 1;
			}
			else { $keep{$orfbuffer[$i]} = 0 }
		}
	}
}

	
	
	
