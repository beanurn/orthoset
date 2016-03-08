#!/usr/bin/perl
use strict; use warnings;

my @inputdata;
my @fields,
my @data_sorted;
my $seqname;


my ($infilename) = @ARGV;
open (my $in, "< $infilename") or die "Can't open $infilename\n";
my $outfilename = $infilename . ".sorted";
$outfilename =~ s/results/fpkm/;
open (my $out, "> $outfilename") or die "Can't create $outfilename\n";

my $headerline = <$in>;

while (my $line = <$in>) {
	chomp($line);	
	@fields = split ('\t',$line);
	$seqname = $fields[0];
	if ($fields[0] =~ /_g/) {
		$seqname =~ s/^c/comp/;
		$seqname =~ s/_g/_c/;
		$seqname =~ s/_i/_seq/;
	}
	else { $seqname = $fields[0] }
	$seqname =~ /comp(\d+)_c(\d+)/;
	push @inputdata, { comp_number => $1, subcomp => $2, name => $seqname, fpkm => $fields[6] };
}

close ($in);

@data_sorted = sort { $a->{comp_number} <=> $b->{comp_number} || $a->{subcomp} <=> $b->{subcomp} } @inputdata;

for (my $i=0; $i<@data_sorted; $i++ ) {
	print $out "$data_sorted[$i]->{name}\t$data_sorted[$i]->{fpkm}\n";
}

close ($out);
	
	
