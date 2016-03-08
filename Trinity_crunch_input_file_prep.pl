#!/usr/bin/perl
use strict; use warnings;
# this takes a list of Trinity contig names in the new contig name format as input, converts these to old style names  and outputs them sorted in numeric order;
my @inputdata;
my @data_sorted;
my $contigname;
my @fields;
my $refseq;

my ($infilename) = @ARGV;
open (my $in, "< $infilename") or die "Can't open $infilename\n";
my $outfilename = $infilename . ".sorted";
open (my $out, "> $outfilename") or die "Can't create $outfilename\n";


while (my $line = <$in>) {
	chomp($line);	
	@fields = split ('\t',$line);
	if ($fields[1] =~ /\d+/) { $refseq = 1 }
	else 			 { $refseq = 0 }
	$fields[0] =~ /c(\d+)_g(\d+)_i(\d+)/;
	push @inputdata, { comp_number => $1, subcomp => $2, seq => $3, ref => $refseq };
}

close ($in);

@data_sorted = sort { $a->{comp_number} <=> $b->{comp_number} || $a->{subcomp} <=> $b->{subcomp} || $a->{seq} <=> $b->{seq} } @inputdata;

for (my $i=0; $i<@data_sorted; $i++ ) {
	$contigname = 'comp' . $data_sorted[$i]->{comp_number} . '_c' . $data_sorted[$i]->{subcomp} . '_seq' . $data_sorted[$i]->{seq};
	print $out "$contigname\t$data_sorted[$i]->{ref}\n";
}

close ($out);
	
	
