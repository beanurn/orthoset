#!/usr/perl/bin/
use strict; use warnings;

my $line;
my @data;
my @data_sorted;
my @fields;
my $i;
my $comp;
my $new_seq;
my $trinity_seq;
my $buffer_454 = 'bla';

my $infile1 = ' Bvv454_blastn_trinity.out';   # file needs to have an .out extension!
# my $infile2 = 'Bb_trinity_single_blastn.out';
my $besthitfile = $infile1;
$besthitfile =~ s/.out/_besthit.out/;   


open (my $in1, "< $infile1") or die "Can't open $infile1\n";
open (my $best,"> $besthitfile") or die "Can't create $besthitfile\n";

while ($line = <$in1>) {
	chomp ($line);
	@fields = split ('\t',$line);
	if (!($fields[0] eq $buffer_454)) {			 # record the best blasthit for any one 454 contig in a separate file
		$comp = $fields[1];
		$comp =~ s/_seq\d+//;
		print $best "$fields[0]\t$comp\n";
	}
	$buffer_454 = $fields[0];
	$fields[1] =~ /comp(\d+)_c(\d+)_seq(\d+)/;
#	push @data, {seq454 => $fields[0], comp => $1, subc => $2, seq => $3, type => 'multiseq'};
	push @data, {seq454 => $fields[0], comp => $1, subc => $2, seq => $3, qstart => $fields[6], qend => $fields[7], s_start => $fields[8], s_end => $fields[9]};
}

close ($in1);
close ($best);
$buffer_454 = 'bla';
	
#open (my $in2, "< $infile2") or die "Can't open $infile2\n";
#while ($line = <$in2>) {
#	chomp ($line);
#	@fields = split ('\t',$line);
#	$fields[1] =~ /comp(\d+)_c(\d+)_seq(\d+)/;
#	push @data, {seq454 => $fields[0], comp => $1, subc => $2, seq => $3, type => 'singleseq'};
#}
	
# close ($in2);

@data_sorted = sort { $a->{seq454} cmp $b->{seq454} || $a->{comp} <=> $b->{comp} || $a->{subc} <=> $b->{subc} || $a->{seq} <=> $b->{seq} } @data;

my $outfile = 'Bvv_new_components.txt';
open (my $out, "> $outfile") or die "Can't create $outfile\n"; 

#my $record;
#foreach $record (@data_sorted) {
#	print $test "$record->{seq454}\t$record->{comp}\t$record->{subc}\t$record->{seq}\n";
#}
#close ($test);

my $compct = 0;
my $seqct = 0;
my $buffer_trin = 'bla';
#my $s_count = 0;
#my $m_count = 0;
my $two_source = 0;
my $reverse;
my $q_orient;
my $s_orient;

for ($i=0;$i<@data_sorted;$i++) {
	if (!($buffer_454 eq $data_sorted[$i]->{seq454})) {
#		if ($s_count + $m_count == 2) { $two_source++ }
		$buffer_454 = $data_sorted[$i]->{seq454};
		$compct++;
		$seqct = 1;
		$new_seq = "comp" . $compct . "_c0_seq" . $seqct;
		print $out "$new_seq\t$data_sorted[$i]->{seq454}\n";
#		$s_count = 0;
#		$m_count = 0;
		$buffer_trin = 'bla';
	}
	$trinity_seq = "comp" . $data_sorted[$i]->{comp} . "_c" . $data_sorted[$i]->{subc} . "_seq" . $data_sorted[$i]->{seq};
	if (!($buffer_trin eq $trinity_seq)) {
		$seqct++;
		$new_seq = "comp" . $compct . "_c0_seq" . $seqct;
		$q_orient = $data_sorted[$i]->{qend} - $data_sorted[$i]->{qstart};
		$s_orient = $data_sorted[$i]->{s_end} - $data_sorted[$i]->{s_start};
		if (($q_orient > 0 and $s_orient < 0) || ($q_orient < 0 and $s_orient > 0)) { $reverse = 1 }
		else 									    { $reverse = 0 }
		print $out "$new_seq\t$trinity_seq\t$reverse\n";
		$buffer_trin = $trinity_seq;
#		if ($data_sorted[$i]->{type} eq 'singleseq') 	{ $s_count = 1}
#		else						{ $m_count = 1}
	}
}

# if ($s_count + $m_count == 2) { $two_source++ }

close ($out);
print "number of new components: $compct\n";
