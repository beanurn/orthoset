#!/usr/bin/perl
# - Trinity_condense_alignments.pl - takes a MUSCLE fasta alignment as input and removes small +-terminal gaps
use strict; use warnings;

my ($infilename) = @ARGV;

open (my $in,"< $infilename") or die "Can't open $infilename\n";

my @headerline;
my @sequence;
my $i = -1;
my $alignment_length;
my @current_gap;
my $seq_count;
my @blocks;		# holds array references ($hashrefs0, $hashrefs1) to arrays of hashes
my @hashrefs0;  	# two empty array references to store the hash references of blocks in
my @hashrefs1;
my $blockcount;
my $a;
my $b;
my $n;
my $index;
my $j;
my $k;

@blocks = (		# blocks now holds the two array references which will be populated with hash references
	\@hashrefs0,
	\@hashrefs1,
);

while (my $line = <$in>) {
	chomp ($line);
	if ($line =~ '>') {
		$i++;
		$headerline[$i] = $line;
	}
	else {
		$sequence[$i] .= $line;
	}
}

close ($in);


$seq_count = $i + 1;
$alignment_length = length($sequence[0]);


for ($i=0;$i<$seq_count;$i++) {
	$current_gap[$i] = 0;
	if (substr ($sequence[$i],0,1) =~ /[ACTG]/) {

		# The use of 'my' on the next line is critical, because new memory is being allocated each time.
		# Without 'my' the hash reference would point to the same memory address each time

		my %hash = (start => 0, length => 1, aligned => 0, id => 0, gap_left => 0, gap_right => 0);
		push @{$blocks[$i]}, \%hash; 
	}
	else	{ $current_gap[$i] = 1 }
}


for ($n=1;$n<$alignment_length;$n++) {
	for ($i=0;$i<$seq_count;$i++) {
		$index = scalar @{$blocks[$i]} - 1;		  # last index in array
		if (substr ($sequence[$i],$n,1) =~ /[ACTG]/) {    # now in block
			if ($current_gap[$i] > 0) {		  # new block
				if ($index >= 0) {
					$blocks[$i][$index]->{gap_right} = $current_gap[$i];
				}
				$index++;
				my %newhash = (start => $n, length => 1, aligned => 0, id => 0, gap_left => $current_gap[$i], gap_right => 0);
				push @{$blocks[$i]}, \%newhash;
				$current_gap[$i] = 0;
			}	
			else {					# existing block
				$blocks[$i][$index]->{length} += 1;
			}
		}
		else {						# now in gap
			if ($current_gap[$i] > 0) {		# existing gap
				$current_gap[$i]++;
			}
			else {					# new gap
				if ($blocks[$i][$index]->{aligned} > 0) {
					$blocks[$i][$index]->{id} = $blocks[$i][$index]->{id} / $blocks[$i][$index]->{aligned};
				}
				$current_gap[$i] = 1;
			}
		}
	}

	if (substr($sequence[0],$n,1) =~ /[ACTG]/ and substr($sequence[1],$n,1) =~ /[ACTG]/) {
		$a = scalar @{$blocks[0]} - 1;
		$blocks[0][$a]->{aligned} += 1;
		$b = scalar @{$blocks[1]} - 1;
		$blocks[1][$b]->{aligned} += 1;
		if (substr($sequence[0],$n,1) eq substr($sequence[1],$n,1)) {
			$blocks[0][$a]->{id} += 1;			
			$blocks[1][$b]->{id} += 1;	
		}
	}
}


for ($i=0;$i<$seq_count;$i++)  {
	$a = scalar @{$blocks[$i]} - 1;
	if ($a >= 0) {
		if ($current_gap[$i] > 0) {       	# aligned sequence ends in a gap
			$blocks[$i][$a]->{gap_right} = $current_gap[$i];
		}	
		else {					# aligned sequence ends in nucleotide			
			if ($blocks[$i][$a]->{aligned} > 0) {	
				$blocks[$i][$a]->{id} = $blocks[$i][$a]->{id} / $blocks[$i][$a]->{aligned};
			}
		}
	}
}


#************************ determine the first and last blocks in the alignment that serves as anchors **********
#************************ because of their length and high sequence identity************************************
my $left_anchor;
my $right_anchor;
my $left_done;

for ($i=0;$i<$seq_count;$i++)  {
	$left_done = 0;
	$blockcount = scalar @{$blocks[$i]};
	for ($j=0;$j<$blockcount;$j++) {
		if 	(($blocks[$i][$j]->{length} > 30 and $blocks[$i][$j]->{id} > 0.95) or 
			 ($blocks[$i][$j]->{length} > 50 and $blocks[$i][$j]->{id} > 0.85)) {
			if (!($left_done)) {	
				$left_anchor = $j;
				$left_done = 1;
			}
			$right_anchor = $j;
		}
	}
	if ($left_done) {
#		print "seq $i, left anchor: $left_anchor, right anchor: $right_anchor, blks: $blockcount\n";
		if ($left_anchor > 0) { defragment_left_side ($i, $left_anchor) }
		if ($right_anchor < $blockcount - 1) { defragment_right_side ($i, $right_anchor) }
	}
}


open (my $out,"> $infilename") or die "Can't open $infilename for writing\n";

for ($i=0;$i<$seq_count;$i++) {
	print $out "$headerline[$i]\n$sequence[$i]\n";
}

close ($out);

#***************************************************************************************************************

sub defragment_left_side {

my $j;
my $condensed = '';
my $totalgaps = 0;
my $newseq;

	my ($a,$anchor) = @_;
	
	$totalgaps = $blocks[$a][0]->{gap_left};
	for ($j=0;$j<$anchor;$j++) {
		$condensed .= substr($sequence[$a],$blocks[$a][$j]->{start},$blocks[$a][$j]->{length});
		$totalgaps += $blocks[$a][$j]->{gap_right};
	}
	$newseq = '-' x $totalgaps;
	$newseq .= $condensed;
	$newseq .= substr ($sequence[$a],$blocks[$a][$anchor]->{start});
	$sequence[$a] = $newseq;
}					 

sub defragment_right_side {

my $j;
my $condensed = '';
my $totalgaps = 0;
my $newseq;
my $anchorend;
my $blks;


	my ($a,$anchor) = @_;
	
	$totalgaps = $blocks[$a][$anchor]->{gap_right};
	$blks = scalar @{$blocks[$a]};
	$anchorend = $blocks[$a][$anchor]->{start} + $blocks[$a][$anchor]->{length};
	for ($j=$anchor+1; $j<$blks; $j++) {
		$condensed .= substr($sequence[$a],$blocks[$a][$j]->{start},$blocks[$a][$j]->{length});
		$totalgaps += $blocks[$a][$j]->{gap_right};
	}
	$newseq = substr ($sequence[$a],0,$anchorend);
	$newseq .= $condensed;
	$newseq .= '-' x $totalgaps;
	$sequence[$a] = $newseq;
}					 


	
