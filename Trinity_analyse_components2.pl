#!/usr/bin/perl
# - Trinity_analyse_components.pl - takes an msf alignment of sequences of a component and identifies isoforms, diverged UTRs, introns
# and the like in order to select suitable sequences for further analysis.
use strict; use warnings;


#************************************************** NUMERIC CONSTANTS *************************************************************************

my $CDS_threshold = 24; # min. number of utr bases for which a search for obscured cds is initiated and min. # bases 
			# with near perfect match to a cds block in another sequence before a block of utr is reclassified
			# as an obscured cds
my $TOTAL_FRAGMENTS = 35.0; # number of Bb Illumina fragments after QC (millions), data for the other taxa: 
				# Bb: 35.0
				# Bvv: 41.0
				# Bo: 34.4
				# Bvs: 25.2
my $MEAN_READLENGTH = 90;  # rough estimate of mean read length after QC  
my $I_PAD = 10;            # padding at either end of the search string around a suspected intron to make sure any 
			   				# existing intron/exon boundary isn't missed
my $MIN_LENGTH_3P = 18;    # the minimum length of the 3' intron signal starting at A and ending on AG (inclusive)
my $MIN_INTERVAL = 48;       # minimum length of an intron sequence between 5' and 3' signals
my $MIN_EXON_LENGTH = 15;
my $AL_WINDOW = 10;			# size of alignment window over which seq id is computed
my $ID_LOW = 0.50;		# upper threshold for low ID alignment windows, which will be counted



#************************************************ GLOBAL VARIABLES **************************************************************************
my $first = 1;
my $raw_seq;

my $component_to_follow = 'bla';

my %coverage;
my $rsem;
my $rsem_line = '';
my $rsemfile = '/users/bnurnberger/Trinity_components/bom/RSEM.isoforms.fpkm.sorted';   # will be opened at the end to add fpkm data to output file
my $rsem_component = '';                      # stores the line where the last match to a given component name occurred
						# upon next entry to sub write_output_to_file, the rsemfile is fast-forwarded to that line
						# and the search for the next matching component continues from there
my $in;
my $component;
my @fields;

my @forward_stops = ('TAA','TGA','TAG');
my @reverse_stops = ('TTA','TCA','CTA');

my $intron_info_file = 'intron_info.txt';
my $iinfo;
my $intron_start;
my $intron_end;

my $fshiftfile = 'frameshift_overlaps.txt';
my $fshift;

my $count = 0;
my $i;
my $j;
my $k;

my $position;

# raw data being read in from file
my @sequences;
my @seq_names;
my @seq_index;

my $seq_count;		# number of contigs within the component (not orfs or xt matches)
my $name;

# sequences sorted in numerial order of their Trinity index, the index itself is stored in @seq_index 
my @seqnames_sorted;
my @sequences_sorted;

my $href; 		# a reference to a hash in the @annotation array; 
my @annotation = ();

my $alignment_length;	# length of the muscle alignment string
my $matches;		# counter of matched sections in a given pairwise sequence comparison
my $in_alignment;	# flag to indicate that 
my $m;


# arrays that summarise the information about each sequence in blocks along the alignment
my %patterns;		# 0 = no sequence at start/end of alignment or gap, 
			# 100 = sequence without xt or orf support, 110 = orf support but not xt, 101 = xt support but not orf, 111 = orf+xt support
			# in case of multiple orfs per contig the second digit would reflect their count, analogous for multiple xt matches (3. digit)
my %changepos;		# position in alignment where next block starts, current block ends at the position before this one
my %status;		# cds, intron, indel, utr, noseq 
my %frames;		# for cds only: reading frame (-3, -2, -1, 1, 2, 3)

my @minusbases;
my @plusbases;

# arrays that hold the information on sequence ID

my @alignment;
my %matchstart;		# these are hashes of pairwise sequence comparisons, each pointing to an array of aligned blocks
my %matchend;
my %matched_bases;

my $sections;
my $endpoint;
my $start;
my $stop;
my $pattern;
my $current_pattern;

# summary info about each sequence
my @seq_stats;

my $seq_name;
my $suffix = 'pattern';

my $comp_count = -1;
my $print_flag = 0;
# ************************************* MAIN PROGRAM BLOCK ********************************************************


die "Usage: Trinity_analyse_components.pl - input file with names of aligned components needed\n" if (@ARGV != 1);
my ($listfilename) = @ARGV;
open (my $listfile, "< $listfilename") or die "Can't open $listfilename.\n";

open ($iinfo, ">> $intron_info_file") or die "Can't write to $intron_info_file.\n";

open ($rsem, "< $rsemfile") or die "Can't open $rsemfile\n"; 

open ($fshift, ">> $fshiftfile") or die "Can't create $fshiftfile\n";

while ($component = <$listfile>) {
	if ($comp_count % 25 == 0) { print "analysing comp $component\n" }

	chomp ($component);           	# This is the component name with an extension for the alignment (e.g. '_35' for seqs 3 and 5)
					# The extension is only needed in order to open the correct alignment file. Thereafter the extenion is  stripped away.

	my $infilename = $component . '.afa2';
	if (!(open ($in, "< $infilename"))) {
		print "$component: Can't open $infilename. $!\n";
		next;
	}

	$component =~ /_\d+/;
	$component = $`;		# Now this is a proper component name (e.g. 'comp34456_c0')

	if ($component =~ /$component_to_follow/) {
		print "opened $infilename\n";
	}	

	if (!($component eq $rsem_component)) {
		collect_coverage_info ();
	}

	read_sequence_input ();
	close ($in);
	analyse_component ();
	system ("rm $infilename");
	$comp_count++;
}

close ($rsem);
close ($fshift);
close ($listfile);
close ($iinfo);

# ********************** Read infile and write data to two arrays: @seq_names and @sequences *******************************

sub read_sequence_input {

	$first = 1;
	@sequences = ();
	@seq_names = ();
	@annotation = ();
	$raw_seq = '';
	while (my $line = <$in>) {
		chomp($line);
		if ($line =~ />/) { 
			$line =~ s/>//;
			if ($line =~ 'xt' or $line =~ 'orf') {
				@fields = split ('\t',$line);
				push @annotation, {match => $fields[0], seqindex => $fields[1], start => $fields[2], end => $fields[3], frame => $fields[4]};
			}
			else {
				push (@seq_names,$line);
				if ($first) {$first = 0 }
				else {
					push (@sequences, $raw_seq);
					$raw_seq = '';
				}
			}
		}
		else {
			$raw_seq = $raw_seq . $line;
		}
	}
	push (@sequences,$raw_seq);
}

sub analyse_component {

# ******************** Count the number of sequences and the number of orfs and xt blast matches per sequence ***********************

@seq_index = ();
@plusbases = ();
@minusbases = ();

foreach $name (@seq_names) {
	$name =~ /seq(\d+)/;
	push @seq_index, $1;
	push @plusbases, 0;
	push @minusbases, 0;
}

$seq_count = @seq_index;
@seq_index = sort {$a <=> $b} @seq_index;

@seqnames_sorted = ();
@sequences_sorted = ();


for ($i=0;$i<$seq_count;$i++) {
	$name = $component . '_seq' . $seq_index[$i];
	$position = get_raw_array_position ($name);
	push @seqnames_sorted, $seq_names[$position];
	push @sequences_sorted, $sequences[$position];
}

$alignment_length = length($sequences_sorted[0]);
initialise_seq_stats ();

for $href ( @annotation ) {
	for ($i=0;$i<$seq_count;$i++) { 
		if ($seq_index[$i] == $href->{seqindex}) {
			if ($href->{match} eq 'orf') 	{ 
				$seq_stats[$i]->{orfcount} += 1;
				$seq_stats[$i]->{orfbases} += $href->{end} - $href->{start} + 1;
			}
			else {
				$seq_stats[$i]->{xtcount} += 1;
				$seq_stats[$i]->{xtbases} += $href->{end} - $href->{start} + 1;
			}
			# note that a decision on the orientation of a component is now made in Trinity_align_components4.pl and all annotations passed to 
			# this programme are in one orientation only
			if ($href->{frame} > 0) { $seq_stats[$i]->{orientation} = '+' }
			else 			{ $seq_stats[$i]->{orientation} = '-' }  
			last;
		}
	}
}

# ****************** write the data in the following sort order: seq1, orf1(seq1), orf2(seq1),..., xt1(seq1), xt2(seq1),..., seq2,... etc *************


if (!(scalar(@sequences_sorted) == scalar(@sequences))) {
	die "raw and sorted sequence array differ in length: sequences ", scalar(@sequences), " sequences_sorted ", scalar(@sequences_sorted), "\n";
}


if ($seqnames_sorted[0] =~ /$component_to_follow/) {
	print "initialised stats for $seqnames_sorted[0]\n";
}	


if (@annotation == 0) { 

	write_output_to_file ();
}
else {

# ******************  do annotation, tidy up alignment and do preliminary interpretation ********************************************************  

my $seq_start;
my $shift_flag;
							
for ($i = 0; $i < $seq_count; $i++) {	

	$seq_start = record_presence_absence ($i);
#	write_patterns_to_screen ($i,0);    # parameter is last seq in array to be written to screen 
	add_orf_and_xt ($i, $seq_start);
	$shift_flag = check_gap_ends ($i);
	if ($shift_flag) {
		$seq_start = record_presence_absence ($i);
		add_orf_and_xt ($i,$seq_start);
	}
	# populate_status_array also records frameshifts
	populate_status_array ($i);

#	if ($minusbases[$i] < 0.1 * $plusbases[$i]) 	{ $seq_stats[$i]->{orientation} = '+' }
#	elsif ($plusbases[$i] < 0.1 * $minusbases[$i]) 	{ $seq_stats[$i]->{orientation} = '-' }
#	else						{ print "mismatch: $component - minus: $minusbases[$i], plus: $plusbases[$i]\n" }
#	print "$component, $i: minus: $minusbases[$i], plus: $plusbases[$i], orientation: $seq_stats[$i]->{orientation}\n";
		
}

# write_patterns_to_screen ($seq_count-1,1);
# output_annotated_alignment ($component);



#  ********************** refine interpretation **********************************************************************************


# then loop over sequences again and search for obscured cds...

for ($i = 0; $i < $seq_count; $i++) {	
	check_for_obscured_orf ($i);
	define_orf_range ($i);
	check_for_internal_stops ($i);
}



# after all obscured orf blocks are identified, a search for introns is triggered by gaps
# each identified intron is marked as such in all sequences having the intron and all associated gaps in the remaining sequences are marked as 'igaps'...

for ($i = 0; $i < $seq_count; $i++) {
#	for ($j=1;$j<@{$patterns{$i}}-1;$j++) {
#		if ($patterns{$i}[$j] == 0) {
#			check_for_intron ($i,$j);
#		}
#	}	
	check_for_intron ($i);
}

# all remaining gaps are then searched for likely splice sites.... - gaps that don't qualify as splice sites are labelled 'indel'
# NB: this classification is done irrespective of 3n length of gap

for ($i = 0; $i < $seq_count; $i++) {
	for ($j=1;$j<@{$patterns{$i}}-1;$j++) {
		if ($status{$i}[$j] eq 'gap') {
			check_for_splice ($i,$j);
		}
	}	
}

# write_patterns_to_screen ($seq_count-1,1);

#  ********************************  compute base identity per category in windows of $AL_WINDOW length for all pairs of sequence ********************


my $block_status;
my %min = (utr5p => 1, mm5p => 1, cds => 1, mm3p => 1, utr3p => 1);
my %max = (utr5p => 0, mm5p => 0, cds => 0, mm3p => 0, utr3p => 0);
my %mean = (utr5p => 0, mm5p => 0, cds => 0, mm3p => 0, utr3p => 0);
my %count = (utr5p => 0, mm5p => 0, cds => 0, mm3p => 0, utr3p => 0);
my %blocks = (utr5p => 0, mm5p => 0, cds => 0, mm3p => 0, utr3p => 0);
my $current = 'null';   # current block status
my $blocktype;
my $id;
my $id_count;
my $id_sum;
my $mean_id;
my @worst_id;
my $cds_one;

for ($i=0;$i<$seq_count;$i++) {
	push @worst_id, 10;
}
 
for ($i=0;$i<$seq_count-1;$i++) {
	for ($j=$i+1;$j<$seq_count;$j++) {
		foreach $blocktype (keys %min) {
			$min{$blocktype} = 1;
			$max{$blocktype} = 0;
			$count{$blocktype} = 0;
		}
		$id_sum = 0;
		$id_count = 0;	
		$cds_one = 0;
		
		my $aseq = $sequences_sorted[$i];
		my $bseq = $sequences_sorted[$j];
		for ($k=0; $k<$alignment_length-$AL_WINDOW; $k++) {
			if (substr($aseq,$k,$AL_WINDOW) =~ /[ACTG]{$AL_WINDOW}/ and substr($bseq,$k,$AL_WINDOW) =~ /[ACTG]{$AL_WINDOW}/) {
				$block_status = get_per_block_status ($i,$j,$k);
				if (!($block_status eq 'null')) {
					if (!($block_status eq $current)) { $blocks{$block_status}++ }
					$id = compute_id ($i,$j,$k);
					$mean{$block_status} += $id;
					if ($id > $max{$block_status}) { $max{$block_status} = $id }
					if ($id < $min{$block_status}) { $min{$block_status} = $id }
					$count{$block_status}++;   	# each window provides one estimate of id and one new base
					$id_count++;
					$id_sum += $id;
					if ($block_status eq 'cds' and abs($id - 1.0) < 0.01) { $cds_one++ }
				}
				$current = $block_status;
			}
			else {  $current = 'null' }
		}
				
		foreach $blocktype (keys %count) {
			# each new block of a given type adds $Al_WINDOW-1 extra bases at the start, which need to be added here
			# to give the total count of bases per type. This is done after the mean has been computed
			if ($count{$blocktype} > 0) { 
				$mean{$blocktype} = $mean{$blocktype} / $count{$blocktype};
				if ($blocktype eq 'cds') { $cds_one = $cds_one / $count{$blocktype} }
				$count{$blocktype} += $AL_WINDOW-1 * $blocks{$blocktype}
			}
		}
		if ($id_count > 0) {
			$mean_id = $id_sum / $id_count;	
			if ($mean_id < $worst_id[$i]) { 
				record_alignment_data ($i, $j, \%min, \%max, \%count, \%mean, $cds_one);
				$worst_id[$i] = $mean_id;
			 }  
			if ($mean_id < $worst_id[$j]) { 
				record_alignment_data ($j, $i, \%min, \%max, \%count, \%mean, $cds_one); 
				$worst_id[$j] = $mean_id;
			} 
		}
	}
}


# ****************** count bases in the different status categories, check cds/splice blocks for length = multiples of 3 *************

my $block_length;
my $var;

for ($i=0; $i < $seq_count; $i++) {
	compact_seq_annotation_arrays ($i,'status');   #  compact based on status only, patterns loses its significance after this step 
	for ($j=0;$j < @{$status{$i}}; $j++) {
		if (!($status{$i}[$j] eq 'igap')) {                   # data on igaps and plain is not part of output
			if ($j == 0) 	{ $block_length = $changepos{$i}[$j] }	
			else		{ $block_length = $changepos{$i}[$j] - $changepos{$i}[$j-1] }
			$var = $status{$i}[$j] . '_count';
			$seq_stats[$i]->{$var} += 1;
			$var = $status{$i}[$j] . '_bases';
			$seq_stats[$i]->{$var} += $block_length;
		}
	}
}

if ($component =~ /$component_to_follow/) {
	print "final sequence analysis\n";
	for ($i=0;$i<$seq_count;$i++) {
		print "seq $seq_index[$i]: ";
		for ($j=0;$j<@{$changepos{$i}};$j++) {
			if (defined($status{$i}[$j]) and defined($changepos{$i}[$j])) {
				print "$status{$i}[$j] $changepos{$i}[$j] ";
			}
		}
		print "\n\n";
	}
}

write_output_to_file ();

} # end of block dealing with components that have at least 1 orf or xt

} # end of sub process_component


# ********************************************* *******************************************+
	
sub collect_coverage_info {

my $line;
my $seqno;
my $found = 0;
my @fields;

	%coverage = ();			# empty coverage hash
	if (length($rsem_line) > 0) {   # check the contents of the last line that was read from file on the previous call to this function
		@fields = split ('\t',$rsem_line);
		if ($fields[0] =~ /$component/) {
			$fields[0] =~ /_seq/;
			$seqno = $';
			$coverage{$seqno} = $fields[1] * $TOTAL_FRAGMENTS * $MEAN_READLENGTH / 1000;
			$rsem_component = $component;
			$found = 1;
		}
	}
	while ($line = <$rsem>) {
		chomp ($line); 
		@fields = split ('\t',$line);
		if ($fields[0] =~ /$component/) {
			$fields[0] =~ /_seq/;
			$seqno = $';
			$coverage{$seqno} = $fields[1] * $TOTAL_FRAGMENTS * $MEAN_READLENGTH / 1000;
#			print "$component\t$seqno\t$coverage{$seqno}\n";
			$rsem_component = $component;
			$found = 1;
		}
		elsif ($found)  { 
			$rsem_line = $line;     # This the first line past the current component. It is stored and checked first when this function is called next.
			last;
		}
	}
}
 			

sub extract_seq_number {
	
	my ($string) = @_;
	$string =~ /_seq/;
	$string = $';
	if ($string =~ /\./) {
		$string = $`;
	}
	return $string;
}

sub get_raw_array_position {
	
	my $i;
	my $pos;

	my ($string) = @_;
	
	my $found = 0;
	for ($i = 0; $i < @seq_names; $i++) {
		if ($string eq $seq_names[$i]) {
			$pos = $i;
			$found = 1;
		}
	}
	if (!($found)) {
		die "Couldn't find sequence name $string in raw array of seq_names\n";
	}
	return $pos;
} 

sub initialise_seq_stats {

	my $i;
	my $j;
  	my $k;
	my $seqbases;

	@seq_stats = ();
	
	for ($i=0; $i < $seq_count; $i++) {
		$seqbases = 0;
		for ($j=0; $j< $alignment_length; $j++) {
			if (substr($sequences_sorted[$i],$j,1) =~ /[ACTG]/) { $seqbases++ }
		}

		push @seq_stats, {
			seqname => $seqnames_sorted[$i],
			basecount => $seqbases,
			orfcount => 0,
			orfbases => 0,
			xtcount => 0,
			xtbases => 0,
			cporf_count => 0,
			cporf_bases => 0,
			start_codon => 0,
			stop_codon => 0,
			int_stop => 0,
			utr_count => 0,
			utr_bases => 0,
			intron_count => 0,
			intron_bases => 0,
			indel_count => 0,
			indel_bases => 0,
			splice_count => 0,
			splice_bases => 0,
			cds_count => 0,
			cds_bases => 0,
			plain_count => 0,
			plain_bases => 0,
			utr5p_aligned => 0,
			utr5p_min => 1.0,
			utr5p_max => 0,
			utr5p_mean => 0,
			mm5p_aligned => 0,
 			mm5p_min => 1.0,
			mm5p_max => 0,
			mm5p_mean => 0,
			cds_aligned => 0,
			cds_min => 1.0,
			cds_max => 0,
			cds_mean => 0,
			cds_100 => 0,
			mm3p_aligned => 0,
			mm3p_min => 1.0,
			mm3p_max => 0,
			mm3p_mean => 0,
			utr3p_aligned => 0,
			utr3p_min => 1.0,
			utr3p_max => 0,
			utr3p_mean => 0,
			worstMatch => 0,
			frameshift => 0,
			orientation => 'na',
			cds_start => $alignment_length,
			cds_end => 0,
			orf_start => $alignment_length,
			orf_end => 0,
		};
	}				
}

sub record_presence_absence  {

	my $start;
	my ($seq) = @_;
	
	@{$patterns{$seq}} = ();
	@{$changepos{$seq}} = ();
	@{$status{$seq}} = ();
	@{$frames{$seq}} = ();
	$start = -1;
	
	# establish the first pattern
	if (substr($sequences_sorted[$seq],0,1) =~ /[ACTG]/) 	{ $current_pattern = 100 }
	else 						   	{ $current_pattern = 0 }

	if ($current_pattern == 100) { $start = 0 }
	# loop over all remaining positions in the alignment and record changes of the pattern
	for ($j = 1; $j < $alignment_length; $j++) {	
		if (substr($sequences_sorted[$seq],$j,1) =~ /[ACTG]/)	{ $pattern = 100 }
		else 							{ $pattern = 0 }


		# at a change in the alignment pattern, record the previous pattern in the three arrays: patterns and changepos

		if (!($pattern eq $current_pattern)) {  
			push @{$patterns{$seq}}, $current_pattern;
			push @{$changepos{$seq}}, $j;
			push @{$frames{$seq}}, 0;
			if ($pattern == 100 and $start == -1) { $start = $j }
			$current_pattern = $pattern;
		}
	}

	push @{$patterns{$seq}}, $pattern;
	push @{$changepos{$seq}}, $j;
	push @{$frames{$seq}}, 0;
	if ($start == -1 and $pattern == 100) { $start = $changepos{$seq}[$j-1] }
	
	my $blockcount = scalar @{$patterns{$seq}};
#	print "completed first scan of sequence $seq, blocks: $blockcount, end of last block: $changepos{$seq}[$blockcount-1]\n";
	if ($start > -1) { return $start }
	else		 { die "record_presence_absence(): Can't find the start of sequence $seq\n"; }

}


sub add_orf_and_xt {

my $gapbases;
my $start_block;
my $end_block;
my $code;
my $block_end;
my $newblocks;
my $endpoint;
my $feature_start;
my $feature_end;

	my ($seq, $seq_start) = @_;

	for my $href (@annotation)  {
		if ($seq_index[$seq] == $$href{seqindex}) {      		# check if annotation belongs to current sequence
			$gapbases = 0;
			$start_block = -1;
			$end_block = -1;

# The coordinates of the start and end of a feature need to be translated into the coordinates of the alignment so that the patterns and changepos arrays can
# be updated accordingly

			for ($j=0;$j<@{$patterns{$seq}};$j++) {	

# translate the position of each block end into the local sequence coordinates and look for the blocks in which
# feature start and end can be found, store these blocks and reset the feature{start} and feature{end} to the alignment coordinates. 

				if ($changepos{$seq}[$j] > $seq_start) {
					if ($patterns{$seq}[$j] == 0) {
						$gapbases += $changepos{$seq}[$j] - $changepos{$seq}[$j-1];
					}
					else {
						$block_end = $changepos{$seq}[$j] - $seq_start - $gapbases;	# sequence position in local sequence coordinates (1. pos of next block)
						if ($block_end > $$href{start} and $start_block == -1) {		
							$start_block = $j;					# identify start block for the feature
							$feature_start = $$href{start} + $seq_start + $gapbases; # replace local with alignment coordinates in feature hash
						}
						if ($block_end > $$href{end} and $end_block == -1) {		# do the same for the end block
							$end_block = $j;
							$feature_end = $$href{end} + $seq_start + $gapbases + 1; # add to get the 1. pos of the next block
														# this is in line with how $$href{end} is used subsequently
						}
					}
				}
			}
			if ($end_block ==  -1) {
				$end_block = @{$patterns{$seq}} - 1;
				$feature_end = $alignment_length;
			}

			if ($start_block == -1) {
				print "$component, feature start not found: length $seq_stats[$seq]->{basecount}, annotation $$href{start}-$$href{end}\n";
			}
			else {
				if ($feature_start < $seq_stats[$seq]->{cds_start}) { $seq_stats[$seq]->{cds_start} = $feature_start }
			}
			if ($feature_end > $seq_stats[$seq]->{cds_end}) { $seq_stats[$seq]->{cds_end} = $feature_end }
#			print "add_orf_and_xt: seq $seq, cds_start $seq_stats[$seq]->{cds_start}\n";
#			print "add_orf_and_xt: seq $seq, cds_end $seq_stats[$seq]->{cds_end}\n";

#			print "new feature data: $$href{match}, $$href{start}-$$href{end}, startblock: $start_block, endblock: $end_block, frame: $$href{frame}\n";

# 3. Once start and end blocks and start and end positions of the feature are known in alignment coordinates, loop over the blocks again
# and insert the feature into the patterns array. Note that the number of elements in the changepos and patterns arrays are likely to increase as a result.
# Whenever that happens, the loop over blocks has to be re-entered from the end of the last modified block onwards with new array length _and_ new
# index for the block in which the feature ends ($end_block)


			$endpoint = 0;
			
#			write_patterns_to_screen ($seq,0); 
		
			while ($endpoint < $alignment_length) {
				for ($j=0;$j<@{$patterns{$seq}};$j++) {
					if ($changepos{$seq}[$j] > $endpoint and $j >= $start_block and $j <= $end_block and $patterns{$seq}[$j] > 0) {
						$endpoint = $changepos{$seq}[$j];
						if ($$href{match} eq 'orf')  	{ $code = $patterns{$seq}[$j] + 10 }
						else		 		{ $code = $patterns{$seq}[$j] + 1 }
						my $block_buffer = $j;
						$newblocks = insert_feature ($code,$feature_start,$feature_end,$seq,$j,$$href{frame});
						if ($newblocks > 0) {
							$end_block += $newblocks;
							last;
						} 
					}
					if ($j == @{$patterns{$seq}}-1) { $endpoint = $alignment_length }
				}
			}
		}  # if correct sequence
		
	}  # loop over annotations
	
#  this is left out here for now, because compacting based on @patterns only risks losing frameshift information
#	compact_seq_annotation_arrays ($seq,'no_status');   #  compact only @patterns and @changepos, not @status (yet to be assigned) 

}

sub check_gap_ends {

my $k;
my $gap_fixed = 0;
my $nleft;
my $nright;
my $n_shift;
my $remainder;
my $shift_bases;

	my ($seq) = @_;

	for ($k=1;$k<@{$patterns{$seq}}-1;$k++){
		if ($patterns{$seq}[$k] == 0) {
			if ($patterns{$seq}[$k-1] > 100) {   
				if ($k-1 == 0) 	{ $nleft = $changepos{$seq}[$k-1] }
				else 		{ $nleft = $changepos{$seq}[$k-1] - $changepos{$seq}[$k-2] };
				$remainder = $nleft % 3;
				if ($remainder > 0) {
					$nright = $changepos{$seq}[$k+1] - $changepos{$seq}[$k];
					if ($nleft > $nright) {						
						$n_shift = 3 - $remainder;
						$shift_bases = substr ($sequences_sorted[$seq],$changepos{$seq}[$k],$n_shift);
						substr ($sequences_sorted[$seq],$changepos{$seq}[$k-1],$n_shift) = $shift_bases;
						substr ($sequences_sorted[$seq],$changepos{$seq}[$k],$n_shift) = "-" x $n_shift;
					}
					else {
						$n_shift = $remainder;
						$shift_bases = substr($sequences_sorted[$seq],$changepos{$seq}[$k-1]-$n_shift,$n_shift); 
						substr ($sequences_sorted[$seq],$changepos{$seq}[$k]-$n_shift,$n_shift) = $shift_bases;
						substr ($sequences_sorted[$seq],$changepos{$seq}[$k-1]-$n_shift,$n_shift) = "-"  x $n_shift;
					}
					$gap_fixed = 1;
				}
			}
		}
	}
	return $gap_fixed;
}
					
			
sub populate_status_array {

my $frame_buffer = 0;

	my ($seq) = @_;	


	#  NB: all unannotated sequences are now labelled as 'plain' here. 'utr' is assigned further down only to blocks outside the (start codon,stop codon) range					

	for ($j=0;$j < @{$patterns{$seq}}; $j++) {
		if ($patterns{$seq}[$j] == 0) {
			if ($j == 0)	{ push @{$status{$seq}}, 'noseq' }
			else		{ push @{$status{$seq}}, 'gap' }
		}
		elsif ($patterns{$seq}[$j] == 100) 	{ push @{$status{$seq}}, 'plain' }
		else 					{ push @{$status{$seq}}, 'cds' }

		if ($patterns{$seq}[$j] > 100) {
			if ($frame_buffer == 0) {  $frame_buffer = $frames{$seq}[$j] }
			else {
				if ($frame_buffer != $frames{$seq}[$j]) { 
					$seq_stats[$seq]->{frameshift} += 1;
					$frame_buffer = $frames{$seq}[$j];
			 	}
			}
		}	 
	}
	if ($status{$seq}[$j-1] eq 'gap') { $status{$seq}[$j-1] = 'noseq' }
}


sub check_for_obscured_orf {

my $utr_start;
my $utr_end;
my $orf_start;
my $orf_end;
my $orf_overlap;
my $j;
my $i;
my $k;
my $m;
my $ol_buffer = $CDS_threshold;  # the focal seq is compared to all other seqs, ol_buffer is length of longest cds block so far found
			     # provided it is > $CDS_threshold  bp long
my $ol_start;
my $ol_end;
my $overlap;
my $start_j;
my $end_j;

my @frame_buffer;
my $orf_frame = 0;
my $template_frame = 0;
my $template;
my $endpoint;

	my ($seq) = @_;


	$endpoint = 0;
	while ($endpoint < $alignment_length) {            # loop over blocks in target sequence
		for ($k=0;$k<@{$status{$seq}}; $k++) {
			if ($changepos{$seq}[$k] > $endpoint and $status{$seq}[$k] eq 'plain') {   # only  process utr blocks downstream of any previously analsed utr in this sequence
				$orf_start = -1;
				$orf_end = -1;
				$ol_buffer = $CDS_threshold;
				if ($k == 0) { $utr_start = 0 }
				else	     { $utr_start = $changepos{$seq}[$k-1] }
				$utr_end = $changepos{$seq}[$k];
				
				if ($utr_end - $utr_start > $CDS_threshold) {        # ********** no point in dealing with tiny utrs ****************
					# loop over all other sequences
					for ($i=0; $i <  $seq_count; $i++) {
						if ($i != $seq) {
							@frame_buffer = ();
							$orf_overlap = 0;
							for ($j=0; $j < @{$status{$i}}; $j++) {    # loop over all blocks in seq i and record start and endpoints
																	   # of cds overlap with utr in target seq  
								if ($j == 0) { $start_j = 0 }
								else 	     { $start_j = $changepos{$i}[$j-1] }	
								$end_j = $changepos{$i}[$j];
								# overlap may stretch over more than one block in seq i (assumed: overlapping orf blocks in i are contiguous)
								if (!($utr_start >= $end_j or $utr_end <= $start_j) and $patterns{$i}[$j] >= 110) {
									if (!($orf_overlap)) {
										if ($utr_start < $start_j) 	  	{$ol_start = $start_j }
										else 				  	{$ol_start = $utr_start }
										$orf_overlap = 1;
									}
									if ($utr_end > $end_j)      	{$ol_end = $end_j}
									else 				{$ol_end = $utr_end }
									push @frame_buffer, $frames{$i}[$j];
								}
							}
							if ($orf_overlap) { 			# if cds overlap exists, check sequence identity
								$overlap = $ol_end - $ol_start;
								my $matches = 0;
#								print "overlap: $ol_start-$ol_end, recipient: $seq, template: $i\n";
								for ($m=$ol_start; $m < $ol_end; $m++) {
									if (substr($sequences_sorted[$seq],$m,1) eq substr($sequences_sorted[$i],$m,1)) {
										$matches++;
									}
								}
								# best to use stringent criteria here:
								if ($overlap > $ol_buffer and abs(($matches / $overlap) - 1)  < 0.01) {
									$orf_start = $ol_start;
									$orf_end = $ol_end;
									$ol_buffer = $overlap;
									$template = $i;
									if ($seq_stats[$seq]->{orientation} eq '+') 	{ $template_frame = pop(@frame_buffer) }
									else					   	{ $template_frame = shift(@frame_buffer) }
								}		
							}
						}
					}  # loop over sequences
					if ($orf_start > -1) {
#						print "processing orf overlap from $template to $seq, $orf_start - $orf_end\n";
						$endpoint = $utr_end;			# about to exit for-loop, so endpoint needs to be reset
						$orf_frame = get_correct_orf_frame ($seq,$template,$template_frame,\$orf_start,\$orf_end);
						insert_feature ('cporf',$orf_start,$orf_end,$seq,$k,$orf_frame);
						check_frame_of_obscured_orf ($seq,$orf_end,$orf_frame);
						last;
					}
				}  # if utr long enough
			} # if utr
			if ($k == @{$status{$seq}} - 1) { $endpoint = $alignment_length }
		}	# loop over blocks in target sequence	
	}
	
}

sub get_correct_orf_frame {

my $i;
my $basecount;
my $startpos;
my $endpos;
my $framepos;


	my ($seq,$t,$t_frame,$orfstart_ref,$orfend_ref) = @_;

	#  determine the positions of the start and end nucleotides of the orfblock in the local coordinates of the template sequence and in the correct direction
	#  then determine the start and end positions of the first and last codon, respectively, within that window (according to the template frame)
	#  reset $$orfstart_ref and $$orfend_ref accordingly
	#  finally determine the frame of these codons in the recipient sequence

	# do all of the above in forward and, alternatively, in reverse orientation

	# NB: $startpos and $endpos are defined in the orf orientation. I.e. for minus orfs $startpos is on the right and
	# $endpos is on the left. This means that $template_frame can be applied uniformly to plus and minus orfs.
	# Only the adjustment of the orfrange has be done differely for plus and minus sequences ($$orfstart_ref and
	# $$orfend_ref).

	if ($seq_stats[$t]->{orientation} eq '+') {
		$basecount = 0;
		for ($i=0;$i<=$$orfend_ref;$i++) {            #  note that $$orfend_ref is the first base *after* the orf block
			if (substr($sequences_sorted[$t],$i,1) =~ /[ACTG]/) {	$basecount += 1 }
			if ($i == $$orfstart_ref ) { $startpos = $basecount }
		}
		# count the 1. position past the end of the orf block even if it is a gap
		if (substr($sequences_sorted[$t],$$orfend_ref,1) =~ /-/) { $basecount += 1 }
		$endpos = $basecount;		# if the last codon in the orfblock ends immediately before $$orfend_ref, then $basecount (= $endpos) is at pos 1 in frame (see below)
	}
	else {
		$basecount = 0;
		for ($i=$alignment_length-1;$i>=$$orfstart_ref;$i--) {
			if (substr($sequences_sorted[$t],$i,1) =~ /[ACTG]/) {	$basecount += 1 }
			if ($i == $$orfend_ref - 1) { $startpos = $basecount }
		}
		$endpos = $basecount + 1;    # move one position beyond end, so that the modulo method works
	}

	# Start- and endpos are based on basecount (starting at 1) and not on coordinates. 
	# So, the modulo method works.

	
	for ($i=$startpos;$i<=$startpos+2;$i++) {
		$framepos = $i % 3;
		if ($framepos == 0) { $framepos = 3 }
		if ($framepos == abs($t_frame)) { 
			if  ($seq_stats[$seq]->{orientation} eq '+') {
				$$orfstart_ref = $$orfstart_ref + ($i - $startpos);
			}
			else {
				$$orfend_ref = $$orfend_ref - ($i - $startpos);
			}
		}
	}
	for ($i=$endpos;$i>=$endpos-2;$i--) {
		$framepos = $i % 3;
		if ($framepos == 0) { $framepos = 3}
		if ($framepos == abs($t_frame)) {
			if  ($seq_stats[$seq]->{orientation} eq '+') {
				$$orfend_ref = $$orfend_ref - ($endpos - $i);
			}
			else {
				$$orfstart_ref = $$orfstart_ref + ($endpos - $i);
			}
		}
	}

	# frame of $$orfstart_ref in recipient sequence

	if ($seq_stats[$seq]->{orientation} eq '+') {
		$basecount = 0;
		for ($i=0;$i<=$$orfstart_ref;$i++) {
			if (substr($sequences_sorted[$seq],$i,1) =~ /[ACTG]/) {	$basecount += 1 }
		}
		$framepos = $basecount % 3;
		if ($framepos == 0) { $framepos = 3}
		return $framepos;
	}
	else {
		$basecount = 0;
		for ($i=$alignment_length-1;$i>=$$orfend_ref-1;$i--) {
			if (substr($sequences_sorted[$seq],$i,1) =~ /[ACTG]/) {	$basecount += 1 }
		}
		$framepos = $basecount % 3;
		if ($framepos == 0) { $framepos = 3}
		return $framepos * (-1);
	}

}		
	

sub check_frame_of_obscured_orf {

my $i;
my $frame_buffer = 0;
my $block_found = 0;

	my ($seq,$end,$frame) = @_;


	# Frameshifts among 'regular' cds blocks have already been noted
	# Here, one needs to just look at the frames of the cds blocks immediately to the left and right of the 
	# inferred orf block. There are three cases to consider...
	# 1. current block is the inferred orf block: compare to neighbouring cds on the left (in $frame_buffer)
	# 2. current block is the next one after the inferred block, compare frames and exit loop
	# 3. current block is cds to the left of the inferred block: store frame in $frame_buffer
	for ($i=0;$i<@{$frames{$seq}};$i++) {
		if ($changepos{$seq}[$i] == $end) {
			$block_found = 1;
			if ($frame_buffer != 0 and $frame_buffer != $frame) {
				$seq_stats[$seq]->{frameshift} += 1;
			}
		}
		elsif ($block_found and $status{$seq}[$i] eq 'cds') {
			if ($frames{$seq}[$i] != $frame) { $seq_stats[$seq]->{frameshift} += 1 }
			last;
		}
		elsif (!($block_found) and $status{$seq}[$i] eq 'cds') {  $frame_buffer = $frames{$seq}[$i] }
	}
}

sub define_orf_range {

my $i;
my $orf_start = -1;
my $orf_end = 0;

my $codon;

	my ($seq) = @_;
	
	#  Irrespective of the direction of the orf, orf_start is always on the left and orf_end on the right
	#  Direction is important, clearly, for identifying start and stop codons and 5' vs 3' utrs

		for ($i=0;$i<@{$patterns{$seq}};$i++) {
			if ($patterns{$seq}[$i] >=  110 or $status{$seq}[$i] eq 'cporf') {
				if ($orf_start == -1 ) {
					if ($i == 0) 	{ $orf_start = 0 }
					else 		{ $orf_start = $changepos{$seq}[$i-1] }
					$seq_stats[$seq]->{orf_start} = $orf_start;
				}
				if ($changepos{$seq}[$i] > $orf_end) { 
						$orf_end = $changepos{$seq}[$i];
						$seq_stats[$seq]->{orf_end} = $orf_end;
				}
			}
		}
		
		if ($seq_stats[$seq]->{orientation} eq '+') {
			if ($orf_start != -1) {					# start codon on the left hand side
				if (substr($sequences_sorted[$seq],$orf_start,3) =~ /ATG/) { $seq_stats[$seq]->{start_codon} = 1 }
			}
			if ($orf_end > 0) {						# stop codon on the right hand side
				foreach $codon (@forward_stops) {
					if (substr($sequences_sorted[$seq],$orf_end-3,3) eq $codon) {
						$seq_stats[$seq]->{stop_codon} = 1;
						last;
					}
				}
			}
			for ($i=0;$i<@{$patterns{$seq}};$i++) {
				if ($changepos{$seq}[$i] <= $orf_start and $status{$seq}[$i] eq 'plain' and $seq_stats[$seq]->{start_codon} ) {
					$status{$seq}[$i] = 'utr';
				}
				elsif ($changepos{$seq}[$i] > $orf_end and $status{$seq}[$i] eq 'plain' and $seq_stats[$seq]->{stop_codon} ) {
					$status{$seq}[$i] = 'utr';
				}
			}
		}	
		else {
			if ($orf_end > 0) {						# start codon on the right hand side
					if (substr($sequences_sorted[$seq],$orf_end-3,3) =~ /CAT/) { $seq_stats[$seq]->{start_codon} = 1 }
			}
			if ($orf_start != -1) {					# stop codon on the left hand side
				foreach $codon (@reverse_stops) {
					if (substr($sequences_sorted[$seq],$orf_start,3) eq $codon) {
						$seq_stats[$seq]->{stop_codon} = 1;
						last;
					}
				}
			}
			for ($i=0;$i<@{$patterns{$seq}};$i++) {
				if ($changepos{$seq}[$i] <= $orf_start and $status{$seq}[$i] eq 'plain' and $seq_stats[$seq]->{stop_codon} ) {
					$status{$seq}[$i] = 'utr';
				}
				elsif ($changepos{$seq}[$i] > $orf_end and $status{$seq}[$i] eq 'plain' and $seq_stats[$seq]->{start_codon} ) {
					$status{$seq}[$i] = 'utr';
				}
			}
		}
			
}

sub check_for_internal_stops {

my $i;
my $end;
my $pos;

	my ($seq) = @_;

	for ($i=0;$i<@{$patterns{$seq}};$i++) {
		if ($status{$seq}[$i] eq 'cds' or $status{$seq}[$i] eq 'cporf') {
			if ($i == 0) 	{ $pos = 0 }
			else			{ $pos = $changepos{$seq}[$i-1] }
			# this will count all stop codons in cds/cporf that are not that last codon of the last orf
			if ($changepos{$seq}[$i] == $seq_stats[$seq]->{orf_end}) {
				$end = $changepos{$seq}[$i] - 6;
			}
			else {
				$end = $changepos{$seq}[$i] - 3;
			}
			while ($pos <= $end) {
				if ($seq_stats[$seq]->{orientation} eq '+') {
					foreach $stop (@forward_stops) {
						if (substr($sequences_sorted[$seq],$pos,3) eq $stop) { $seq_stats[$seq]->{int_stop}++ }
					}
				}
				else {
					foreach $stop (@reverse_stops) {
						if (substr($sequences_sorted[$seq],$pos,3) eq $stop) { $seq_stats[$seq]->{int_stop}++ }
					}
				}
				$pos += 3;
			}
		}
	}
}
	
	
			
sub check_for_intron {

#	my $i;
	my $i_start;
	my $i_end;
	my $teststring;
	my $intron_length = 0;
	my $fiveprimesignal;
	my $threeprime_A_pos;
	my $Ynumber;
	my $search_length;   # the length of the teststring in which any 3' signal must be located
	my $intron_found = -1;
	my $intron_seq;
	my $search_start;
	my $search_end;
	my $window_length;
	my $window_start;
	my $threeprimeseq;

	my $max_length = 0;    # This is used only for the global intron search

#	my ($seq,$k) = @_;
	my ($i) = @_;

	# the parameters ensure that gaps that are too small for an intron will be left unchecked here
#	for ($i=0;$i<$seq_count;$i++) {
#		if ($i != $seq) {
#			if ($seq_stats[$seq]->{orientation} eq '+') {
			if ($seq_stats[$i]->{orientation} eq '+') {
#				$search_start = $changepos{$seq}[$k-1] - $I_PAD;
#				$search_end = $changepos{$seq}[$k] + $I_PAD - $MIN_INTERVAL - $MIN_LENGTH_3P;      # the last position where an intron could start										 # and still fit into the search window

				$search_start = 0;
				$search_end = $alignment_length - $MIN_INTERVAL - $MIN_LENGTH_3P;


				for ($j = $search_start; $j < $search_end; $j++) {
					if (substr($sequences_sorted[$i],$j,2) =~ /GT/) {
#						$search_length = ($changepos{$seq}[$k] + $I_PAD) - ($j + $MIN_INTERVAL);
						$search_length = $alignment_length - ($j + $MIN_INTERVAL);
						$teststring = substr($sequences_sorted[$i],$j + $MIN_INTERVAL,$search_length);
						if ($teststring =~ /(A[CT]{2}[CTG][CT][CTG][CT]{10,20}[ACTG][CT]AG)/) {
							if ($MIN_INTERVAL + length ($`) + length ($1) > $intron_length) {
								$intron_length = $MIN_INTERVAL + length ($`) + length ($1);
#								if ($intron_length > $max_length) { $max_length = $intron_length }
								$i_start = $j;
								$i_end = $j + $intron_length + 1;
								$fiveprimesignal = substr($sequences_sorted[$i],$j,6);
								$threeprime_A_pos = $MIN_INTERVAL + length ($`);
								$Ynumber = length($1) - 10; 	# only the pure polypyrmidine stretch is considered here
#								$intron_seq = $i;
			print $iinfo "$seqnames_sorted[$i]\t+\t$i_start\t$i_end\t$fiveprimesignal\t$threeprime_A_pos\t$Ynumber\t$1\n";
							}  # if intron_length longer than previous
						}  # if 3' splice_site found
					} # if 5' GT found
				} # loop over search window
			}
			else {
#				$search_start = $changepos{$seq}[$k] + $I_PAD;
#				$search_end = $changepos{$seq}[$k-1] - $I_PAD + $MIN_INTERVAL + $MIN_LENGTH_3P;     

				$search_start = $alignment_length - 1;
				$search_end = $MIN_INTERVAL + $MIN_LENGTH_3P;
				for ($j = $search_start; $j >= $search_end; $j--) {
					if (substr($sequences_sorted[$i],$j-1,2) =~ /AC/) {
#						$window_start = $changepos{$seq}[$k-1] - $I_PAD;
#						$search_length =  ($j - $MIN_INTERVAL) - $window_start;

						$window_start = 0;
						$search_length = $j - $MIN_INTERVAL; 

						$teststring = substr($sequences_sorted[$i],$window_start,$search_length);
						if ($teststring =~ /(CT[GA][ACTG][GA]{10,20}[GAC][GA][GAC][GA]{2}T)/) {
							if ($MIN_INTERVAL + length ($') + length ($1) > $intron_length) {
								$intron_length = $MIN_INTERVAL + length ($') + length ($1);
#								if ($intron_length > $max_length) { $max_length = $intron_length }
								$i_start = $window_start + length ($`);
								$i_end = $j;
								$fiveprimesignal = substr($sequences_sorted[$i],$j-5,6);
								$fiveprimesignal = reverse_complement ($fiveprimesignal);
								$threeprime_A_pos = $j - $MIN_INTERVAL - length ($');
								$Ynumber = length($1) - 10; 	# only the pure polypyrmidine stretch is considered here
#								$intron_seq = $i;
								$threeprimeseq = reverse_complement ($1);
			print $iinfo "$seqnames_sorted[$i]\t-\t$i_start\t$i_end\t$fiveprimesignal\t$threeprime_A_pos\t$Ynumber\t$threeprimeseq\n";
							}  # if intron_length longer than previous
						}  # if 3' splice_site found
					} # if 5' GT found
				} # loop over search window
			}
#		} # if i is not gap seq
#	}  # loop over sequences
	if ($intron_length > 0) {
		$seq_stats[$i]->{intron_bases} = $intron_length;
		$seq_stats[$i]->{intron_count} += 1;
#		print $iinfo "$seqnames_sorted[$intron_seq]\t$$i_start\t$fiveprimesignal\t$threeprime_A_pos\t$Ynumber\n";
#		insert_intron_annotation ($intron_seq,$i_start,$i_end);
	}
}


sub reverse_complement {

my $i;
my $revcomp = '';
my $base;
my %comphash = (
	A => 'T', 
	C => 'G', 
	G => 'C', 	
	T => 'A',
); 

	my ($seqstring) = @_;

	for ($i=length($seqstring)-1;$i>=0;$i--) {
		$base = substr($seqstring,$i,1);
		if (exists $comphash{$base}) 	{$revcomp .= $comphash{$base} }
		else				{ print "reverse_complement: comphash for key $base not found, $component\n" }
	}

	return $revcomp;
}

sub insert_intron_annotation {

my $i;
my $j;
my $k;
my $intron_length;
my $overlap;
my $matches;
my $blockstart;
my $blockend;
my $ol_start;
my $ol_end;

	my ($intron_seq,$i_start,$i_end) = @_;

	# loop over all sequences and and over all blocks to apply info of the identified intron
	# three different situation:
	# 1. sequence is intro_seq -> set all overlapping block segments to intron
	# 2. sequence block is gap -> set to igap (this case includes the sequence that triggered the intron search)
	# 3. sequence block is not gap -> add to cumulative count of bases that overlap with intron and of matching bases with intron_seq
	# 	if sequence has complete intron with > 0.98 id -> set all overlapping block segments to 'intron'

	$intron_length = $i_end - $i_start;

	for ($i=0;$i<$seq_count;$i++) {
		$overlap = 0;
		$matches = 0;
		for ($j=0;$j<@{$patterns{$i}};$j++) {
			if ($j == 0) 	{ $blockstart = 0 }
			else 		{ $blockstart = $changepos{$i}[$j-1] }
			$blockend = $changepos{$i}[$j];
			if (!($i_start >= $blockend or $i_end <= $blockstart)) {
				if ($i == $intron_seq) { insert_feature ('intron',$i_start,$i_end,$i,$j,0) }   		# CASE 1
				elsif ($status{$i}[$j] == 'gap') { insert_feature ('igap',$i_start,$i_end,$i,$j,0) }  # CASE 2
				else {																				# CASE 3
					if ($i_start > $blockstart) 	{ $ol_start = $i_start }
					else  				{ $ol_start = $blockstart }
					if ($i_end < $blockend)		{ $ol_end = $i_end }
					else				{ $ol_end = $blockend }
					for ($k=$ol_start;$k<$ol_end;$k++) {  
						if (substr($sequences_sorted[$intron_seq],$k,1) eq substr($sequences_sorted[$i],$k,1)) { $matches++ }
					}
					$overlap += $ol_start - $ol_end;
				}
			}
		}
		if ($overlap == $intron_length and abs($matches/$overlap - 1) < 0.02) {    # CASE 3
			for ($j=0;$j<@{$patterns{$i}};$j++) {
				if ($j == 0) 	{ $blockstart = 0 }
				else 		{ $blockstart = $changepos{$i}[$j-1] }
				$blockend = $changepos{$i}[$j];
				if (!($i_start >= $blockend or $i_end <= $blockstart)) {
					insert_feature ('intron',$i_start,$i_end,$i,$j,0);
				}
			}
		}
	}
}

sub check_for_splice {

my $splice_found = 0;
my $i;
my $j;
my $gapstart;
my $gapend;
my $blockstart;
my $blockend;
my $remainder;

	my ($seq,$k) = @_;
	
	#  a splice gap is identified by the following criteria:
	#  1. cds adjoining at least on one side of the gap
	#  2. at least one overlapping block in at least one other sequence has status "cds"
	
	$gapstart = $changepos{$seq}[$k-1];
	$gapend = $changepos{$seq}[$k];
	if (($status{$seq}[$k-1] eq 'cds' or $status{$seq}[$k+1] eq 'cds' or $status{$seq}[$k-1] eq 'cporf' or 
		$status{$seq}[$k+1] eq 'cporf') and $gapend - $gapstart > $MIN_EXON_LENGTH) {
		for ($i=0;$i<$seq_count;$i++) {
			if (!($i == $seq)) {
				for ($j=0;$j<@{$patterns{$i}};$j++) {
					if ($j == 0) 	{ $blockstart = 0 }
					else			{ $blockstart = $changepos{$i}[$j-1] }
					$blockend = $changepos{$i}[$j];
					if (!($blockend <= $gapstart or $blockstart >= $gapend) and ($status{$i}[$j] eq 'cds' or $status{$i}[$j] eq 'cporf')) {
						if (($gapend - $gapstart) % 3 != 0) {
							if ($seq_stats[$seq]->{frameshift} > 0) { $splice_found = 1 }
						}
						else { $splice_found = 1 }
						last;
					}
				}
				if ($splice_found) { last }
			}
		}
	}
	if ($splice_found) {
		insert_feature ('splice',$gapstart,$gapend,$seq,$k,0);
	}
	else {
		insert_feature ('indel',$gapstart,$gapend,$seq,$k,0);
	}
}	
				
					
sub insert_feature {

	# this sub has two functions: 1. to insert annotation patterns into the basic sequence alignment ($feature = numeric code to be inserted into @patterns)
	# or 2. to add @status info inferred from comparison with other sequences ($feature is a string such as 'cds', 'intron' etc.)
 
 	my $insert_pos;
	my $lower_bound;
	my $upper_bound;
	my $current_pattern;
	my $current_status;
	my $current_frame;
	my $blocks_inserted = 0;
	my $new_pattern;
	my $new_frame;
	my $endbuffer;
 
	my ($feature, $fstart, $fend, $seq, $block,$frame) = @_;

#	print "insert_feature: $component, seq $seq, block $block, frame $frame, current frame: $frames{$seq}[$block]\n";


	$current_pattern = $patterns{$seq}[$block];
	$current_frame = $frames{$seq}[$block];

	if ($feature =~ /\d+/)  {$new_pattern = $feature }     # if the update pertains to @status info, the @patterns info does not change
	else			{$new_pattern = $current_pattern }


	if ($frame != 0) {
		if ($current_frame != 0 and $current_frame != $frame) {
			$new_frame = resolve_frame ($seq,$block,$fstart,$fend,$frame);
		}
		else {
			$new_frame = $frame;
		}
	}
	else { $new_frame = $current_frame }


	if ($block == 0) { $lower_bound = 0 }
	else 		 { $lower_bound = $changepos{$seq}[$block-1] }


	# feature overlaps with or starts at lower bound of block and ends within the block
	if ($fstart <= $lower_bound and $fend < $changepos{$seq}[$block]) { 
		$insert_pos = $block;
		splice @{$changepos{$seq}}, $insert_pos, 0, $fend;
		splice @{$patterns{$seq}}, $insert_pos, 0, $new_pattern;
		splice @{$frames{$seq}}, $insert_pos,0, $new_frame;
		if (!($feature =~ /\d+/)) {
			splice @{$status{$seq}}, $insert_pos, 0, $feature;
		}
		$blocks_inserted = 1;
#		print "1 block inserted, $changepos{$i}[$j-1] - $fend\n";
	}
	
	# feature extends past block on either side
	elsif ($fstart <= $lower_bound and $fend >= $changepos{$seq}[$block]) {
		if ($feature =~ /\d+/) 	{ 
			if ($patterns{$seq}[$block] > 0 ) { 
				$patterns{$seq}[$block] = $new_pattern;
				$frames{$seq}[$block] = $new_frame;
			}
		}	
		else 	{ $status{$seq}[$block] = $feature }
#		print "0 block inserted, pattern change only\n";
	}

	# feature start within block and either ends at blockend or extends beyond it
	elsif ($fstart > $lower_bound and $fend >= $changepos{$seq}[$block]) {
		$endbuffer = $changepos{$seq}[$block];
		$changepos{$seq}[$block] = $fstart;
		$insert_pos = $block+1;                  # if $block is last element in current array then this will add an element to the array
		splice @{$changepos{$seq}}, $insert_pos, 0, $endbuffer;
		splice @{$patterns{$seq}}, $insert_pos, 0, $new_pattern;
		splice @{$frames{$seq}}, $insert_pos, 0, $new_frame;
		if (!($feature =~ /\d+/)) {
			splice @{$status{$seq}}, $insert_pos, 0, $feature;
		}
		$blocks_inserted = 1;
#		print "1 block inserted, $changepos{$seq}[$block] - $endbuffer\n";

	}

	# feature is contained within block
	elsif ($fstart > $lower_bound and $fend < $changepos{$seq}[$block]) {
		$upper_bound = $changepos{$seq}[$block];
		$changepos{$seq}[$block] = $fstart;
		$insert_pos = $block+1;
		splice @{$changepos{$seq}}, $insert_pos, 0, $fend, $upper_bound;
		splice @{$patterns{$seq}}, $insert_pos, 0, $new_pattern, $current_pattern;
		splice @{$frames{$seq}}, $insert_pos, 0, $new_frame, $current_frame;
		if (!($feature =~ /\d+/)) {
			$current_status = $status{$seq}[$block];
			splice @{$status{$seq}}, $insert_pos, 0, $feature, $current_status;
		}
		$blocks_inserted = 2;
#		print "2 blocks inserted: endpoints $fstart - $fend - $upper_bound\n";
#		print "2 blocks inserted: patterns $current_pattern - $feature - $current_pattern\n";
#		print "insert_feature: total blocks now: ",scalar @{$patterns{$seq}},"\n";
	}
	return $blocks_inserted;
}


sub resolve_frame {

	my $correct_frame;
	my $current_frame;
	my $current_start;
	my $j;
	my $blockstart;
	my $blockend;
	my $ol_start;
	my $ol_end;
	my $overlap;

	my ($seq,$block,$fstart,$fend,$frame) = @_; 

	$current_frame = $frames{$seq}[$block];

	if (($current_frame > 0 and $frame < 0) or ($current_frame < 0 and $frame > 0)) {           	#  This won't happen any more because all annotations now have the same orientation
		print "$component - seq $seq, block $block, current frame: $current_frame, new frame: $frame\n";
		$correct_frame = 0;
	}
	elsif ($current_frame > 0) {									
		# both annotations in positive direction: correct frame is the one that starts further to the right
		# find the start of the current frame to the left of the current block
		if ($block == 0) { $current_start = 0 }
		else 		 { $current_start = $changepos{$seq}[$block-1] }
		for ($j=$block-1;$j >= 0;$j--) {
			if ($patterns{$seq}[$j] > 100) {  	# consider only blocks that have an annotation (ingnore gaps and non-coding sequence)
				if ($frames{$seq}[$j] == $current_frame) {
					if ($j == 0) 	{ $current_start = 0 }
					else 		{ $current_start = $changepos{$seq}[$j-1] }
				}
				else { last }
			}
		}
		if ($current_start <= $fstart) 	{ $correct_frame = $frame }
		else 				{ $correct_frame = $current_frame }
	}
	else {												
		# both annotations in negative directions: correct frame is the one that starts further to the left
		# find the start of the current frame to the right of the current block
		$current_start = $changepos{$seq}[$block];
		for ($j=$block+1;$j < @{$changepos{$seq}};$j++) {
			if ($patterns{$seq}[$j] > 100) {     # consider only blocks that have an annotation (ingnore gaps and introns)
				if ($frames{$seq}[$j] == $current_frame) {  $current_start = $changepos{$seq}[$j] }
				else { last }
			}
		}
		if ($current_start >= $fend) 	{ $correct_frame = $frame }
		else 				{ $correct_frame = $current_frame }
	}
	
	if ($block == 0) 	{ $blockstart = 0 }
	else			{ $blockstart = $changepos{$seq}[$block-1] }
	if ($blockstart < $fstart) 	{ $ol_start = $fstart }
	else 				{ $ol_start = $blockstart }
	if ($changepos{$seq}[$block] < $fend) 	{ $ol_end = $changepos{$seq}[$block] }
	else 					{ $ol_end = $fend }
	$overlap = $ol_end - $ol_start;
	print $fshift "$component\t$seq\t$block\t$overlap\n";

	return $correct_frame;
}
		
			
			

sub get_per_block_status {

	my $i;
	my $j;
	my $done;
	my $status_a = 'null';
	my $status_b = 'null';

	my ($seqa, $seqb, $pos) = @_;
	$done = 0;
	
	# NB gaps do not need to be considered here because the sequence windows passed to this function are always gapfree
	
	# seqa
	
	if ($pos < $seq_stats[$seqa]->{cds_start}) {
		if ($pos+$AL_WINDOW < $seq_stats[$seqa]->{cds_start})  { $status_a = 'utr' }
	}
	elsif ($pos >= $seq_stats[$seqa]->{cds_end}) {  $status_a = 'utr' }
	else { 
			for ($i=0;$i<@{$patterns{$seqa}};$i++) {
				if ($pos < $changepos{$seqa}[$i]) {
					if (!($status{$seqa}[$i] eq 'plain')) {
						for ($j=$i;$j<@{$patterns{$seqa}};$j++) {
							if ($pos+$AL_WINDOW < $changepos{$seqa}[$j]) {
								if (!($status{$seqa}[$j] eq 'plain') and $pos+$AL_WINDOW < $seq_stats[$seqa]->{cds_end}) {
									$status_a = 'cds';
								}
								last;  # make sure the loop over blocks is exited if the correct block has been found
							}  # 2. correct block
						}  # 2. for
					} # 1. if !plain
					last;   # same here
				} # 1. correct bloc
			} # 1. for
	}
	
	#  seqb
	
	if ($pos < $seq_stats[$seqb]->{cds_start}) {
		if ($pos+$AL_WINDOW < $seq_stats[$seqb]->{cds_start})  { $status_b = 'utr' }
	}
	elsif ($pos >= $seq_stats[$seqb]->{cds_end}) {
		if  ($pos+$AL_WINDOW >= $seq_stats[$seqb]->{cds_end}) {  $status_b = 'utr' }
	}
	else { 
			for ($i=0;$i<@{$patterns{$seqb}};$i++) {
				if ($pos < $changepos{$seqb}[$i]) {
					if  (!($status{$seqb}[$i] eq 'plain')) {
						for ($j=$i;$j<@{$patterns{$seqb}};$j++) {
							if ($pos+$AL_WINDOW < $changepos{$seqb}[$j]) {
								if (!($status{$seqb}[$j] eq 'plain') and $pos+$AL_WINDOW < $seq_stats[$seqb]->{cds_end}) {
									$status_b = 'cds';
								}
								last;   
							}  # 2. corect block
						}  # 2. for
					}  # 1 if !plain
					last;
				} # 1. correct block
			}  # 1. for
	}
	

	# the first of these cases is the default in case one of the sequences doesn't have any annotation
	if ($seq_stats[$seqa]->{cds_start} == $alignment_length or $seq_stats[$seqb]->{cds_start} == $alignment_length) {
		return 'utr5p';
	}
	elsif ($status_a eq 'cds' and  $status_b eq 'cds')  { 
		return 'cds' 
	}
	elsif ($status_a eq 'utr' and $status_b eq 'utr' ) {
		if ($seq_stats[$seqa]->{orientation} eq '+') {
			if ($pos < $seq_stats[$seqa]->{cds_start}) 	{ return 'utr5p' }
			else 						{ return 'utr3p' }
		}
		else {
			if ($pos > $seq_stats[$seqa]->{cds_end}) 	{ return 'utr5p' }
			else 		   				{ return 'utr3p' }
		}
	}
	elsif ($status_a eq 'cds' and $status_b eq 'utr') {
		if ($pos < $seq_stats[$seqb]->{cds_start}) { 
			if  ($seq_stats[$seqa]->{orientation} eq '+') 	{ return 'mm5p' }
			else 						{ return 'mm3p' }
		}
		else {
			if  ($seq_stats[$seqa]->{orientation} eq '+') 	{ return 'mm3p' }
			else 						{ return 'mm5p' }
		}
	}
	elsif ($status_a eq 'utr' and $status_b eq 'cds') {
		if ($pos < $seq_stats[$seqa]->{cds_start}) { 
			if  ($seq_stats[$seqa]->{orientation} eq '+') 	{ return 'mm5p' }
			else 						{ return 'mm3p' }
		}
		else {
			if  ($seq_stats[$seqa]->{orientation} eq '+') 	{ return 'mm3p' }
			else 						{ return 'mm5p' }
		}
	}
	else { return 'null' }
}


sub compute_id {


	my $j;
	my ($a,$b,$k) = @_;
	my $matches = 0;

	for ($j=$k;$j<$k+$AL_WINDOW;$j++) {
		if (substr($sequences_sorted[$a],$j,1) eq substr($sequences_sorted[$b],$j,1)) {
			$matches++;
		}
	}  
	return $matches / $AL_WINDOW;
}


sub record_alignment_data {
	
my $var;

	my ($seq,$worst,$minref,$maxref,$countref,$meanref, $cds_one) = @_;

	$seq_stats[$seq]->{worstMatch} = $worst;
	
	# keys in the input hashes are: utr5p, mm5p, cds, mm3p, utr3p, 
	# each of these is appended with correct extension to form a key in the seq_stats array of hashes

	foreach my $key (keys %{ $minref }) {  # note the code for accessing keys using a hash ref
		$var = $key .'_min';
		$seq_stats[$seq]->{$var} = $minref->{$key};
		$var = $key .'_max';
		$seq_stats[$seq]->{$var} = $maxref->{$key};
		$var = $key .'_aligned';
		$seq_stats[$seq]->{$var} = $countref->{$key};
		$var = $key .'_mean';
		$seq_stats[$seq]->{$var} = $meanref->{$key};
	}
	$seq_stats[$seq]->{cds_100} = $cds_one;
}	

	

sub compact_seq_annotation_arrays {

	my $b;
	my $last_block;
	my ($seq,$option) = @_; 

	$last_block = @{$changepos{$seq}} - 1;
	for ($b = $last_block;$b > 0; $b--) {
		if ($option eq 'no_status') {            # ignore status, because this array hasn't been assigned yet
			if ($patterns{$seq}[$b-1] == $patterns{$seq}[$b]) {
				splice @{$changepos{$seq}}, $b-1, 1;
				splice @{$patterns{$seq}}, $b-1, 1;
				splice @{$frames{$seq}}, $b-1, 1;		
			}
		}
		elsif ($option eq 'status') {		# compact based on status ony, this is the final step before output, where @patterns and @frames lose their meaning
			if ($status{$seq}[$b] eq $status{$seq}[$b-1]) {
				splice @{$changepos{$seq}}, $b-1, 1;
				splice @{$patterns{$seq}}, $b-1, 1;
				splice @{$frames{$seq}}, $b-1, 1;			
				splice @{$status{$seq}}, $b-1, 1;
			}
		}
		else {					# compact based on patterns and status
			if ($patterns{$seq}[$b-1] == $patterns{$seq}[$b] and $status{$seq}[$b] eq $status{$seq}[$b-1]) {
				splice @{$changepos{$seq}}, $b-1, 1;
				splice @{$patterns{$seq}}, $b-1, 1;
				splice @{$status{$seq}}, $b-1, 1;
				splice @{$frames{$seq}}, $b-1, 1;			
			}
		}
	}
}


sub output_annotated_alignment {

my $outname;
my $out;
my $line_length = 80;
my $line;
my $seq0 = '';
my $seq1 = '';
my $xt0 = '';
my $xt1 = '';
my $orf0 = '';
my $orf1 = '';
my $matches = '';
my $m;
my $n;
my $apat;
my $bpat;
my $a;
my $b;
my $i;
my $spaces = 0;
my $pict = 0;
my $frame;
my $direction;
my $findex;
my $align_score = '';

	my ($compname) = @_;


	for ($i=0;$i<$seq_count;$i++) {
		$align_score .= $seq_index[$i];
	}
	
	$outname = $compname . '_' .$align_score . '_aligned.txt';
	open ($out,"> $outname") or die "Can't create $outname\n";

	for ($m=0;$m<$alignment_length;$m++) {
		$a = substr($sequences_sorted[0],$m,1);
		$seq0 .= $a;
		if ($pict) 	{ $pict = 0 }     # toggle between the two components of the frame representation
		else 		{ $pict = 1 }
		for ($n=0;$n<@{$patterns{0}}; $n++) {
			if ($m < $changepos{0}[$n]) {
				$apat = $patterns{0}[$n];
				$frame = $frames{0}[$n]; 
				if ($frame > 0) { $direction = '>' }
				else 		{ $direction = '<' } 
				$findex = abs($frame);
				if ($apat >= 110) { 
					if ($pict) 	{ $orf0 .= $direction }
					else		{ $orf0 .= $findex }
				}
				else				{ $orf0 .= ' ' }
				if ($apat % 10 > 0) { 
					if ($pict) 	{ $xt0 .= $direction }
					else		{ $xt0 .= $findex }
				}
				else				{ $xt0 .= ' ' }
				last;
			}
		}
		$b = substr($sequences_sorted[1],$m,1);
		$seq1 .= $b;
		for ($n=0;$n<@{$patterns{1}}; $n++) {
			if ($m < $changepos{1}[$n]) {
				$bpat = $patterns{1}[$n];
				$frame = $frames{1}[$n]; 
				if ($frame > 0) { $direction = '>' }
				else 		{ $direction = '<' } 
				$findex = abs($frame);
				if ($bpat >= 110) { 
					if ($pict) 	{ $orf1 .= $direction }
					else		{ $orf1 .= $findex }
				}
				else			{ $orf1 .= ' ' }
				if ($bpat % 10 > 0) { 
					if ($pict) 	{ $xt1 .= $direction }
					else		{ $xt1 .= $findex }
				}
				else			{ $xt1 .= ' ' }
				last;
			}
		}
		if ($a =~ /[ACTG]/ and $b =~ /[ACTG]/ and $a eq $b) { $matches .= '*' }
		else						{ $matches .= ' ' }
		
		if (($m+1) % 3 == 0) {
			$seq0 .= ' ';
			$seq1 .= ' ';
			$xt0  .= ' ';
			$xt1  .= ' ';
			$orf0 .= ' ';
			$orf1 .= ' ';
			$matches .= ' ';
			$spaces++;
		}	
	}	
	
	my $line_count = 1;
	$n = 0;
	while ($n < $alignment_length + $spaces) {
		print $out "seq 0, xt\t",substr($xt0,$n,$line_length), "\n";
		print $out "seq 0, orf\t",substr($orf0,$n,$line_length), "\n";		
		print $out "seq 0, seq\t",substr($seq0,$n,$line_length), "\n";
		print $out "       id\t",substr($matches,$n,$line_length)," ",$line_count * 60, "\n";
		print $out "seq 1, seq\t",substr($seq1,$n,$line_length), "\n";
		print $out "seq 1, orf\t",substr($orf1,$n,$line_length), "\n";
		print $out "seq 1, xt\t",substr($xt1,$n,$line_length), "\n";
		print $out "\n";	
		$n += $line_length;
		$line_count++;
	}
	close ($out);
}

sub write_patterns_to_screen {

my $i;
my $j;

my ($seq,$status_flag) = @_;

print "$component\t";
for ($i=0;$i<=$seq;$i++) {

	print "seq $seq_index[$i], no of blocks: ", scalar @{$patterns{$i}},"\n";
	for ($j=0;$j<@{$patterns{$i}};$j++) {
		print "$j\t$changepos{$i}[$j]\t";
	}
	print "\n";

	for ($j=0;$j<@{$patterns{$i}};$j++) {
	print "$j\t$patterns{$i}[$j]\t";
	}
	print "\n";
	if ($status_flag) {
		for ($j=0;$j<@{$patterns{$i}};$j++) {
			print "$j\t$status{$i}[$j]\t";
		}
		print "\n";
	}
	for ($j=0;$j<@{$patterns{$i}};$j++) {
		print "$j\t$frames{$i}[$j]\t";
	}
	print "\n";

}

}

					 	


sub write_output_to_file {		

my $i;
my $print_flag = 0;
my $aligncode = '';

my $outfilename = 'bom_component_crunch_5p3pmm.txt';
open (my $out, ">> $outfilename") or die "Can't append $outfilename\n";


for ($i=0;$i<$seq_count;$i++) {
	$seqnames_sorted[$i] =~ /_seq/;
	$aligncode .= $';
}

if ($comp_count % 10 == 0) {
	print $out ">seq\t";
	print $out "align\tlen\torient\torf#\to_le\txt#\txt_le\tcpo#\tcpo_le\tstart\tstop\tintern\t";
	print $out "utr#\tutr_le\tintr#\tintr_le\tindl\tindl_le\tspl#\tspl_le\tcds#\tcds_le\tplain#\tpln_le\t";
	print $out "utr5p\t\t\t\tmm5p\t\t\t\tcds\tone\t\t\tmm3p\t\t\t\tutr3\t\t\t\tworst\tfshft\tcovr\n";
}

print $out ">$component\n";
for ($i=0;$i<$seq_count;$i++) {
	print $out "seq$seq_index[$i]\t$aligncode\t$seq_stats[$i]->{basecount}\t$seq_stats[$i]->{orientation}\t$seq_stats[$i]->{orfcount}\t";
	print $out "$seq_stats[$i]->{orfbases}\t$seq_stats[$i]->{xtcount}\t$seq_stats[$i]->{xtbases}\t";
	print $out "$seq_stats[$i]->{cporf_count}\t$seq_stats[$i]->{cporf_bases}\t$seq_stats[$i]->{start_codon}\t";
	print $out "$seq_stats[$i]->{stop_codon}\t$seq_stats[$i]->{int_stop}\t";
	print $out "$seq_stats[$i]->{utr_count}\t$seq_stats[$i]->{utr_bases}\t$seq_stats[$i]->{intron_count}\t";
	print $out "$seq_stats[$i]->{intron_bases}\t$seq_stats[$i]->{indel_count}\t$seq_stats[$i]->{indel_bases}\t";
	print $out "$seq_stats[$i]->{splice_count}\t$seq_stats[$i]->{splice_bases}\t$seq_stats[$i]->{cds_count}\t";
	print $out "$seq_stats[$i]->{cds_bases}\t$seq_stats[$i]->{plain_count}\t$seq_stats[$i]->{plain_bases}\t";
	print $out "$seq_stats[$i]->{utr5p_aligned}\t";
	printf $out "%.2f\t%.2f\t%.2f\t",$seq_stats[$i]->{utr5p_min},$seq_stats[$i]->{utr5p_mean},$seq_stats[$i]->{utr5p_max};
	print $out "$seq_stats[$i]->{mm5p_aligned}\t";
	printf $out "%.2f\t%.2f\t%.2f\t",$seq_stats[$i]->{mm5p_min},$seq_stats[$i]->{mm5p_mean},$seq_stats[$i]->{mm5p_max};
	print $out "$seq_stats[$i]->{cds_aligned}\t";
	printf $out "%.2f\t%.2f\t%.2f\t%.2f\t",$seq_stats[$i]->{cds_100},$seq_stats[$i]->{cds_min},$seq_stats[$i]->{cds_mean},$seq_stats[$i]->{cds_max};
	print $out "$seq_stats[$i]->{mm3p_aligned}\t";
	printf $out "%.2f\t%.2f\t%.2f\t",$seq_stats[$i]->{mm3p_min},$seq_stats[$i]->{mm3p_mean},$seq_stats[$i]->{mm3p_max};
	print $out "$seq_stats[$i]->{utr3p_aligned}\t";
	printf $out "%.2f\t%.2f\t%.2f\t",$seq_stats[$i]->{utr3p_min},$seq_stats[$i]->{utr3p_mean},$seq_stats[$i]->{utr3p_max};
	print $out "$seq_stats[$i]->{worstMatch}\t$seq_stats[$i]->{frameshift}\t";
	printf $out "%.2f\n",$coverage{$seq_index[$i]};	
#	printf  "$i\t%.2f\n",$i,$coverage{$seq_index[$i]};	

}


close ($out);

}



