#!/usr/bin/perl
# Trinity_align_components.pl
# all input files need to be specified in the code with paths etc (see below), supply the name of the target component on the command line
# this is a slightly modified version of Trinity_align_componets.pl
# it uses a different format for input parameters: instead of just the component name, you need to supply the component name followed by a list 
# of numbers that are identify the sequences to be used (i.e. only the ones with read support)

# NB: Annotations are written to the .afa2 file in local sequence coordinates with first base = pos 0. Coordinates indicate
# the first and the last annotated base, respectively.

use strict; use warnings;

#  USAGE: command line arguments (required in this order): component name, two or more Trinity sequence indices, (optional): after the last seq index 'msf'
# The required minimal command line is used by Trinity_alanyse_comp_batch.pl and generates an .afa alignment
# The addition of 'msf' generates an msf file instead of an afa file. In this case, the sequences appear in the same order as in the input file

# needed: BLASTX tabular output. The following columns are required: fields[0] = seqname, fields[4] = frames (3/0 format), fields[10] = qstart, fields[11] = qend;

 
# important identifiers and their format:
# $component: comp[0-9]+_c[0-9]+
# $seq_name: $component . '_seq' . [0-9]+

my $component_to_follow = '';


my $OVERLAP_THRESHOLD = 0.5;   #  a new blast match for a given sequence will only be output if its overlap with other matches is < this fraction of its length
my $ORF_THRESHOLD = 0.4;       # an orf has to overlap with a blast match for the same sequence by at least this fraction in order to be kept
							   # (and it has to have the same orientation)
my $line;
my $revline;
my @output = ();   # a collection of header and sequence lines that will be written to the output file
my $flag;
my $contigseq;
my %seqlength; 		# records the length of each sequence, needed for determining the frame of orfs in (-) direction
my $buffer;
my @fields;
my $seqname;
my $seqno;	
my $ol_start;
my $ol_end;
my $matches;
my $start;
my $end;
my $frame;
my %xt_start;			# hashes that hold all the xt match data for each sequence in a component
my %xt_end;
my %xt_frame;
my %orf_start;
my %orf_end;
my %orf_frame;

my @original_name;
my @headerline;
my @seqline;
my $i;
my $j;
my $found;
my $in;
my $out;
my $afa2;
my $revin;

my @original_names = ();
my $index;



# NB: sequence indices need to be in ascending though not necessarily consecutive order

die "Input component and original sequence names, please.\n" if (@ARGV == 0);
my $component = shift (@ARGV);
push @original_name, shift (@ARGV);
push @original_name, shift (@ARGV);
my $reverse = shift (@ARGV);


if ($component eq $component_to_follow) {
	print "align_components: parameters read in, current component = $component. Number of sequences: 2 \n";
}


# **************** SPECIFY INPUT FILES HERE *************************

my $trinity_file = '/users/bnurnberger/Trinity_components/var/trin454/var_nmt_crunch.fasta';
open (my $trin, "< $trinity_file") or die "Can't open $trinity_file.\n";

# orf finder: TransDecoder in Trinity using the script 'transcripts_to_best_scoring_ORFs.pl'
my $orffile = '/users/bnurnberger/Trinity_components/var/trin454/var_nmt_crunch_orf.cds';
open (my $orf, "< $orffile") or die "Can't open $orffile.\n";

# exp blast file columns:
# qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore

my $blastfile = '/users/bnurnberger/Trinity_components/var/trin454/var_nmt_crunch_blastx.out';
open (my $blast, "< $blastfile") or die "Can't open $blastfile.\n";

for ($i=0;$i<2;$i++) {
	@{$xt_start{$i}} = ();
	@{$xt_end{$i}} = ();
	@{$xt_frame{$i}} = ();
	@{$orf_start{$i}} = ();
	@{$orf_end{$i}} = ();
	@{$orf_frame{$i}} = ();

	if ($component eq $component_to_follow) {
		print "align_components: ", $component . '_seq' . $i . "\n";
	}

}

# *********************  collect contigs from Trinity_file *****************************************

$flag = 0;
$contigseq = '';
my $found_count = 0;


# print "looking now for $original_name[0] and $original_name[1]\n";
while ($line = <$trin>) {
	chomp ($line);
	if ($flag) {                                            # flag=1: sequence to be output found
		if ($line =~ m/>/) {				# start of new contig: add sequence to output array and reset flag 
			push (@seqline,$contigseq);
			$seqlength{$i} = length($contigseq);
			$contigseq = '';
			$flag = 0;
			$found_count++;
		}
		else 	{
			$contigseq .= uc($line);
		}			
	}
	if ($line =~ />/) {                              	# check whether header line contains one of the sought-after sequence names
		for ($i=0;$i<2;$i++) {			# if so that flag=1 so that the next lines of sequence are stored in $contigseq
			if ($line =~ /(>$original_name[$i])$/) {	 # this regex is for 454 data: line ends at end of seqname
				push @headerline, $1;
				$flag = 1;
				last;
			}
			if ($line =~ /(>$original_name[$i]) /) {	# this regex is for trinity data: there is a space after the sequence name
				push @headerline, $1;
				$flag = 1;
				last;
			}
		}
	}
	if ($found_count == 2) {
		last;
	}
}

# print "$found_count sequences found\n";

if ($found_count == 1) {		
	push (@seqline,$contigseq);
	$seqlength{$found_count} = length($contigseq);   # $found_count holds both the current array length and the index of the last, missing element
}

close ($trin);
my $dummy;
if ($reverse) {
	open ($out,">contig.fasta") or die "Can't create fasta file\n";
	print $out "$headerline[0]\n$seqline[0]\n";
	close ($out);
	system ("revseq contig.fasta contig.rev 1> /dev/null 2> /dev/null"); 	 
	open($revin, "<contig.rev") or die "Can't open revseq file\n";
	$dummy = <$revin>;	
	$seqline[0] = '';
	while ($line = <$revin>) { 
		chomp($line);
		$seqline[0] .= $line;
	 }
	close ($revin);
	system ("rm contig.fasta");
	system ("rm contig.rev");
}


open ($out,">$component.fasta") or die "Can't create fasta file\n";
if ($reverse) { print $out "$headerline[0]_rev\n$seqline[0]\n$headerline[1]\n$seqline[1]" }
else 	      { print $out "$headerline[0]\n$seqline[0]\n$headerline[1]\n$seqline[1]" }
close ($out);

system ("muscle -in $component.fasta -fastaout $component.afa2 1> /dev/null 2> /dev/null");

# print "returned from muscle after aligning $component\n";
system ("rm $component.fasta");
# system ("ls $component.afa2\n");
print "now defragmenting $component.afa2...\n";
system ("perl /users/bnurnberger/Trinity_components/Trinity_defragment_alignments_2seq.pl $component.afa2");

$headerline[0] =~ s/>//;   # remove > from headerline, which is now the array of sought-after sequences to get annotation
$headerline[1] =~ s/>//;

# print "header0: $headerline[0], header 2: $headerline[1]\n";

open ($afa2, ">>$component.afa2") or die "Can't append $component.afa2\n";

# ***************************collect blast hits ****************************************************** 

$buffer = '';
my $comp;
my $subc;
my @drop_array_elements;
my $match_no;
my @matched_nucleotides;
my $matchlength;
my @startbuffer = ();
my @endbuffer = ();
my @frame_buffer = ();
my $k;

$headerline[1] =~ /comp(\d+)_c(\d+)/;    # headerline[1] always holds a trinity sequence 
my $ref_comp = $1;
my $ref_subc = $2;
my $seq1_hit = 0;	# records whether >= 1 match to seq 1 has been found (important for storing match data from the buffers in the correct hashes)
# print "ref component is $ref_comp\n";
 
while ($line = <$blast>) {
	chomp ($line);
	@fields = split ('\t',$line);
#	if ($fields[0] =~ /comp(\d+)_c(\d+)/) {             	# extract the component name and compare to $component:
#		$comp = $1;
#		$subc = $2;
#		if ($comp > $ref_comp) {
#			$component_to_follow =~ /^comp(\d+)/; 
#			if ($ref_comp eq $1) { print "$fields[0] - stopping to look for $component_to_follow\n"; }				
#			last;
#		 }
#		elsif ($comp == $ref_comp and $subc > $ref_subc) {		# if the current component has a larger index number than than $component 
#			$component_to_follow =~ /^comp(\d+)/; 
#			if ($ref_comp eq $1) { print "$fields[0] - stopping to look for $component_to_follow\n"; }				
#			last;
#		}
						                			# then quit parsing the blast file. This saves a lot of time.
			
											# important that the trinity contigs are listed after the 454 in the fasta 
#	}

	$found = 0;

	for ($i=0;$i<2;$i++) {
		if ($fields[0] eq $headerline[$i]) {
			$found = 1;
			if ($i == 1) { $seq1_hit = 1}
		}
	} 
	next unless ($found);

	# proceed if query 0 or 1 have been matched
			
	if ($fields[10] > $fields[11]) {         # determine start and end of the blast match in the original sequence
		$start = $fields[11];
		$end = $fields[10];	
	}
	else {
		$start = $fields[10];
		$end = $fields[11];
	}

	if ($fields[4] =~ /^([\-]*\d)/) {        # matches the frame index before the / either with or without a preceding minus sign, examples: -2/0, 3/0
		$frame = $1;
	}			
	else {
		$frame = 4;			# error flag: non-sensical frame, format error
	}
		
	# This approach assumes that matches to a given query sequence are consecutive in the blast file (which is true)
	# @startbuffer, @endbuffer and @frame_buffer are reset as soon as a blast line with a new sought-after seqname is encountered

	if ($fields[0] eq $buffer){                  # if this is a new match to the buffered sequence, determine its overlap with existing matches
	
		@matched_nucleotides = ();             # need to reset array before pushing another set of zeros for the next match

		$matchlength = $end - $start + 1;      
		for ($j=0;$j<$matchlength;$j++) {
			push (@matched_nucleotides,0);
		}

		@drop_array_elements = ();
		$matches = @startbuffer;
		$flag = 0;			

		for ($j=0; $j<$matches; $j++) {						# check for overlap with existing blast matches for this query
			if (!($end < $startbuffer[$j] || $start > $endbuffer[$j])) {

		 		if ($startbuffer[$j] >= $start and $endbuffer[$j] <= $end) {    # new match eclipses existing match
					push (@drop_array_elements, $j);
				}
				else {
					if ($start < $startbuffer[$j])  {$ol_start = $startbuffer[$j]}   # determine start- and endpoints of the overlap
					else 				{$ol_start = $start}
					if ($end < $endbuffer[$j])  	{$ol_end = $end}
					else 				{$ol_end = $endbuffer[$j]}

				
					$ol_start -= $start;			# shift $ol_start and $ol_end down to fit them into @matched_nucleotides   
					$ol_end -= $start;
							
					for ($k=$ol_start;$k<=$ol_end;$k++) {
						$matched_nucleotides[$k] = 1;
					}
				}
			}  # match with overlap
		}  # loop over already stored matches
			

		my $sum = 0;
		for ($j=0;$j<$matchlength;$j++) {
			$sum += $matched_nucleotides[$j];
		}
		if ($sum < $OVERLAP_THRESHOLD * $matchlength) { $flag = 1 }
					
	}   # match to sequence in buffer
	else {						# this condition holds only twice: at the first match to the 0th and the 1st query.
		$flag = 1;
#		print "1. match to new sequence $fields[0]\n";
		if (@startbuffer > 0) {			# this conditions is fulfilled only at the first match to the 1st query
			@{$xt_start{0}} = @startbuffer;	# store the final set of match statistics for the 0th query in three hashes
			@{$xt_end{0}} = @endbuffer;		
			@{$xt_frame{0}} = @frame_buffer;
		}
		$buffer = $fields[0];
		@startbuffer = ();
		@endbuffer = ();
		@frame_buffer = ();
		@drop_array_elements = ();
	}
			
	if ($flag) {   # new match has passed all tests: any eclipsed previous matches are removed 
							# and then the new match is added to output, 
				
		if (@drop_array_elements > 0) {
			for ($j=scalar(@drop_array_elements)-1;$j>=0; $j--) {
				splice (@startbuffer,$drop_array_elements[$j],1); 	# remove start and end pos of the eclipsed match from arrays
				splice (@endbuffer,$drop_array_elements[$j],1);
				splice (@frame_buffer,$drop_array_elements[$j],1);
			}
		}

		push (@startbuffer,$start);
		push (@endbuffer,$end);
		push (@frame_buffer,$frame);

				
	}


}  # loop over blast lines
close ($blast);

# After processing the last sequence with read support, the loop over blast lines may continue
# (-> other sequences w/o read support for the same component), and eventually stops when the next component or eof is reached.
# The three arrays (start,end,frame) then still hold the data from the last sequence to be considered. 
# They are copied here into the appropriate hashes, provided that there have been any matches at all ($eq_buffer > -1):

if ($seq1_hit) 	{ $i = 1 }  # if seq1 has blast hits then any remaining content of the buffers belongs to seq1
else 		{ $i = 0 }  # if seq1 doesn't have any blast hits, then whatever is left in the buffers belongs to seq0

if (@startbuffer > 0) {
	@{$xt_start{$i}} = @startbuffer;
	@{$xt_end{$i}} = @endbuffer;		
	@{$xt_frame{$i}} = @frame_buffer;
}

#for ($i=0;$i<2;$i++) {
#	unless (defined @{$xt_start{$i}}) {
#		@{$xt_start{$i}} = ();
#		@{$xt_end{$i}} = ();		
#		@{$xt_frame{$i}} = ();
#	}
#	print "i = $i, $headerline[$i], number of xt matches: ", scalar(@{$xt_start{$i}}), "\n";
#}

		

# ****************** collect open reading frames *****************************************************


while ($line = <$orf>) {
	chomp ($line);
	if ($line =~ />/) {	
		# the following doesn't work because the components are not in numerical order in the orf file		
#		if ($line =~ /comp(\d+)_c(\d+)_seq/) {             	# extract the component name and compare to $component:
#			$comp = $1;
#			$subc = $2;
#			if ($comp > $ref_comp) { last }
#			elsif ($comp == $ref_comp and $subc > $ref_subc) { last }	
#		}

		for ($i=0;$i<2;$i++) {
			if ($line  =~ /$headerline[$i]:(\d+)-(\d+)/) {
				$start = $1;  		# NB: the orf finder counts residues starting at 1, best to keep it that way until frame is determined
				$end = $2;

				if ($headerline[$i] eq '$component_to_follow') {
					print "orf for $headerline[$i] found\n;"
				}
				if ($line =~ /\+/) {
					$frame = $start % 3;                # if the first nuc of the sequence is the start of the first codon, then the remainder is 1,
										# and so the translation is in the first frame (same coding as the blast output)
					if ($frame == 0) { $frame = 3 }
				}
				else {
					$frame = (($seqlength{$i} - $end + 1) % 3);
					if ($frame == 0) 	{ $frame = -3 }
					else 			{ $frame *= -1 }
				}

				push @{$orf_start{$i}}, $start;
				push @{$orf_end{$i}}, $end;		
				push @{$orf_frame{$i}}, $frame;	
				
			}		
		}
	}
}

close ($orf);

#for ($i=0;$i<2;$i++) {
#	unless (defined @{$orf_start{$i}}) {
#		@{$orf_start{$i}} = ();
#		@{$orf_end{$i}} = ();		
#		@{$orf_frame{$i}} = ();
#	}
#}

if ($reverse) { reverse_annotation () }

my $orientation;
my $plus = 0;
my $minus = 0;
my $overlap;
my $keep;

# First look for overlaps between orf pred and blast matches to decide on the correct orientation. Count the number of overlap bases in each direction...

for ($i=0;$i<2;$i++) { 	
	for ($j=0;$j<@{$xt_start{$i}};$j++) {
		for ($k=0;$k<@{$orf_start{$i}};$k++) {
			$overlap = check_for_matched_overlap ($i,$j,$k);
			if ($overlap > 0) { $plus += $overlap }
			elsif ($overlap < 0) { $minus += abs($overlap) }
		}
	}
}

# ... Then see which of the counts is larger. If they are both 0, then no such overlap has been found and the direction is decided based on a staight count
# of annotated bases over sequences and over all xt/orf. If plus bases are more numerous, then the orientation is '+', otherwise it is '-'.

if ($plus > $minus) { $orientation = '+' }
elsif ($minus > $plus) { $orientation = '-' }
elsif ($minus == 0) {
	$orientation = count_plus_and_minus_bases ();
}
else {  print "$component: orientation unclear - minus bases: $minus, plus bases: $plus\n" }

# print "plusbases: $plus, minusbases: $minus\n";

# $orientation = '-';
# print "orientation: $orientation\n";

for ($i=0;$i<2;$i++) {
	for ($j=0;$j<@{$xt_start{$i}};$j++) {
		$keep = 0;
		if ($xt_frame{$i}[$j] > 0 and $orientation eq '+') { $keep = 1}
		elsif ($xt_frame{$i}[$j] < 0 and $orientation eq '-') { $keep = 1}
		if ($keep) {
			$xt_start{$i}[$j] -= 1; 	# leave the rescaling of annotation to start at base 0 til the very end
			$xt_end{$i}[$j] -= 1;
			print $afa2 ">xt\t$i\t$xt_start{$i}[$j]\t$xt_end{$i}[$j]\t$xt_frame{$i}[$j]\n";
		}
	}
	for ($j=0;$j<@{$orf_start{$i}};$j++) {
		$keep = 0;
		if ($orf_frame{$i}[$j] > 0 and $orientation eq '+') { $keep = 1}
		elsif ($orf_frame{$i}[$j] < 0 and $orientation eq '-') { $keep = 1}
		if ($keep) {
			$orf_start{$i}[$j] -= 1; 	# leave the rescaling of annotation to start at base 0 til the very end
			$orf_end{$i}[$j] -= 1;
			print $afa2 ">orf\t$i\t$orf_start{$i}[$j]\t$orf_end{$i}[$j]\t$orf_frame{$i}[$j]\n";
		}
	}
}

close ($afa2);

#  -------------------------------------------------------- END OF MAIN PROGRAM BLOCK -----------------------------------------------------


sub check_for_matched_overlap {

my $overlap = 0;
my $ol_start;
my $ol_end;

	my ($seq,$xt,$orf) = @_;
	if (!($orf_end{$seq}[$orf] <= $xt_start{$seq}[$xt] || $xt_end{$seq}[$xt] <= $orf_start{$seq}[$orf])) {
		if (($orf_frame{$seq}[$orf] < 0 and $xt_frame{$seq}[$xt] < 0) or ($orf_frame{$seq}[$orf] > 0 and $xt_frame{$seq}[$xt] > 0)) {
			if ($xt_start{$seq}[$xt] < $orf_start{$seq}[$orf])  	{$ol_start = $orf_start{$seq}[$orf] }   # determine start- and endpoints of the overlap
			else 							{$ol_start = $xt_start{$seq}[$xt] }
			if ($xt_end{$seq}[$xt] < $orf_end{$seq}[$orf])  	{$ol_end = $xt_end{$seq}[$xt]}
			else 							{$ol_end = $orf_end{$seq}[$orf]}
			$overlap = $ol_end - $ol_start;
			if ($orf_frame{$seq}[$orf] <  0) {  $overlap *= -1 }
		}
	}
	return $overlap;
}


sub count_plus_and_minus_bases {

my $i;
my $j;
my $sum = 0;

	for ($i=0;$i<2;$i++) {
		for ($j=0;$j<@{$xt_start{$i}};$j++) {
			if ($xt_frame{$i}[$j] > 0) 	{ $sum += $xt_end{$i}[$j] - $xt_start{$i}[$j] }
			else 				{ $sum -= $xt_end{$i}[$j] - $xt_start{$i}[$j] }
		}
		for ($j=0;$j<@{$orf_start{$i}};$j++) {
			if ($orf_frame{$i}[$j] > 0) 	{ $sum += $orf_end{$i}[$j] - $orf_start{$i}[$j] }
			else 				{ $sum -= $orf_end{$i}[$j] - $orf_start{$i}[$j] }
		}
	}	
	if ($sum > 0) 	{ return '+' }
	else 		{ return '-' }
}
	


sub reverse_annotation {

my $new_orientation;
my $j;
my $buffer;

	for ($j=0;$j<@{$xt_start{0}};$j++) {
		$buffer = $xt_start{0}[$j];
		$xt_start{0}[$j] = $seqlength{0} - $xt_end{0}[$j] + 1;
		$xt_end{0}[$j] = $seqlength{0} - $buffer + 1;
		if ($xt_frame{0}[$j] > 0) 	{ $new_orientation = '-' }
		else				{ $new_orientation = '+' }
		if ($new_orientation eq '+') {
			$xt_frame{0}[$j] = $xt_start{0}[$j] % 3;
			if ($xt_frame{0}[$j] == 0) {$xt_frame{0}[$j] = 3}
		}
		else {	 
			$xt_frame{0}[$j] = ($seqlength{0} - $xt_end{0}[$j] + 1) % 3;
			if ($xt_frame{0}[$j] == 0) {$xt_frame{0}[$j] = 3}
			$xt_frame{0}[$j] *= -1;
		}
	}
	for ($j=0;$j<@{$orf_start{0}};$j++) {
		$buffer = $orf_start{0}[$j];
		$orf_start{0}[$j] = $seqlength{0} - $orf_end{0}[$j] + 1;
		$orf_end{0}[$j] = $seqlength{0} - $buffer + 1;
		if ($orf_frame{0}[$j] > 0) 	{ $new_orientation = '-' }
		else				{ $new_orientation = '+' }
		if ($new_orientation eq '+') {
			$orf_frame{0}[$j] = $orf_start{0}[$j] % 3;
			if ($orf_frame{0}[$j] == 0) {$orf_frame{0}[$j] = 3}
		}
		else {	 
			$orf_frame{0}[$j] = ($seqlength{0} - $orf_end{0}[$j] + 1) % 3;
			if ($orf_frame{0}[$j] == 0) {$orf_frame{0}[$j] = 3}
			$orf_frame{0}[$j] *= -1;
		}
	}
}



__END__





sub rearrange_msf_file {

my $i;
my $j;
my $msf;
my $msf2;
my $line;
my $header_end = 0;
my @linebuffer = ();
my $first = 1;

open ($msf,"< $component.msf") or die "Can't open $component.msf\n";
open ($msf2, "> $component.msf2") or die "Can't create $component.msf2\n";

while ($line = <$msf>) {
	chomp($line);

	if ($line =~ /\/\//) {                 # copy lines to output file until '//' is reached
		print $msf2 "$line\n";
		$header_end = 1;
		next;
	}

	if ($header_end == 0) {
		print $msf2 "$line\n";
	}
	elsif ($line =~ /^$/) {			# after that look for empty lines that mark the beginning of a new alignment block
		print $msf2 "$line\n";
		if ($first) { $first = 0 }
		else {
			for ($i=0;$i<@indices;$i++) {
				for ($j=0;$j<@linebuffer;$j++) {   	# first loop looks for the sequence itself, NB: space in regex
#					undef ($1);                    	# make sure $1 holds the seqindex of the current and not of a previous match
					if ($linebuffer[$j] =~ /_seq(\d+) / and $1 == $indices[$i]) {
						print $msf2 "$linebuffer[$j]\n";
					}
				}
				for ($j=0;$j<@linebuffer;$j++) {    # second loop looks for any orf or xt of that sequence, NB: period in regex
					if ($linebuffer[$j] =~ /_seq(\d+)\./ and $1 == $indices[$i]) {
						print $msf2 "$linebuffer[$j]\n";
					}
				}
			}
			@linebuffer = ();
		}
	}
	else {
		push (@linebuffer,$line);
	}
}
	print $msf2 "\n";
	for ($i=0;$i<@indices;$i++) {
		for ($j=0;$j<@linebuffer;$j++) {   # first loop looks for the sequence itself, NB: space in regex
			if ($linebuffer[$j] =~ /_seq(\d+) / and $1 == $indices[$i]) {
				print $msf2 "$linebuffer[$j]\n";
			}
		}

		for ($j=0;$j<@linebuffer;$j++) {    # second loop looks for any orf or xt of that sequence, NB: period in regex
			if ($linebuffer[$j] =~ /_seq(\d+)\./ and $1 == $indices[$i]) {
				print $msf2 "$linebuffer[$j]\n";
			}
		}
	}

close ($msf);
close ($msf2);

}						
				 	
	
	
	
	
	
	
	
	

