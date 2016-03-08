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

my $component_to_follow = 'bla';


my $OVERLAP_THRESHOLD = 0.5;   #  a new blast match for a given sequence will only be output if its overlap with other matches is < this fraction of its length

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

my @headerline;
my @seqline;
my $i;
my $j;
my $found;
my $in;
my $out;
my $afa2;

my @indices;
my $index;
my $msf_flag = 0;
my $component_to_match;
my $outname;

# NB: sequence indices need to be in ascending though not necessarily consecutive order

die "Input component name and sequence numbers\n" if (@ARGV == 0);
my $component = shift (@ARGV);

$component_to_match = $component;
if ($component =~ /^comp/) {
	$component_to_match =~ s/^comp/c/;
	$component_to_match =~ s/_c/_g/;
}


while (@ARGV > 0) {
	if ($ARGV[0] =~ /\d+/) 	{ push (@indices,shift (@ARGV)) }
	else 	{ 
		$msf_flag = 1;
		@ARGV = ();
	}
}


if ($component eq $component_to_follow) {
	print "align_components: parameters read in, current component = $component. Number of sequences: ", scalar(@indices), "\n";
}


# **************** SPECIFY INPUT FILES HERE *************************

#  All input files with new format contig names

my $trinity_file = '/exports/projects/bombina/bombina_bombina/Bb_trinity_140915/Bbom_trinity_140915_clean.fasta';
open (my $trin, "< $trinity_file") or die "Can't open $trinity_file.\n";

# orf finder: TransDecoder in Trinity using the script 'transcripts_to_best_scoring_ORFs.pl'
my $orffile = '/exports/projects/bombina/bombina_bombina/Bb_trinity_140915/transdecoder.tmp.2388/Bbom_trinity_140915_clean.fasta.transdecoder.cds';
open (my $orf, "< $orffile") or die "Can't open $orffile.\n";

# exp blast file columns:
# qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore

my $blastfile = '/exports/projects/bombina/bombina_bombina/Bb_trinity_140915/Bb_blastx_xt.out';
open (my $blast, "< $blastfile") or die "Can't open $blastfile.\n";

for ($i=0;$i<@indices;$i++) {
	@{$xt_start{$indices[$i]}} = ();
	@{$xt_end{$indices[$i]}} = ();
	@{$xt_frame{$indices[$i]}} = ();
	@{$orf_start{$indices[$i]}} = ();
	@{$orf_end{$indices[$i]}} = ();
	@{$orf_frame{$indices[$i]}} = ();

	if ($component eq $component_to_follow) {
		print "align_components: ", $component . '_seq' . $indices[$i], "\n";
	}

}





# *********************  collect contigs from Trinity_file *****************************************

$flag = 0;
$contigseq = '';
my $found_count = 0;

while ($line = <$trin>) {
	chomp ($line);
	if ($flag) {                                            # flag=1: sequence to be output found
		if ($line =~ m/>/) {				# start of new contig: add sequence to output array and reset flag 
			push (@seqline,$contigseq);
			$seqlength{$indices[$i]} = length($contigseq);
			$contigseq = '';
			$flag = 0;
			$found_count++;
		}
		else 	{
			$contigseq .= $line
		}			
	}
	if ($line =~ />/) {                              	# check whether header line contains one of the sought-after sequence names
		for ($i=0;$i<@indices;$i++) {			# if so that flag=1 so that the next lines of sequence are stored in $contigseq
			$seqname = $component_to_match . '_i' . $indices[$i];
			if ($line =~ /(>$seqname) /) {	
				$outname = '>' . $component . '_seq' . $indices[$i];
				push(@headerline,$outname);
				$flag = 1;
				last;
			}
		}
	}
	if ($found_count == @indices) {
		last;
	}
}

if ($found_count == @indices - 1) {		
	push (@seqline,$contigseq);
	$seqlength{$indices[$found_count]} = length($contigseq);   # $found_count holds both the current array length and the index of the last, missing element
}

close ($trin);

# print "$component_to_match, $found_count sequences found\n";

open ($out,">$component.fasta") or die "Can't create fasta file\n";
print $out "$headerline[0]\n$seqline[0]\n$headerline[1]\n$seqline[1]";
close ($out);
	
if (!($msf_flag)) {
	system ("muscle -in $component.fasta -fastaout $component.afa2 1> /dev/null 2> /dev/null");
}

# print "returned from muscle after aligning $component\n";
system ("rm $component.fasta");
# system ("ls $component.afa2\n");
print "now defragmenting $component.afa2...\n";
system ("perl /users/bnurnberger/Trinity_components/Trinity_defragment_alignments_2seq.pl $component.afa2");

open ($afa2, ">>$component.afa2") or die "Can't append $component.afa2\n";

# ***************************collect blast hits ****************************************************** 

$buffer = '';
my $current_comp;
my $current_number;
my @drop_array_elements;
my $match_no;
my @matched_nucleotides;
my $matchlength;
my $matched_seq;
my @startbuffer = ();
my @endbuffer = ();
my @frame_buffer = ();
my $seq_buffer = -1;   # important to keep track of the sequence number for which data are collected, because last
			 # transfer of data into the definitive arrays after the blast loop can't assume that
			 # the last sequence (at position @indices-1) has any matches 
my $k;

$component =~ /comp(\d+)/;
my $ref_number = $1;

if ($component eq $component_to_follow) {
	print "starting search for blast hits...\n";
}


while ($line = <$blast>) {
	chomp ($line);
	@fields = split ('\t',$line);
	$fields[0] =~ /c(\d+)/;             # extract the component name and compare to $component:
	$current_number = $1;
	if ($current_number > $ref_number) {		# if the current component has a larger index number than than $component 
		last;			                # then quit parsing the blast file. This saves a lot of time.
	}
	$fields[0] =~ /(c\d+_g\d+)/; 
	$current_comp = $1;
	next unless ($current_comp eq $component_to_match);

	for ($i=0;$i<@indices;$i++) {
		$seqname = $component_to_match . '_i' . $indices[$i];   # These are the original Trinity indices of contigs within a component
		if ($line =~ /$seqname\t/) {
			if ($component eq $component_to_follow) { print "match for $seqname\t" }
			
			if ($fields[10] > $fields[11]) {         # determine start and end of the blast match in the original sequence
				$start = $fields[11] - 1;
				$end = $fields[10] - 1;	
			}
			else {
				$start = $fields[10] - 1;
				$end = $fields[11] - 1;
			}


			if ($component eq $component_to_follow) { print "$start-$end\t" }

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

				if ($component eq $component_to_follow) {
					print "$component: match: $start-$end\n";
				}
				@drop_array_elements = ();
				$matches = @startbuffer;
				$flag = 0;			

				for ($j=0; $j<$matches; $j++) {						# check for overlap with existing blast matches for this query
					if (!($end < $startbuffer[$j] || $start > $endbuffer[$j])) {

						if ($startbuffer[$j] >= $start and $endbuffer[$j] <= $end) {    # new match eclipses existing match
							if ($component eq $component_to_follow) {
								print "$component: drop match $j\n";
							}
							push (@drop_array_elements, $j);
						}
						else {
							if ($start < $startbuffer[$j])  {$ol_start = $startbuffer[$j]}   # determine start- and endpoints of the overlap
							else 				{$ol_start = $start}
							if ($end < $endbuffer[$j])  	{$ol_end = $end}
							else 				{$ol_end = $endbuffer[$j]}

				
 if ($component eq $component_to_follow) {
	print "overlap found: $ol_start - $ol_end\n";
 }
							$ol_start -= $start;			# shift $ol_start and $ol_end down to fit them into @matched_nucleotides   
							$ol_end -= $start;
 if ($component eq $component_to_follow) {
       	print "overlap found: rescaled $ol_start - $ol_end, length of blast match: ",scalar @matched_nucleotides,"\n";
 }
							
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
if ($component eq $component_to_follow) {
	print "total overlap: sum = $sum, length of blast hit = $matchlength, theshold: = ", $OVERLAP_THRESHOLD * $matchlength, "\n";
	if ($flag) { print "keep this match\n" }
}
		
					
			}   # match to sequence in buffer
			else {						# first match to the next query sequence
				if ($component eq $component_to_follow) { print "1. match stored\n" }
				$buffer = $fields[0];
				$flag = 1;
				if ($seq_buffer > -1) {
					@{$xt_start{$seq_buffer}} = @startbuffer;	# store the final set of match statistics for a given sequence in three hashes
					@{$xt_end{$seq_buffer}} = @endbuffer;		
					@{$xt_frame{$seq_buffer}} = @frame_buffer;
				}
				$seq_buffer = $indices[$i];    # the sequence number  for which data are currently being collected
				@startbuffer = ();
				@endbuffer = ();
				@frame_buffer = ();
				@drop_array_elements = ();
			}
			
			if ($flag) {   # new match has passed all tests: any eclipsed previous matches are removed 
							# and then the new match is added to output, 
				
				if (@drop_array_elements > 0) {

					if ($component eq $component_to_follow) {
						print "$component: number of elements to drop:", scalar (@drop_array_elements)," size of startbuffer: ",scalar (@startbuffer),"\n";
					}
					for ($j=scalar(@drop_array_elements)-1; $j>=0; $j--) {
						splice (@startbuffer,$drop_array_elements[$j],1); 	# remove start and end pos of the eclipsed match from arrays
						splice (@endbuffer,$drop_array_elements[$j],1);
						splice (@frame_buffer,$drop_array_elements[$j],1);
					}
				}

				push (@startbuffer,$start);
				push (@endbuffer,$end);
				push(@frame_buffer,$frame);

			if ($component eq $component_to_follow) {
				print "stored in buffers: ", scalar @startbuffer," matches\n";
				for ($k=0;$k<@startbuffer;$k++)  {
					print "$k: $startbuffer[$k]-$endbuffer[$k], frame $frame_buffer[$k]\n";
				}
			}

				
			}
		}  # process match to a particular sequence within the current component
	}  # loop over sequence indices 


}  # loop over blast lines
close ($blast);

# After processing the last sequence with read support, the loop over blast lines may continue
# (-> other sequences w/o read support for the same component), and eventually stops when the next component or eof is reached.
# The three arrays (start,end,frame) then still hold the data from the last sequence to be considered. 
# They are copied here into the appropriate hashes, provided that there have been any matches at all ($eq_buffer > -1):
if ($seq_buffer > -1) {
	@{$xt_start{$seq_buffer}} = @startbuffer;
	@{$xt_end{$seq_buffer}} = @endbuffer;		
	@{$xt_frame{$seq_buffer}} = @frame_buffer;
}


# ****************** collect open reading frames *****************************************************


while ($line = <$orf>) {
	chomp ($line);
	if ($line =~ />/) {			
		for ($i=0;$i<@indices;$i++) {
			$seqname = $component_to_match . "_i" . $indices[$i];			
			if ($line =~ /$seqname\:(\d+)-(\d+)/) {   		# matching sequence found: extract orfstart and -end and orientation from header line
				$start = $1;  		# NB: the orf finder counts residues starting at 1, best to keep it that way until frame is determined
				$end = $2;

				if ($line =~ /\+/) {
					$frame = $start % 3;                # if the first nuc of the sequence is the start of the first codon, then the remainder is 1,
										# and so the translation is in the first frame (same coding as the blast output)
					if ($frame == 0) { $frame = 3 }
				}
				else {
					$frame = (($seqlength{$indices[$i]} - $end + 1) % 3);
					if ($frame == 0) 	{ $frame = -3 }
					else 			{ $frame *= -1 }
				}

				push @{$orf_start{$indices[$i]}}, $start - 1;	# now rescale the coordinates to start at 0 and store in arrays
				push @{$orf_end{$indices[$i]}}, $end - 1;		
				push @{$orf_frame{$indices[$i]}}, $frame;	
				
			}		
		}
	}
}

close ($orf);


my $orientation;
my $plus = 0;
my $minus = 0;
my $overlap;
my $keep;

# First look for overlaps between orf pred and blast matches to decide on the correct orientation. Count the number of overlap bases in each direction...

for ($i=0;$i<@indices;$i++) { 	
	$seqno = $indices[$i];
	for ($j=0;$j<@{$xt_start{$seqno}};$j++) {
		for ($k=0;$k<@{$orf_start{$seqno}};$k++) {
			$overlap = check_for_matched_overlap ($seqno,$j,$k);
			if ($overlap > 0) { $plus += $overlap }
			elsif ($overlap < 0) {$minus += abs($overlap) }
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
	
for ($i=0;$i<@indices;$i++) {
	$seqno = $indices[$i];
	for ($j=0;$j<@{$xt_start{$seqno}};$j++) {
		$keep = 0;
		if ($xt_frame{$seqno}[$j] > 0 and $orientation eq '+') { $keep = 1}
		elsif ($xt_frame{$seqno}[$j] < 0 and $orientation eq '-') { $keep = 1}
		if ($keep) {
			print $afa2 ">xt\t$seqno\t$xt_start{$seqno}[$j]\t$xt_end{$seqno}[$j]\t$xt_frame{$seqno}[$j]\n";
		}
	}
	for ($j=0;$j<@{$orf_start{$seqno}};$j++) {
		$keep = 0;
		if ($orf_frame{$seqno}[$j] > 0 and $orientation eq '+') { $keep = 1}
		elsif ($orf_frame{$seqno}[$j] < 0 and $orientation eq '-') { $keep = 1}
		if ($keep) {
			print $afa2 ">orf\t$seqno\t$orf_start{$seqno}[$j]\t$orf_end{$seqno}[$j]\t$orf_frame{$seqno}[$j]\n";
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
my $seq;
my $j;
my $sum = 0;

	for ($i=0;$i<@indices;$i++) {
		$seq = $indices[$i];
		for ($j=0;$j<@{$xt_start{$seq}};$j++) {
			if ($xt_frame{$seq}[$j] > 0) 	{ $sum += $xt_end{$seq}[$j] - $xt_start{$seq}[$j] }
			else 				{ $sum -= $xt_end{$seq}[$j] - $xt_start{$seq}[$j] }
		}
		for ($j=0;$j<@{$orf_start{$seq}};$j++) {
			if ($orf_frame{$seq}[$j] > 0) 	{ $sum += $orf_end{$seq}[$j] - $orf_start{$seq}[$j] }
			else 				{ $sum -= $orf_end{$seq}[$j] - $orf_start{$seq}[$j] }
		}
	}	
	if ($sum > 0) 	{ return '+' }
	else 		{ return '-' }
}
	
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
				 	
	
	
	
	
	
	
	
	

