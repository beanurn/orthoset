#!/usr/bin/perl/
use warnings; use strict;

my $line;
my $buffer;
my $in;
my $fastaout;
my $outlist;
my %last_element;
my $old_name;
my $done;
my @fields;
my $component;
my $counter = 0;
my %besthit;
my $firstmatch = 0;

my $infilename = "var_454_paralogues.txt";
my $outfilename = "Bvv454_paralogues.fasta";
my $outlistname = "Bvv454_paralogues_names.txt";
my $besthitfile = "Bvv454_blastn_trinity_besthit.out";

open ($in,"< $infilename") or die "Can't open $infilename\n";
open ($fastaout,"> $outfilename") or die "Can't create $outfilename\n";
open ($outlist,"> $outlistname") or die "Can't create $outlistname\n";
open (my $best,"< $besthitfile") or die "Can't open $besthitfile\n";

# read in besthit file to select the correct component 

while ($line = <$best>) {
	chomp ($line);
	@fields = split ('\t',$line);
	$besthit{$fields[0]} = $fields[1];
}

close ($best);

my $first = 1;

while ($line = <$in>) {
	chomp($line);
	@fields = split ('\t',$line);
	if ($fields[1] eq 'i0') {		# new paralogue to be added to a trinity component
		if ($first) { $first = 0 }
		elsif ($done == 0) { add_contig_to_component ($old_name,$buffer) }  # previous 454 contig: if the blast besthit component does not have paralogue status
											# then add 454 contig to the first paralogue component in the list
		$fields[2] =~ s/_rev//;		# set up the variables for the next 454 contig
		$old_name = $fields[2];
		$done = 0;
		$firstmatch = 1;
	}
	elsif ($done == 0) {
		$component = $fields[2];	#  a trinity contig
		$component =~ s/_i\d+//;
		$component =~ s/_g/_c/;
		$component =~ s/^c/comp/;
		if ($firstmatch) {		# store the first paralogue component in the list in a buffer and use it in case
						# none of the paralogue components has besthit status
			$buffer = $component;
			$firstmatch = 0;
		}
		if ($besthit{$old_name} =~  /$component/) {   # if the trinity contig is from the same component that scores the best blasthit,
							   # add the 454 contig to that component
			add_contig_to_component ($old_name,$component);
			$done = 1;		# skip over any subsequent components with paralogue status
		}
	}
	$counter++;
#	die if ($counter == 20);
}

if  ($done == 0) { add_contig_to_component ($old_name,$buffer) }

close ($in);
close ($fastaout);
close ($outlist);

sub get_sequence {

	my ($seqname) = @_;

my $flag;
my @headerfields;
my $sequence = '';
my $fastafile = '/users/bnurnberger/Trinity_components/var/trin454/var_nmt_crunch.fasta';
my $seqline;

open (my $fasta,"< $fastafile") or die "Can't open $fastafile\n";

while (my $seqline = <$fasta>) {
      	chomp($seqline);
    	$seqline =~ s/\r//g;
        if ($seqline =~ />/) {
     		if ($flag) {
              		last;
            	}
           	else {                                
			$seqline =~ s/>//;
			$seqline =~ s/\r//g;
			if ($seqline =~ / /) {
				@headerfields = split (' ',$seqline);
				$seqline = $fields[0];
			}
                    	if ($seqline =~ /^$seqname$/) {
                 		$flag = 1;
                        }
                }
    	}
    	elsif ($flag) { $sequence .= $seqline }
}

close ($fasta);
return $sequence;
}

sub add_contig_to_component {

	my ($name454,$comp) = @_;
	

my $sequence;
my $new_name;
my $seqno;

	if (exists $last_element{$comp}) 	{ $seqno = $last_element{$comp} + 1 }
	else 					{ $seqno = 100 }
	$last_element{$comp} = $seqno;
	$new_name = $comp . '_seq' . $seqno;
	print $outlist "$new_name\t$name454\n";    #  write the new name (trinity) and the old (454) name of the 454 contig to a file
	$sequence = get_sequence ($name454);
	print $fastaout ">$new_name\n$sequence\n";   # write the sequence of the 454 contig with its new name to a fasta file

}
