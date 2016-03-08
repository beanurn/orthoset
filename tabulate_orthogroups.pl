#!/usr/bin/perl
# tabulate_orthogrops.pl
# Counts for every orthogroup the number of distinct sequences for all taxa, the number of components represented by 
# more than one sequence (= paralogues), the number of sequences represented by more than one orf (= frameshifts)
# and writes all these data to a tab file ("groups....tab"). Also, a file of "good" orthogroups is generated which are 
# defined as those that have no Bombina paralogues and at least one sequence of a certain subset of Bombina taxa. This is hard-coded
# at the end of the main programme block. CHANGE AS REQUIRED.
# NB: one sequence of each taxon does not insure overlapping sequence for all the sought after taxa!
# This decision can really only be made after the sequences have been aligned.  

use strict; use warnings;

my @fields;
my $taxon;
my $seqname;
my @bbo_seqs;
my @bor_seqs;
my @bvv_seqs;
my @bvs_seqs;

my %bbo_refs;  # hashes that hold the name of the reference sequence for each component
my %bvv_refs;
my %bvs_refs;
my %bor_refs;

my %seqcount_trin;
my %seqcount_454;

my @taxa = ('bbo','bvv','bor','bvs');
my @zeroes = (0,0,0,0);

@seqcount_trin{@taxa} = @zeroes;  # initialise the hashes that hold counts of distinct sequences for each taxon
@seqcount_454{@taxa} = @zeroes;

my %seqsum;

my $goodGroup;
my $para = 0;
my $multiorf = 0;  # counts the number of contigs that are represented by more than one ORF
my $seqcount_trin = 0;
my $true_para;    # set to 1 as soon as a component with more than 1 sequence including the ref sequence is found, leads to exlclusion of group from orthoset
my $seqcount_454 = 0;  
my $group = 1;
my $para454_count = 0; # counts the number of cases in which a group is excluded due to a 454 paralogue
my $paratrin_count = 0; # counts the number of cases in which a group is excluded due to trinity paralogues
my $line;
my $in;
my $out;
my $xtr = 0;
my $dummy;
my $bom_group_ct = 0;
my $taxcount = 0;
my $compcount = 0;  # counts the number of distinct Trinity components for a given taxon

die "tabulate_orthogroups.pl: filename with  mcl output groups needed (.txt file extention)\n" if (@ARGV != 1); 
my ($orthofile) = @ARGV;

my $tabfile = $orthofile;
$tabfile =~ s/txt/tab/;
my $GoodGroupsFile = $orthofile;
$GoodGroupsFile =~ s/\.txt/_GoodGroups.txt/;


open ($in, "< $orthofile") or die "Can't open $orthofile.\n";
open ($out, "> $tabfile") or die "Can't create $tabfile.\.\n";
open (my $out2, "> $GoodGroupsFile") or die "Can't create $GoodGroupsFile.\n";

foreach $taxon (@taxa) {
	read_in_ref_sequence_names ($taxon);
}

print $out "group\tbbo_t\t454\tpara\tm_orf\tm_comp\tbvv_t\t454\tpara\tm_orf\tm_comp\tbor_t\tpara\tm_orf\tm_comp\tbvs_t\tpara\tm_orf\tm_comp\tbomTaxa\txtr\ttrue_para\n";
while ($line = <$in>) {
	chomp ($line);
	@fields = split(" ",$line);
	$dummy = shift(@fields);  # MCL groups designation
	print $out "RG$group\t";
	while (@fields) {
		$seqname = shift(@fields);
		$taxon = lc(substr($seqname,0,3));
		$seqname =~ s/\w{3}\|//;
		if (!($taxon eq 'xtr')) {
			$seqname =~ /\:/;
			$seqname = $`;
		}
#		print "$taxon,$seqname\t";

		# these arrays hold ORFs. So the count of array members does not give the counts of contigs per taxon
		if    ($taxon eq 'bbo') { push @bbo_seqs,$seqname }
		elsif ($taxon eq 'bor') { push @bor_seqs,$seqname }
		elsif ($taxon eq 'bvv') { push @bvv_seqs,$seqname }
		elsif ($taxon eq 'bvs') { push @bvs_seqs,$seqname }
		else 		        { $xtr++ }
 	}

	# determine, how many groups actually have more than one Bombina taxon in them
	$taxcount = 0;
	if (@bbo_seqs > 0) { $taxcount++ }
	if (@bor_seqs > 0) { $taxcount++ }
	if (@bvv_seqs > 0) { $taxcount++ }
	if (@bvs_seqs > 0) { $taxcount++ }
	if ($taxcount >= 2) { $bom_group_ct ++ }

#	print "\n";
	$goodGroup = 1;
	$true_para = 0;   # this is now the only criterion for exclusion of a group due to paralogues. $true_para is set to 1 whenever the 
			  # following is true: more than one seq per component including the refererence sequence. This needs to occur
			  # only once for any orthogroup in any Bombina taxon to trigger exclusion

	paralogues_multiorfs (\%bbo_refs,\@bbo_seqs,\$seqcount_trin{bbo},\$seqcount_454{bbo},\$para,\$multiorf,\$compcount,\$true_para);
#	print "back in main: para $para, multi $multiorf\n";
	print $out "$seqcount_trin{bbo}\t$seqcount_454{bbo}\t$para\t$multiorf\t$compcount\t";
	$seqsum{bbo} = $seqcount_trin{bbo} + $seqcount_454{bbo};
	if ($para > 0) { $goodGroup = 0 }
	paralogues_multiorfs (\%bvv_refs,\@bvv_seqs,\$seqcount_trin{bvv},\$seqcount_454{bvv},\$para,\$multiorf,\$compcount,\$true_para);
#	print "back in main: para $para, multi $multiorf\n";
	print $out "$seqcount_trin{bvv}\t$seqcount_454{bvv}\t$para\t$multiorf\t$compcount\t";
	$seqsum{bvv} = $seqcount_trin{bvv} + $seqcount_454{bvv};
	if ($para > 0) { $goodGroup = 0 }
	paralogues_multiorfs (\%bor_refs,\@bor_seqs,\$seqcount_trin{bor},\$seqcount_454{bor},\$para,\$multiorf,\$compcount,\$true_para);
#	print "back in main: para $para, multi $multiorf\n";
	print $out "$seqcount_trin{bor}\t$para\t$multiorf\t$compcount\t";
	$seqsum{bor} = $seqcount_trin{bor};
        if ($para > 0) { $goodGroup = 0 }
	paralogues_multiorfs (\%bvs_refs,\@bvs_seqs,\$seqcount_trin{bvs},\$seqcount_454{bvs},\$para,\$multiorf,\$compcount,\$true_para);
#	print "back in main: para $para, multi $multiorf\n";
	print $out "$seqcount_trin{bvs}\t$para\t$multiorf\t$compcount\t";
	$seqsum{bvs} = $seqcount_trin{bvs};
        if ($para > 0) { $goodGroup = 0 }
#	COMMENT OUT THE FOLLOWING LINE TO REMOVE BVS FROM THE PARALOGUE CRITERION
	print $out "$taxcount\t$xtr\t$true_para\n";

# 	choose the desired definition of a good group
	if ($goodGroup == 1) { 
		# at least one sequence of each bbo, bvv, bor, no more than 3 of xtr, no paralogues within components of any Bombina taxon (incl Bvs)
#		if (@bbo_seqs > 0 and @bvv_seqs > 0 and @bor_seqs > 0 and $xtr <= 3) { print $out2 "$line\n" }
		# at least one sequence of each bbo, bvv, bor, no paralogues within components of any Bombina taxon, no restrictions on xtr
#		if (@bbo_seqs > 0 and @bvv_seqs > 0 and @bor_seqs > 0) { print $out2 "$line\n" }
		# at least one sequence of each bbo, bvv, bor, no paralogues within components of these three Bombina taxa, no more than 3 seqs of xtr
#		if (@bbo_seqs > 0 and @bvv_seqs > 0 and @bor_seqs > 0 and $xtr <= 3) { print $out2 "$line\n" }
		# at least one sequence of each bvv and bor, no paralogues within components of any Bombina taxon, no more than 3 xtr seqs
		if ($seqsum{bbo} > 1 and $seqsum{bvv} > 1 and $seqsum{bor} > 1 and $seqsum{bvs} > 1 and $xtr < 4) { print $out2 "$line\n" }
#		if ($seqsum{bor} == 1 and $seqsum{bvv} == 1 and $xtr <= 3) { print $out2 "$line\n" }
	}

	$xtr = 0;
	@bbo_seqs = ();
	@bor_seqs = ();
	@bvv_seqs = ();
	@bvs_seqs = ();
	$group++;
}

print "454 paralogues caused exclusion of a group in $para454_count cases.\n";
print "Pairs of Trinity paralogues caused exclusion of a group in $paratrin_count cases.\n";
print "$bom_group_ct groups have at least 2 Bombina contigs in them \n";

close ($in); 
close ($out);
close ($out2);
	

sub paralogues_multiorfs {

my $seqname;
my $component;
my %orfcount = ();
my %compcount = ();
my %addon_454seq = ();
my $sequences;
my $distinct_components = 0;
my $ref_found;
my $c;
my $s;

	my ($refseq_ref,$array_ref, $trin_ref, $four54_ref,$para_ref, $multi_ref, $comp_ref, $true_para_ref) = @_;

	$$para_ref = 0;
	$$multi_ref = 0;
	$$trin_ref = 0;
	$$four54_ref = 0;
	$$comp_ref = 0;
#	print "para $$para_ref, multi = $$multi_ref\n";
	foreach $seqname (@{$array_ref}) {
		$orfcount{$seqname}++;  		# count the number of occurrences of a given seqname
		if ($seqname =~ /comp/) {
			$component = $seqname;
			$component =~ s/_seq\d+//;
			$compcount{$component}++;	# count the number of occurrences of a given component
			$seqname =~ /_seq(\d+)/;
			if ($1 >= 100) { $addon_454seq{$component} = 1 }
		}
	}

	# identify cases in which there are more than one distinct contig in the same component: paralogues
	# does not apply to 454 ref contigs
	foreach $c (keys %compcount) {
#		print "comp $c\t$compcount{$c}\n";		
		$$comp_ref++;
		$sequences = 0;
		$ref_found = 0;
		foreach $s (keys %orfcount) {	# for each component check whether they are represented by more than one sequence (= paralogues)
			my $long_c = $c . '_';      # underscore distinguishes between comp1_c1_seq1 and comp1_c11_seq1 (seqname) when compared to comp1_c1 (component)
			if ($s =~ /$long_c/) {     
				$sequences++;
				if ($s eq $refseq_ref->{$c}) { $ref_found = 1 }
			 }
		}
		if ($sequences > 1) { 
			($$para_ref)++;
			if ($ref_found) { $$true_para_ref = 1 }   # if the set of seqs (> 1) of a component includes the refseq then this is true paralogy
			if (exists $addon_454seq{$component})   { $para454_count++ }
			else					{ $paratrin_count++ }
		}
	}


	# identify cases of multiple orfs per contig: likely frameshifts
	foreach $s (keys %orfcount) {	
#		print "seq $s\t$orfcount{$s}\n";	
		if ($s =~ /comp/) 	{ $$trin_ref++ }      # count the trinity and 454 contigs separately
		else			{ $$four54_ref++ }	
		if ($orfcount{$s} > 1) { ($$multi_ref)++ }
	}
#	print "before leaving: para $$para_ref, multi = $$multi_ref\n";
}


sub read_in_ref_sequence_names {

my $infilename;
my $infile;
my $seqname;
my $component;

	my ($taxon) = @_;
	
	if ($taxon eq 'bbo') {
		$infilename = '/users/bnurnberger/Trinity_components/bom/full_orthomcl_prep/bom_trinity_refseqs_new.txt';
	}
	elsif ($taxon eq 'bvv') {
		$infilename = '/users/bnurnberger/Trinity_components/var/var_refseqs.txt';
	}
	elsif ($taxon eq 'bor') {
		$infilename = '/users/bnurnberger/Trinity_components/orient/full_orthomcl_prep/orient_refseqs.txt';
	}
	elsif ($taxon eq 'bvs') {
		$infilename = '/users/bnurnberger/Trinity_components/scabra/full_orthomcl_prep/scabra_refseqs.txt';
	}
	else { die "read_in_ref_sequence_names: taxon not recognised\n" }

	open ($infile, "< $infilename") or die "Can't open $infilename.\n";

	while ($seqname = <$infile>) {
		chomp ($seqname);
		$seqname =~ s/^c/comp/;   # convert name to old format
		$seqname =~ s/_g/_c/;	
		$seqname =~ s/_i/_seq/;
		$component = $seqname;
		$component =~ s/_seq\d+//;		
		if ($taxon eq 'bbo') {
			$bbo_refs{$component} = $seqname;
		}			
		elsif ($taxon eq 'bvv') {
			$bvv_refs{$component} = $seqname;
		}			
		elsif ($taxon eq 'bor') {
			$bor_refs{$component} = $seqname;
		}			
		elsif ($taxon eq 'bvs') {
			$bvs_refs{$component} = $seqname;
		}			
	}

	close ($infile);
}
