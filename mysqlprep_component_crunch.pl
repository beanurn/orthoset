#!/usr/bin/perl/
use strict; use warnings;

my $line;
my $count = 0;
my $i;
my $buffer;
my @fields;
my $flag = 'trinity';

my $infile = shift (@ARGV);
if (@ARGV > 0) {
	$flag = shift (@ARGV);
}

open (my $in,"< $infile") or die "Can't open $infile\n";

my $outfile = $infile;
$outfile =~ s/crunch/crunch_ed/;
open (my $out,"> $outfile") or die "Can't create $outfile\n";

my $component = '';

while ($line = <$in>) {
	
	chomp($line);
	if ($line =~ />seq/) { next }
	elsif ($line =~ />comp/) { 
		$line =~ s/>comp/c/;  
		$line =~ s/_c/_g/;
		$component = $line;
	}
	else { 
		@fields = split ('\t',$line);
		$fields[0] =~ s/seq/i/;
		$line = "$component\t$fields[0]\t$fields[1]\t";
		if ($flag eq '454') {
			if ($fields[2] =~ /comp/) { 
				$fields[2] =~ s/seq/i/;
				$fields[2] =~ s/comp/c/;
				$fields[2] =~ s/_c/_g/;
				$line .= "$fields[2]\tNULL";
			}
			elsif ($fields[2] =~ /rev/) {
				$fields[2] =~ s/_rev//;
				$line .= "$fields[2]\trev";
		 	}
			else { $line .= "$fields[2]\tfor" }
		}
		else { $line .= "$fields[2]" }

		for ($i=3;$i<@fields;$i++) {
			$line = $line . "\t$fields[$i]";
		}
		print $out "$line\n"; 
	}
}

close ($in);
close ($out);




	



