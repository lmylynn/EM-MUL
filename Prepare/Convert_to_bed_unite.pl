#!/usr/bin/perl
#--------------------------------------------------------------------------------------------------------------------------------------------------------------
#This program is to convert ambiguous alignment or unique alignment from SAM format to BED format
#Input: ambiguous alignments or unique alignments (in SAM format)
#Output: BED format
#--------------------------------------------------------------------------------------------------------------------------------------------------------------

#use warnings;
use strict;
use Data::Dumper;
use Storable;
use Getopt::Long qw(GetOptions);
 
my $amb_file; 
my $unique_file;

GetOptions(
    'ambiguous=s' => \$amb_file,
    'unique=s' => \$unique_file,
) or die "Usage: $0 --ambiguous filemane OR $0 --unique filemane\n";
 
if ($amb_file) {
    #print $amb_file;
	unless(open(IN, $amb_file)) {
		print "Could not open ambiguous read file $amb_file!\n";
		exit;
	}
	print "Processing $amb_file...\n";
	open(OUT,">$amb_file.bed"); #output

	while(my $line = <IN>)
	{
		my @array = split(/\t/,$line);
		#print $array[10];
		my $end = $array[3] + length($array[9]) - 1; #Compute end position
		print OUT "$array[2]\t$array[3]\t$end\t$array[0]\t$array[1]\n"; 
	}
}
elsif ($unique_file) {
	#print $unique_file;
	unless(open(IN, $unique_file)) {
			print "Could not open unique read file $unique_file!\n";
			exit;
	}
	print "Processing $unique_file...\n";
	open(OUT,">$unique_file.bed"); #output
	
	while(my $line = <IN>)
	{
		my @array = split(/\t/,$line);
		#print $array[10];
		my $end = $array[3] + length($array[9]) - 1; #Compute end position
		print OUT "$array[2]\t$array[3]\t$end\t$array[0]\t$array[1]\t$array[4]\t$array[5]\t$array[6]\t$array[7]\t$array[8]\t$array[9]\t$array[10]";
	}
}