#!/usr/local/bin/perl
# Mike McQuade
# Mismatch-motifs.pl
# Takes integers k and d, as well as DNA strings,
# and outputs all (k,d)-motifs that appear in each
# string with at most d mismatches.

use strict;
use warnings;

# Initialize variables
my (@nums,$k,$d,@dna,@patterns);
my @alphabet = ('A','C','G','T');

# Open the file to read
open(my $fh,"<ba2a.txt") or die $!;

# Define variables with respective string and integers
@nums = split / /, <$fh>;
$k = $nums[0];
$d = $nums[1];
while (my $temp = <$fh>){
	chomp($temp);
	push @dna,$temp;
}

# Generate all possible k-mers
my @kmers = allKmers();

# Iterate through each possible k-mer
foreach my $kmer (@kmers) {
	# Find the d-Neighbors of the k-mer
	my @dArr = dNeighbors($kmer);

	# Check how many of the DNA strings the substring
	# or any of its d-Neighbors occur in.
	my $matches = 0;
	foreach my $dnaStr (@dna){
		foreach my $motif (@dArr) {
			if (index($dnaStr,$motif) != -1) {
				$matches++;
				last;
			}
		}
		# If the k-mer or its d-Neighbors are found
		# in every DNA string, the k-mer is added
		# to the patterns list ot be returned.
		if ($matches == scalar(@dna)){
			push @patterns,$kmer;
		}
	}
}

# Close the file
close($fh) || die "Couldn't close file properly";

# Print out the patterns
print "@patterns\n";

# Calculates all d-neighbors for the DNA string
sub dNeighbors {
	my (@dneighbors,@kmerArray);
	my $str = $_[0];
	for (my $iteration = 1; $iteration <= $d; $iteration++) {
		if ($iteration == 1) {push @kmerArray,$str}
		else {@kmerArray = @dneighbors;}
		
		# Iterate through each k-mer in the array,
		# changing one letter in the string
		# and adding it to the array
		foreach my $kmer (@kmerArray){
			for (my $i = 0; $i < length($kmer); $i++){
				foreach my $letter (@alphabet){
					my $distinct = $kmer;
					substr($distinct,$i,1) = $letter;
			
					# Only distinct strings are added to
					# the array. If the string is already
					# in the array, it does not get added
					# again.
					if (!grep(/^$distinct$/,@dneighbors)){
						push @dneighbors,$distinct;
					}
				}
			}
		}
	}
	return @dneighbors;
}

# Calculates all possible k-mers
sub allKmers {
	# Initialize variables
	my @words;
	my $currentPower = my $currentChar = 0;

	# Iterate through each place in the k-mer,
	# adding a value in each spot
	for (my $i = $k-1; $i >= 0; $i--) {
		$currentPower = scalar(@alphabet) ** $i;
		for (my $j = 0; $j < scalar(@alphabet) ** $k; $j++) {
			if ($i == $k-1) {push(@words,$alphabet[$currentChar]);}
			else {$words[$j] .= $alphabet[$currentChar];}

			# Decrement current power each time
			$currentPower--;

			# Increment current power each time
			if ($currentPower == 0){
				$currentChar = ($currentChar+1) % scalar(@alphabet);
				$currentPower = scalar(@alphabet) ** $i;
			}
		}
	}
	return @words;
}