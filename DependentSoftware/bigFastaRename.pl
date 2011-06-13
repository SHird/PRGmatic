#!/usr/bin/perl 
use warnings;
use strict;


#bigFastaRename  -  renames fasta file for use with PRGmatic 
#Copyright (C) 2011 (Sarah M. Hird)
#    This program is free software: you can redistribute it and/or modify
 #   it under the terms of the GNU General Public License as published by
  #  the Free Software Foundation, either version 3 of the License, or
   # (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
 #   but WITHOUT ANY WARRANTY; without even the implied warranty of
  #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   # GNU General Public License for more details.
#    You should have received a copy of the GNU General Public License
 #   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 

my @tags;
my @AoA_tags;


print "Drag the tagfile here: ";
my $tagfile = <STDIN>;
open TAG, "$tagfile" or die "Can't open tagfile: $!";
while (<TAG>) {
	my @tmp = split;
	push @AoA_tags, [@tmp];
}

foreach (0..$#AoA_tags){
	push (@tags, $AoA_tags[$_][1]);
}

print "@tags";

for my $i (0..$#tags){
my $dirname = $tags[$i];
#open (DIR, $dirname) or die "Can't open directory: $!";
#print "opened the $dirname directory\n";

my $fasta = "$dirname/$tags[$i]"."_trimmed.fasta";
my $qual = "$dirname/$tags[$i]"."_trimmed.qual";
my $outputF = "$tags[$i].fasta";
my $outputQ = "$tags[$i].qual";

if (-e "$fasta"){
	open FASTA, "$fasta" or die "Can't open fasta: $!";	
	open OUTPUTF, ">>", "$outputF" or die "Can't open outputF: $!";
	my $count = 1;
	while (<FASTA>) {
		my ($line) = $_; 
		if ($line =~ m/^>(.*)$/){
			print $line, "\t";
			my $countFormat = sprintf("%05s", $count);
			$line= ">$tags[$i]"."_$countFormat";
			print OUTPUTF $line, "\n";
			$count++;
		}#close if
		else {
		print OUTPUTF;
		}#close else
	}#close while
}#close if
else {
	print "file $fasta does not exist\n";
}




if (-e "$qual"){
	open QUAL, "$qual" or die "Can't open qual: $!";
	open OUTPUTQ, ">>", "$outputQ" or die "Can't open outputQ: $!";
	my $count = 1;
	while (<QUAL>) {
		my ($line) = $_; 
		if ($line =~ m/^>(.*)$/){
			print $line, "\t";
			my $countFormat = sprintf("%05s", $count);
			$line= ">$tags[$i]"."_$countFormat";
			print OUTPUTQ $line, "\n";
			$count++;
		}#close if
		else {
			print OUTPUTQ;
		}#close else
	}#close while
}#close if
else {
	print "file $qual does not exist\n";
}#close else

}#close fastaRename