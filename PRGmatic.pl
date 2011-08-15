#!/usr/bin/perl 
use strict;
use warnings;
#S.Hird 08/08/2011
#PRGmatic.v1.6


#PRGmatic: an efficient pipeline for collating NGS data
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
 #   along with this program (in DependentSoftware). 
  #	 If not, see <http://www.gnu.org/licenses/>.
 

my @params;

####get user input####
print "Enter the dataset nickname: ";
my $dataset = <STDIN>;
chomp $dataset;
open SETTINGS, ">", "$dataset.settings.txt" or die $!;

print "############################################\nPARAMETER SETTINGS (default in parentheses):\n############################################\n";
print SETTINGS "############################################\nPARAMETER SETTINGS (default in parentheses):\n############################################\n";

print "Minimum number of reads to call high confidence alleles (5) : ";
$params[0] = <STDIN>;
print SETTINGS "Minimum number of reads to call high confidence alleles (5) : $params[0]";

print "Minimum % identity to call a locus (90): ";
$params[1] = <STDIN>;
print SETTINGS "Minimum % identity to call a locus (90): $params[1]";

print "Minimum overlapping sequence to cluster reads into alleles or loci (100 for 454 data, 50 for Illumina): ";
$params[5] = <STDIN>;
print SETTINGS "Minimum overlapping sequence to cluster reads into alleles or loci (100 for 454 data, 50 for Illumina): $params[5]";

print "Minimum coverage for calling consensus sequence in individual (6): ";
$params[2] = <STDIN>;
print SETTINGS "Minimum coverage for calling consensus sequence in individual (6): $params[2]";

print "Minimum coverage for calling a SNP (3): ";
$params[3] = <STDIN>;
print SETTINGS "Minimum coverage for calling a SNP (3): $params[3]";

print "Minimum % of reads for calling SNP (20): ";
$params[4] = <STDIN>;
print SETTINGS "Minimum % of reads for calling SNP (20): $params[4]";

chomp @params;

#foreach (@params){
#	print "$_ \n";
#	}


####read in fasta files

my $dirname = "inputFASTA";
opendir(DIR, $dirname) or die "can't opendir $dirname: $!";
my @allFilesIF = readdir(DIR);
close (DIR);
my $stringOfNames;
my @names;
my @namesFasta;

foreach my $f (@allFilesIF){
	if ( ($f ne ".") && ($f ne "..") && ($f !~ /^\.DS_Store/i) && ($f =~ /.fasta$/)){
   		$stringOfNames .= "$dirname/$f " ;
   		my $name = $f;
   		$name =~ s/.fasta//;
   		print "$name ";
   		push (@names, $name);
   		push (@namesFasta, $f);
   	}
}

my $numIndividuals = @names;
print "numInds is $numIndividuals \n";
my $totalFasta = "total_$dataset.fasta";
print "$totalFasta\n";


###execute subroutines####

cap3Execute();
chromoFormat();
autoBWAandSamtools();
pileupTables();
writeAlleles();
calcHeterozygosity();
countIndividualsInLoci();
muscleAlignments();
#multiCompute();



##########################################################################################################
################  BEGIN SUBROUTINES  #####################################################################
##########################################################################################################

sub cap3Execute {

my @goodContigNames;

foreach (@namesFasta){
		my $current = $_;
		system"./cap3 inputFASTA/$current -p 99 -o $params[5] -z 5 -c 20";
}#close foreach

my @contigArray;
my @qualArray;

for my $i (0..$#names){
#	@contigArray = ();
	my $contigFile = "$names[$i].fasta.cap.contigs";
	my $qualContig = "$names[$i].fasta.cap.contigs.qual";
	print "contigFile is $contigFile\n";
	open CONTIG, "inputFASTA/$contigFile" or die $!;
	open QUAL, "inputFASTA/$qualContig" or die $!;
	my @contigTemp = <CONTIG>;
	my @qualTemp = <QUAL>;
	chomp @contigTemp;
	chomp @qualTemp;
#	pop @contigTemp;
	my $sequence;
	my $currentLine;

	for $currentLine (@contigTemp){
		if ($currentLine =~ /^[ACGTN]/){
			$sequence .= $currentLine ;
		}#close if
		else {
			$currentLine .= "_$names[$i]";
			# $sequence .= "";
			push (@contigArray, $sequence);
			push (@contigArray, $currentLine);
			$sequence = "";
		}	#close else
	}#close $currentLine
	push (@contigArray, $sequence);

	my $quality;
	my $currentQual;

	for $currentQual (@qualTemp){
		if ($currentQual =~ /^[0-9]/){
			$quality .= $currentQual;
		} #close if
		else {
			$currentQual .= "_$names[$i]";
			push (@qualArray, $quality);
			push (@qualArray, $currentQual);
			$quality = "";
		} #close else
	}#close currentQual	
	push (@qualArray, $quality);
}#close $i
#print "qualArray: @qualArray\n";

####this picks out the sequences with greater than 5 coverage and puts them into a good contig names array####
for my $j (0..$#names){
	my $aceFile = ("$names[$j]".".fasta.cap.ace");
	open ACE, "<", "inputFASTA/$aceFile" or die "Can't open ace file: $!";
	while (<ACE>){
    	chomp;
    	if (/^CO/){ # finds the line with contig info
			m/CO\s+\S+\s+(\d+)\s+(\d+)\s+\d+/;
			if ($1 < 700 && $2 >= $params[0]){
			   # print "$_";
			    my @split = split(/ /, $_);
			    push (@goodContigNames, ("$split[1]"."_$names[$j]")) ;
			}#close if $1 $2
		}#close if ^CO
	}#close while ACE
	close ACE;
}#close foreach @names
#print "good contig names: @goodContigNames\n";
#print "contig array1: $contigArray[0]     $contigArray[1] \n";
shift @contigArray;
shift @qualArray;
#print "contig array1: $contigArray[0]     $contigArray[1] \n";
foreach (@contigArray){
	print $_, "\n";
}
foreach (@qualArray){
	print $_, "\n";
}

#print @contigArray;
my %contigHash;
my %qualHash;

####put all contig seqs and qualities into a hashes
for (my $i=0; $i<$#contigArray; $i++){
	if ($contigArray[$i] =~ /^>/){
	$contigHash{$contigArray[$i]} = "$contigArray[$i+1]";
	} #close if
}#close for

for (my $ii=0; $ii<$#qualArray; $ii++){
	if ($qualArray[$ii] =~ /^>/){
	$qualHash{$qualArray[$ii]} = "$qualArray[$ii+1]";
	}
}

while( my ($k, $v) = each %qualHash ) {
        print "key: $k, value: $v.\n";
    }

open OUTPUT, ">>", "$totalFasta" or die $!;
for (my $i = 0; $i<=$#goodContigNames; $i++){
	print OUTPUT ">$goodContigNames[$i]\n";
	print OUTPUT $contigHash{">$goodContigNames[$i]"}, "\n";
}

open Q_OUT, ">>", "$totalFasta.qual" or die $!;
for (my $i = 0; $i<=$#goodContigNames; $i++){
	print "$goodContigNames[$i]\n";
	print Q_OUT ">$goodContigNames[$i]\n";
	print Q_OUT $qualHash{">$goodContigNames[$i]"}, "\n";
}
print "cap3 $totalFasta p $params[1] -o $params[5] -z 5 -c 20\n";

system"./cap3 $totalFasta -p $params[1] -o $params[5] -z 5 -c 20";

} #close cap3Execute


###########################################

sub chromoFormat {
#converts the .cap.contigs file from cap3 into the pseudo-reference genome

system"cat $totalFasta.cap.contigs $totalFasta.cap.singlets >$dataset.alleles.fasta";
my $fastaFile = ("$dataset.alleles.fasta");


open FASTA, "$fastaFile" or die "Can't open FASTA file";
my @fastaTemp = <FASTA>;
chomp @fastaTemp;
my $sequence;
my $current;
my @fastaArray;

for my $current (@fastaTemp){
	
	if ($current =~ /^[ACGT]/){
		$sequence .= $current ;
	}
	else {
		$sequence .= "\n";
		push (@fastaArray, $sequence);
		push (@fastaArray, $current);
		$sequence = "";	
	}
}
push (@fastaArray, $sequence);
$fastaArray[0] = "0";

print "0 is $fastaArray[0], 1 is $fastaArray[1], 2 is $fastaArray[2]\n";

for (my $j=1; $j<(@fastaArray-1); $j++){

my $length = ((length $fastaArray[$j+1])-1);

if ($j % 2 != 0) {

my $locusNumber = (($j + 1)/2);
$fastaArray[$j] = ">gi|$length|$dataset|$dataset"."_$locusNumber\n";
}
}
shift @fastaArray;
print @fastaArray;
my $outputFile = "$dataset".".fna";
open OUTPUT, ">>", "$outputFile" or die "Can't open file for output";
print OUTPUT @fastaArray;

} #close chromoFormat;

##########################################################################################################


sub autoBWAandSamtools{
#Automates blasting of fragments to pseudo ref genome using BWA and Samtools


#set the p-rg name and index in bwa
my $prgSeqName = ("$dataset".".fna");
system "./bwa index -p $dataset $prgSeqName";

#blast each input fasta file to p-rg
#my $dirname = "inputFASTA";
#opendir(DIR, $dirname) or die "can't opendir $dirname: $!";
#my @allFiles = readdir(DIR);
#close (DIR);

foreach my $f (@names){   		
  	my $fastaFileToBlast = ("$dirname/$f.fasta"); 
  	my $outputPrefix = $f;
	print $fastaFileToBlast;
	print $outputPrefix;

	system "./bwa aln $dataset $fastaFileToBlast >$outputPrefix.sai";
	system "./bwa samse $dataset $outputPrefix.sai $fastaFileToBlast >$outputPrefix.sam";
	system "./samtools view -bS -o $outputPrefix.bam $outputPrefix.sam";
	system "./samtools sort $outputPrefix.bam $outputPrefix.sorted";
	system "./samtools faidx $prgSeqName";
	system "./samtools index $outputPrefix.sorted.bam";
	system "./samtools view -u $outputPrefix.sorted.bam | ./samtools pileup -vcf $prgSeqName -> raw_$outputPrefix.txt";
	system "./samtools pileup -f $prgSeqName $outputPrefix.sorted.bam >$outputPrefix.pileup";
#	system "java -jar VarScan.v2.2.2.jar pileup2snp $outputPrefix.pileup --min-coverage $params[2] >$outputPrefix.snp";
#	system "java -jar VarScan.v2.2.2.jar pileup2cns $outputPrefix.pileup --min-coverage $params[3] >$outputPrefix.cns";
	system "java -jar VarScan.v2.2.2.jar pileup2indel $outputPrefix.pileup --min-coverage $params[3] >$outputPrefix.indel";
		
	}#close $f
}#close autoBWAandSamtools

##########################################################################################################


sub pileupTables {
#rewrites pileup files as summary table of counts for each base

foreach my $nameCurrent (@names){
#my $nameCurrent = $_;
my $pileupFile .= ("$nameCurrent".".pileup");
my @pileupAoA = ();
my $linesPILEUP = 0;

open PILEUP, "$pileupFile" or die "Can't open pileup $!";
while (<PILEUP>) {
	my @tmp = split;
	push @pileupAoA, [ @tmp ];
} #close while
close PILEUP;	

$linesPILEUP = ($#pileupAoA+1);
print $linesPILEUP, "\t";
print $#pileupAoA, "\t"; 
print $pileupAoA[0][0];

my @currentLocus = split(/\|/, $pileupAoA[0][0]);
my $locus = $currentLocus[3];
my $length = $currentLocus[1];

print "locus is $locus ; length is $length\n";

my @counts = qw(loc pos base . a c g t - ins del coverage);
open OUTPUT, ">>$pileupFile.Counts.txt" or die $!;
print OUTPUT "@counts\n";
@counts = ();

for my $i (0..$#pileupAoA){
	print "this pileup # $i \n";
	if ($pileupAoA[$i][1] > $pileupAoA[($i-1)][1]){
		$counts[0] = "$locus";
		$counts[1] = "$pileupAoA[$i][1]";
		$counts[2] = "$pileupAoA[$i][2]";
		$counts[3] = "0";
		$counts[4] = "0";
		$counts[5] = "0";
		$counts[6] = "0";
		$counts[7] = "0";
		$counts[8] = "0";
		$counts[9] = "0";
		$counts[10] = "0";
		$counts[11] = "$pileupAoA[$i][3]";
		
		my @basesFromPileup = split(//, $pileupAoA[$i][4]);
		
		for my $b (0..$#basesFromPileup) {
			#my $baseCurrent = $_;
			#print "now on $b\n";
			if($basesFromPileup[$b] =~ m/[0-9]/){
				print "found a number!  $basesFromPileup[$b] at $b\n";
				my $delete1 = ($basesFromPileup[$b] + $b);
				print $delete1, "\n";
				for (my $delete = $delete1; $delete >= $b; $delete--){
					print "$delete ";
					$basesFromPileup[$delete] = "";
					print "@basesFromPileup\n";
				}
			}		
			elsif ($basesFromPileup[$b] =~ m/\./){
				$counts[3]++;}
			elsif($basesFromPileup[$b] =~ m/A/){
				$counts[4]++;}
			elsif($basesFromPileup[$b] =~ m/C/){
				$counts[5]++;}
			elsif($basesFromPileup[$b] =~ m/G/){
				$counts[6]++;}
			elsif($basesFromPileup[$b] =~ m/T/){
				$counts[7]++;}
			elsif($basesFromPileup[$b] =~ m/\*/){
				$counts[8]++;}
			elsif($basesFromPileup[$b] =~ m/\+/){
				$counts[9]++;}
			elsif($basesFromPileup[$b] =~ m/\-/){
				$counts[10]++;}
			
				#delete next N bases}
		}#close @basesFromPileup

	print OUTPUT "@counts\n";
	#print OUTPUT "c0 $counts[0] c1 $counts[1] c2 $counts[2] c3 $counts[3] c4 $counts[4] c5 $counts[5] c6 $counts[6] c7 $counts[7] c8 $counts[8] c9 $counts[9]\n";
	}#close if
	
	else {
		@currentLocus = split(/\|/, $pileupAoA[$i][0]);
		$locus = $currentLocus[3];
		$length = $currentLocus[1];
		$counts[0] = "$locus";
		$counts[1] = "$pileupAoA[$i][1]";
		$counts[2] = "$pileupAoA[$i][2]";
		$counts[3] = "0";
		$counts[4] = "0";
		$counts[5] = "0";
		$counts[6] = "0";
		$counts[7] = "0";
		$counts[8] = "0";
		$counts[9] = "0";
		$counts[10] = "0";
		$counts[11] = "$pileupAoA[$i][3]";

		my @basesFromPileup = split(//, $pileupAoA[$i][4]);
		
		for my $b (0..$#basesFromPileup) {
			#my $baseCurrent = $_;
			if($basesFromPileup[$b] =~ m/[0-9]/){
				print "found a number!  $basesFromPileup[$b] at $b\n";
				my $delete1 = ($basesFromPileup[$b] + $b);
				print $delete1, "\n";
				for (my $delete = $delete1; $delete >= $b; $delete--){
					print "$delete ";
					$basesFromPileup[$delete] = "";
					print "@basesFromPileup\n";
				}
			}			
			elsif ($basesFromPileup[$b] =~ m/\./){
				$counts[3]++;}
			elsif($basesFromPileup[$b] =~ m/A/){
				$counts[4]++;}
			elsif($basesFromPileup[$b] =~ m/C/){
				$counts[5]++;}
			elsif($basesFromPileup[$b] =~ m/G/){
				$counts[6]++;}
			elsif($basesFromPileup[$b] =~ m/T/){
				$counts[7]++;}
			elsif($basesFromPileup[$b] =~ m/\*/){
				$counts[8]++;}
			elsif($basesFromPileup[$b] =~ m/\+/){
				$counts[9]++;}
			elsif($basesFromPileup[$b] =~ m/\-/){
				$counts[10]++;}
			
		}#close @basesFromPileup
	print OUTPUT "@counts\n";
	#print OUTPUT "c0 $counts[0] c1 $counts[1] c2 $counts[2] c3 $counts[3] c4 $counts[4] c5 $counts[5] c6 $counts[6] c7 $counts[7] c8 $counts[8] c9 $counts[9]\n";
	}#close else
}#close for $i
print @names, "\n";

} #close foreach @names
}#close pileupTables

#############################################################

sub writeAlleles {

print @names;
foreach my $nameForWrite (@names){
	#my $nameForWrite = $_;
	my @firstAllele = ();
	my @secondAllele = (); 
	my @PTableAoA = ();
	my $pileupTable = ("$nameForWrite".".pileup.Counts.txt");

	open PTABLE, "$pileupTable" or die "Can't open Pileup Table $!";
	while (<PTABLE>) {
		my @tmp = split;
		push @PTableAoA, [ @tmp ];
	}#close while
	close PTABLE;

my $chromoNumber;

print $PTableAoA[0][0];
for my $i (1..$#PTableAoA) {
$chromoNumber = $PTableAoA[$i][0];
$chromoNumber =~ s/$dataset//;
$chromoNumber =~ s/_//;
print "this is chromoNumber $chromoNumber\n\n\n\n\n\n\n";
#	$outputPrefix =~ s/.fasta//;
print "num rows $#PTableAoA\n";
	if ($PTableAoA[$i][11]>=$params[2]){
	
	my $multiHits = 0;
	for my $counter (3..8) {
		if ($PTableAoA[$i][$counter] > 0) {
			$multiHits++;
		}
	}
	if ($multiHits >2) {
		open MHL, ">>", "MultiHitLoci.txt" or die "Can't open MultiHitLoci.txt";
		print MHL "Individual $nameForWrite at locus $PTableAoA[$i][0] has >2 bp at position $PTableAoA[$i][1]\n";
		close MHL;
	}
	
	
	
	#	if ( ($PTableAoA[$i][3] > $params[3]) && (($PTableAoA[$i][3]/$PTableAoA[$i][11]) >= $params[4]) ) {  #consensus base in more than cutoff reads
			
			if ($PTableAoA[$i][3] == $PTableAoA[$i][11]){ #if there is NO discrepencies, just call both alleles
				$firstAllele[$chromoNumber][($PTableAoA[$i][1])] = "$PTableAoA[$i][2]";
				print "$PTableAoA[$i][0] $PTableAoA[$i][1] is $PTableAoA[$i][2] \n";
				$secondAllele[$chromoNumber][($PTableAoA[$i][1])] = "$PTableAoA[$i][2]";
			}
				
			else {
				print "ELSE! 0 $PTableAoA[$i][0] , 1 $PTableAoA[$i][1] , 3 $PTableAoA[$i][3]\n";
				if ( ($PTableAoA[$i][3] >= $params[3]) && (($PTableAoA[$i][3]/$PTableAoA[$i][11]) >= ($params[4]/100)) ) {				
					$firstAllele[$chromoNumber][($PTableAoA[$i][1])] = "$PTableAoA[$i][2]";  #subject to change
					print 	"firstAllele $PTableAoA[$i][0] ($PTableAoA[$i][1]) = $PTableAoA[$i][2] \n";  #subject to change
					$secondAllele[$chromoNumber][($PTableAoA[$i][1])] = "$PTableAoA[$i][2]";  #subject to change			
					for my $j (4..8) {
						print "$PTableAoA[$i][$j] \n";
						if ( ($PTableAoA[$i][$j] >= $params[3]) && (($PTableAoA[$i][$j]/$PTableAoA[$i][11]) >= ($params[4]/100) ) ) {
							
							if ( ($PTableAoA[$i][$j] <= $PTableAoA[$i][3]) ){
								$secondAllele[$chromoNumber][($PTableAoA[$i][1])] = "$PTableAoA[0][$j]";
								print "ELSE1 secondAllele $PTableAoA[$i][0] $PTableAoA[$i][1] = $PTableAoA[0][$j]\n";

							} #if
							else { 
								$firstAllele[$chromoNumber][($PTableAoA[$i][1])] = "$PTableAoA[0][$j]";
								print "ELSE2 firstAllele $PTableAoA[$i][0] ($PTableAoA[$i][1]) = $PTableAoA[0][$j]\n";


							}#close else
						}#close if	
					}#close $j
				}#close if
		#}		
	
				else { #when there are fewer . than the snp cutoff
					for my $j (3..8) {
						if ( ($PTableAoA[$i][$j] >= $params[3]) && (($PTableAoA[$i][$j]/$PTableAoA[$i][11]) >= ($params[4]/100) ) ) {
							if ( ($PTableAoA[$i][$j]/$PTableAoA[$i][11]) >= 0.5) {
								$firstAllele[$chromoNumber][($PTableAoA[$i][1])] = "$PTableAoA[0][$j]";
								$secondAllele[$chromoNumber][($PTableAoA[$i][1])] = "$PTableAoA[0][$j]";
								print "they were fewer . called first and second $PTableAoA[0][$j]\n";
								next;
							}
							else {
								if ($j == 3){
									$secondAllele[$chromoNumber][($PTableAoA[$i][1])] = "$PTableAoA[$i][2]";  #subject to change
									print "Called second $PTableAoA[$i][2]\n";
								}
								else{
									$secondAllele[$chromoNumber][($PTableAoA[$i][1])] = "$PTableAoA[0][$j]";
									print "Called second $PTableAoA[0][$j]\n";
								}							
							}#close else
						}#close if
						elsif ( ($PTableAoA[$i][3] < $params[3]) && ($PTableAoA[$i][4] < $params[3]) && ($PTableAoA[$i][5] < $params[3]) && ($PTableAoA[$i][6]< $params[3]) && ($PTableAoA[$i][7] < $params[3]) && ($PTableAoA[$i][8] < $params[3]) ) {
							$firstAllele[$chromoNumber][($PTableAoA[$i][1])] = "N";
							$secondAllele[$chromoNumber][($PTableAoA[$i][1])] = "N";
							print "not enough coverage anywhere. both strands called N \n";
						}#close elsif
					}#close for	
				}#close else
		}#close else
	}#close if coverage
}#close for $i			
	
	
####INSERT INSERTIONS INTO ALLELES#######
my @indelAoA = ();

my $indelTable = ("$nameForWrite".".indel");
open INDEL, "$indelTable" or die "Can't open indel table $!";
while (<INDEL>) {
	my @tmp = split;
	push @indelAoA, [ @tmp ];
}#close while
close INDEL;

for my $i (1..$#indelAoA) {
	my $basesToInsert;
	if ( ($indelAoA[$i][3] =~ m/^\+/) && ($indelAoA[$i][5] >= $params[3]) && (($indelAoA[$i][5])/($indelAoA[$i][4]+$indelAoA[$i][5]) >= ($params[4]/100) ) ) {	
		my @splitName = split(/\|/, $indelAoA[$i][0]);
		
		#make a variable of the bases to insert by removing the + from the 4th column of indel table
		my @splitInsert = split(//, $indelAoA[$i][3]);
		shift(@splitInsert);
		$basesToInsert = join('', @splitInsert);
		
		if ( $indelAoA[$i][5]/($indelAoA[$i][4]+$indelAoA[$i][5]) >= 0.50 ) { #if the insertion is more than half the reads, put on strand 1
			$firstAllele[$splitName[3]][$indelAoA[$i][1]] .= "$basesToInsert";
			print "greater than 50, first allele now $firstAllele[$splitName[3]][$indelAoA[$i][1]] at firstAllele[$splitName[3]][$indelAoA[$i][1]]\n";
			if ( ($indelAoA[$i][4] < $params[3]) || ($indelAoA[$i][4]/($indelAoA[$i][4]+$indelAoA[$i][5]) <= ($params[4]/100)) ) {  #if var1 doesn't meet snp cutoffs, put insertion on both strands
				$secondAllele[$splitName[3]][$indelAoA[$i][1]] .= "$basesToInsert";
				print "both alleles, now $firstAllele[$splitName[3]][$indelAoA[$i][1]] $secondAllele[$splitName[3]][$indelAoA[$i][1]] at $firstAllele[$splitName[3]][$indelAoA[$i][1]] secondAllele[$splitName[3]][$indelAoA[$i][1]]\n";
			}#close if
		}#close if
		elsif ( $indelAoA[$i][5]/($indelAoA[$i][4]+$indelAoA[$i][5]) < 0.50 ) { #if the strand is less than half the reads, put on strand 2
			$secondAllele[$splitName[3]][$indelAoA[$i][1]] .= "$basesToInsert";
			print "less than 50, second allele now $secondAllele[$splitName[3]][$indelAoA[$i][1]] at secondAllele[$splitName[3]][$indelAoA[$i][1]]\n";
		}#close elsif
		else{
		print "nothing significant at $i\n";}
	} #close if	
	else{
		print "still nothing significant at $i\n";}
}#close for $i	
	

####PRINT FIRST ALLELE TO OUTPUT######
		for my $row (1..$#firstAllele){
			if ($firstAllele[$row][1]) {
				open LOCI_OUT, ">>", ("calledAlleles/"."$dataset"."_$row.fasta") or die "Can't write output.";
				print LOCI_OUT ">$nameForWrite.01\n";
				for my $column (1..$#{$firstAllele[$row]}){
					if ($firstAllele[$row][$column] ne "-"){  #delete deletions
						print LOCI_OUT "$firstAllele[$row][$column]";
						print "firstAllele[$row][$column]  $firstAllele[$row][$column]\n";
					}
				} #close $column
				print LOCI_OUT "\n";
				close LOCI_OUT;
			} #close if
		} #close $row
	
####PRINT SECOND ALLELE TO OUTPUT######
		for my $row (1..$#secondAllele){
			if ($secondAllele[$row][1]) {
				open LOCI_OUT, ">>", ("calledAlleles/"."$dataset"."_$row.fasta") or die "Can't write output.";
				print LOCI_OUT ">$nameForWrite.02\n";
				for my $column (1..$#{$secondAllele[$row]}){ #delete deletions
					if ($secondAllele[$row][$column] ne "-"){
						print LOCI_OUT "$secondAllele[$row][$column]";
					}
				} #close $column
				print LOCI_OUT "\n";
				close LOCI_OUT;
			}#close if
		}#close $row
	}#close @names


}#close writeAlleles



##########################################################################################################


sub calcHeterozygosity{
#calculates Ho and He for each locus in the calledAlleles folder


open PRINT, ">>", "$dataset.HoHe.txt" or die "Can't write to HoHe.txt";
print PRINT "locusFile\tHeterozygotes\tnumber_individuals\tnumber_alleles\tHo\tHe\n";
close PRINT;

my $dirnameCA = "calledAlleles";
opendir(DIR, $dirnameCA) or die "can't opendir $dirnameCA: $!";
my @allFilesCA = readdir(DIR);
close (DIR);

foreach my $f (@allFilesCA){
		if ( ($f ne ".") && ($f ne "..") && ($f !~ /^\.DS_Store/i) && ($f !~ /^aligned/) ){
	my @AoA= (); #AoA of the alleles
	my @AoA_alleles = ();
	my $countHeterozy = 0;
	my $linesLOCUS = 0;
	my $sumFreqSq = 0;
#	my $locusFile = "$dataset$i.fasta";
	open LOCUS, "$dirnameCA/$f" ;#or die "Can't open .fasta file: $!";
		while (<LOCUS>) {
			my @tmp = split;
			push @AoA, [ @tmp ];
			$linesLOCUS++;
		} #close while
	close LOCUS;	
	print $linesLOCUS, "linesLOCUS \n";
	my $inds = $linesLOCUS/4;
	my $total_alleles = $linesLOCUS/2;
	
	$AoA_alleles[0][0] = $AoA[1][0];
	$AoA_alleles[0][1] = 1;
	my $heterozygotes = 0;
	my $numberAlleles = 1;
	for (my $j = 3; $j<=($linesLOCUS-1); $j+=2){
		
		
		
		for my $k (0..($numberAlleles-1)){
		
			if ($AoA[$j][0] eq $AoA_alleles[$k][0]){	
				print "k is $k j is $j \n";
				$AoA_alleles[$k][1]++;
				print "AoA[ $j ][0] eq AoA_alleles[ $k ][0] so AoA_alleles[ $k ][1] is $AoA_alleles[$k][1] \n";
				last;
			} #close if
			elsif ($k==($numberAlleles-1)){
				$AoA_alleles[$numberAlleles][0] = $AoA[$j][0];
				print "k is $k j is $j \n";
#				print "AoA_alleles[numbAllelPlus1][0] = $AoA_alleles[$numberAlleles+1][0]\n";
				$AoA_alleles[($numberAlleles)][1] = 1;
				$numberAlleles++;
				print "NumbAlleles is $numberAlleles\n";
				last;
			}#close elsif
			else { 
				$k++;
				next;
			} #close else
		}#close my $k
	
	
	}#close for $j
	
#####SUMMARIZE DATA/COUNTS AND CALC He#####

#count heterozygotes
	for (my $j = 1; $j<=($linesLOCUS-1); $j+=4){
		if ($AoA[$j][0] ne $AoA[$j+2][0]){
			$countHeterozy++;
		}#close if
	}#close $j
	print "heterozygotes = $countHeterozy \n";
	
#determine allele frequencies
	for my $row (0..($numberAlleles-1)){
		$AoA_alleles[$row][2] = ($AoA_alleles[$row][1]/($total_alleles));
		$sumFreqSq += ($AoA_alleles[$row][2]*$AoA_alleles[$row][2]);
		for my $column (0..2){
			print $AoA_alleles[$row][$column], "\t";
			}
		print "\n";
	}
my $obsHet = $countHeterozy/$inds;
my $expHet = 1 - $sumFreqSq;
print "For $f, observed heterozygosity is $obsHet and expected heterozygosity is $expHet.\n";
open PRINT, ">>", "$dataset.HoHe.txt" or die "Can't write to HoHe.txt";
print PRINT "$f\t$countHeterozy\t$inds\t$numberAlleles\t$obsHet\t$expHet\n";
close PRINT;
	
	}#close if exists
	
} #close $i

}#close calcHeterozygosity

##########################################################################################################


sub countIndividualsInLoci {
# script for summarizing how many individuals are in each locus after being called by the pipeline

#my @names = ();
my @lociNumbers;
my $loci;
my $totalLoci = 0;


####calculate how many loci in prg (IS THAT IT? GLOBAL VARIABLE???)
my $dirnameCA = "calledAlleles";
opendir(DIR, $dirnameCA) or die "can't opendir $dirnameCA: $!";
my @allFilesCA = readdir(DIR);
close (DIR);
foreach my $f (@allFilesCA){
	if ( ($f ne ".") && ($f ne "..") && ($f !~ /^\.DS_Store/i) && ($f !~ /^clustalw/) ){
		my $outputPrefix = $f;
		$outputPrefix =~ s/.fasta//;
		my @splitF = split(/_/, $outputPrefix);
		print $splitF[0], "\t", $splitF[1], "\n";
		push (@lociNumbers, $splitF[1]);
		
		if ($splitF[1]>$totalLoci){
			print "totalLoci update; totalLoci is $totalLoci , splitF[1] is $splitF[1], ";
			$totalLoci = $splitF[1];
			print "totalLoci is NOW $totalLoci. \n";
		}#close if
	} #close if
} #close f

push (@names, "total");
unshift (@names, "name");
unshift (@lociNumbers, "name");
my $numColumns = @names;
#chomp (@names);
my $prefix = ("$dataset"."_");  #the names of the loci files minus the number and .fasta
print "numColumns $numColumns , prefix $prefix , names0 $names[0] , names1 $names[1] , numRows $totalLoci, lociNumbers0 $lociNumbers[0]\n";

my @AoA = ();
my @AoAcount = ();

$AoAcount[$totalLoci+1][0] = "total";

foreach my $i (1..($numColumns-1)) {
	$AoAcount[0][$i] = $names[$i];
	foreach my $j (1..$totalLoci) {
		$AoAcount[$j][0] = "Loc$j";
		$AoAcount[$j][$i] = "0";
	}
}
	my @loci = ();
	my $currentLine = 0;
	
	foreach my $f (@allFilesCA){
		if ( ($f ne ".") && ($f ne "..") && ($f !~ /^\.DS_Store/i) && ($f !~ /^clustalw/) ){

		my $outputPrefix = $f;
			$outputPrefix =~ s/.fasta//;
			my @splitF = split(/_/, $outputPrefix);
			$currentLine = $splitF[1];
			
	open LOCI, "$dirnameCA/$f" or print "Can't open locus file: $f.\n";
		@loci = <LOCI>;
		chomp (@loci);
	my $numLinesInLOCI = @loci;
	
	print "@loci\n";
	print "Opened $f. \n";
	print $numLinesInLOCI, "\n";
	
	for (my $j=0; $j<$numLinesInLOCI; $j+=4){
		print "Unrevised locus is $loci[$j] \t";
		$loci[$j] =~ s/>//;
		$loci[$j] =~ s/\.01//;
		print "revised locus is $loci[$j] \n";
		for my $m (1..($numColumns-2)){
			if ($loci[$j] eq $names[$m]){
			print "$loci[$j] equals $names[$m]\n";
				$AoAcount[$currentLine][$m] = "1";
				next;
			}
		}
	}
close LOCI;
}	
}

####sum each row
for my $row (1..$totalLoci){
	for my $column (1..($numColumns-2)){
	$AoAcount[$row][($numColumns-1)] += $AoAcount[$row][$column];
	}
}
####sum each column
for my $column (1..($numColumns-1)){
	for my $row (1..$totalLoci){
	$AoAcount[$totalLoci+1][$column] += $AoAcount[$row][$column];
	}
}
$AoAcount[0][0] = "ind/loc";
####print output
for my $row (0..($totalLoci+1)){
	print "row is $row \t";
	open LOCI_OUT, ">>", "$dataset.counts.txt";
	for my $column (0..($numColumns-1)){
		print LOCI_OUT "$AoAcount[$row][$column]\t";
	}
	print LOCI_OUT "\n";
}
	
} #close countIndividualsInLoci

##########################################################################################################


sub muscleAlignments {

my $dirnameCA = "calledAlleles";
opendir(DIR, $dirnameCA) or die "can't opendir $dirnameCA: $!";
my @allFilesCA = readdir(DIR);
close (DIR);

foreach my $f (@allFilesCA){
	if ( ($f ne ".") && ($f ne "..") && ($f !~ /^\.DS_Store/i) && ($f !~ /^clustalw/) ){
	print "$f\t";
	my $nameOut = $f;
	$nameOut =~ s/.fasta/.aln.fasta/;
	print "$nameOut \n";
	system"./muscle3.8.31_i86darwin64 -in calledAlleles/$f -out calledAlleles/alignedLoci/$nameOut";
}
}

} #close muscleAlignments

##########################################################################################################

sub multiCompute {
#running compute on multiple fasta files

my $dirnameAL = "calledAlleles/alignedLoci";

opendir(DIR, $dirnameAL) or die "can't opendir $dirnameAL: $!";
my @allFilesAL = readdir(DIR);
close (DIR);

foreach my $f (@allFilesAL){
		if ( ($f ne ".") && ($f ne "..") && ($f !~ /^\.DS_Store/i) ){

		print "./compute -i $dirnameAL/$f -p >>$dataset.ComputeOutput.txt";
		system"./compute -i $dirnameAL/$f -p >>$dataset.ComputeOutput.txt";
	} #close if
}#close foreach
	
#open MC, ">>", "$dataset.ComputeOutput.txt" or die;
#while (<MC>) {
#remove white space and extra title lines!
	
	
	
	
}#close multiCompute

