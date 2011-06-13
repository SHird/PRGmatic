PRGmatic [v.1] README	
Copyright (C)  2011  Sarah M. Hird 
    Permission is granted to copy, distribute and/or modify this document
    under the terms of the GNU Free Documentation License, Version 1.3
    or any later version published by the Free Software Foundation;
    with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.
    A copy of the license is included in the document entitled "GNU
    Free Documentation License".

1. INSTALLATION
2. EASY INSTALL
3. MORE DETAILED INSTALL
	a. BWA
	b. CAP3
	c. MUSCLE
	d. SAMTOOLS
	e. VARSCAN
	f. COMPUTE
4. TO TURN OFF THE COMPUTE ANALYSIS/MUSCLE ALIGNMENTS
5. PRGmatic FOLDER
6. GENERATING FASTA FILES (RDP)
7. TO RUN bigFastaRename.pl
8. TO RUN PRGmatic
9. INPUT PARAMETERS
10. OUTPUT
11. TEST DATA
12. CONTACT
13. RECOMMENDED READING
14. CITATIONS


1. INSTALLATION
The pipeline is dependent on several other pieces of software. These are in the DependentSoftware folder. The packages have been included with PRGMATIC, but if you have difficulties installing them, please see the original webpages. The pipeline was written on MacOSX and requires that the Developer Tools be installed on your machine. These are on the supplementary disc that comes with Macs.

2. EASY INSTALL
1. Open a terminal window
2. cd into the DependentSoftware folder. 
3. Type “chmod 755 setup.pl” 
4. Type  “./setup.pl”. 
You should see some output to the screen and when the script has finished, there should be three executables (bwa, cap3, samtools) and a script (samtools.pl) in the PRGMATIC folder. 

3. MORE DETAILED INSTALL:
a. BWA:
Included is MacOSX version, which should also work on some other platforms as well. To install, double click bwa-0.5.8a.tar.bz2. Open a terminal window and cd into the bwa-0.5.8a folder. Type “make” (without quotation marks). When this finishes, make a copy of the bwa executable and put the copy in the PRGMATIC folder. This is correctly installed if you double click the executable and a screen with a parameter list opens.
http://sourceforge.net/projects/bio-bwa/files/bwa-0.5.8a.tar.bz2/download
http://bio-bwa.sourceforge.net/

b. CAP3:
Included is the 64-bit, MacOSX version of this program. To install CAP3, double-click cap3.macosx.intel64.tar. When it unpacks into a folder, make a copy of the cap3 executable and put the copy in the PRGMATIC folder. This is correctly installed if you double click the executable and a screen with a parameter list opens.
http://seq.cs.iastate.edu/cap3.html.

c. MUSCLE
Muscle is the sequence aligner for after the loci have been called. It is executable from the file included with PRGMATIC and shouldn’t require anything to make it run as long as it’s in the correct folder. You can double check that it is correctly installed if you type “chmod 755 muscle3.8.31_i86darwin64” then “./muscle3.8.31_i86darwin64” from a terminal window inside the PRGMATIC directory. The included version is for 64-bit Intel MacOSX machines – if you have different hardware or need to get an original copy, see http://www.drive5.com/muscle/.

d. SAMTOOLS:
Included is the MacOSX version, which should also work on some other platforms as well. To install, double click samtools-0.1.8.tar.bz2.  Open a terminal window and cd into samtools-0.1.8. Type “make” (without quotation marks). When this finishes, make a copy of the samtools executable and put the copy in the PRGMATIC folder. This is correctly installed if you double click the executable and a screen with a parameter list opens.  Also make a copy of samtools.pl from the misc folder and put that copy in the PRGMATIC folder. (samtools.pl is not an executable, so as long as it is in the correct folder, it is correctly “installed”.)
http://sourceforge.net/projects/samtools/

e. VARSCAN
The included VarScan should already be installed upon downloading the PRGMATIC package.  You can check this by opening a terminal, cd into PRGMATIC and type “java –jar VarScan.v2.2.2.jar”. This is correctly installed if you see a screen with a parameter list.
http://varscan.sourceforge.net/

f. COMPUTE
The compute package is more difficult to install than the other software. I recommend following the instructions at http://molpopgen.org/software/libsequence.html to get the necessary library (libsequence) and the instructions at http://molpopgen.org/software/lseqsoftware.html to get the analysis package (which contains compute). Once you have this installed, make a copy of the compute executable and put it in the PRGMATIC folder. This is optional software, so if it is not installed on your machine, PRGMATIC will still run, you just need to turn off the part of PRGMATIC.pl that calls compute. 

4. TO TURN OFF THE COMPUTE ANALYSIS OR MUSCLE ALIGNMENT
Open PRGMATIC.pl in a text editor. On line 74, it should say “multiCompute();”. Put a pound sign (#) in front of this line and save the file. That should effectively turn off calling the compute package and there should be no errors. To turn off the MUSCLE multi-sequence alignment, open the PRGmatic.pl in a text editor and put a pound sign (#) in from of line 
73, which should read “muscleAlignments()”. Save the file.
5. PRGMATIC FOLDER
In the PRGMATIC folder you should now have 10 things:

      1	 BWA executable
      2	 calledAlleles folder
      3	 CAP3 executable
      4	 COMPUTE executable (optional)
      5	 DependentSoftware folder
      6	 inputFASTA folder (with alignedLoci folder inside)
      7	 MUSCLE3.8.31_I86DARWIN64
      8	 PRGMATIC.pl
      9	 SAMTOOLS executable
      10 SAMTOOLS.PL
      11 VARSCAN.V2.2.2.JAR

Inside the inputFASTA folder you should place your tag separated fasta files. 

6. (MY PREFERRED METHOD FOR) GENERATING FASTA FILES – USING RDP WEBSITE
Off the 454 machine, you should have gotten at least one .fna file, one .qual file and a folder of .sff files.  To quality control the reads, I use the Ribosomal Database Project’s Pyrosequencing Pipeline (at http://pyro.cme.msu.edu/ ). Their “Pipeline Initial Process” is easy to use, fast and on its own server, so there’s nothing to download. 

RDP Pyrosequencing Pipeline Initial Process Parameters:
Sequence file in FASTA format: (upload the .fna file here)
Quality file in FASTA format (optional): (upload the .qual file here)
Upload a tag file: (upload your tag file here – this is the file that says which tag sequence belongs to which individual. The format is very easy; on a new line for each individual: TagSequence (tab) IndividualName. **When you run the bigFastaRename.pl script to rename the files, the tag file will need to have UNIX line breaks.**)
Gene name: Other
Forward Primers: (paste your forward primer here)
Reverse Primers: (paste your reverse primer here)
*FILTERS*
Forward primer max edit distance (0 to 2): (this refers to how many errors you’ll allow in your forward primer sequence. I use 2 but if you’re being conservative, 0 or 1 will weed out more sequences.)
Reverse primer max edit distance (0 to 2): (same as above but for reverse primer. Again, I use 2.)
Max number of N’s: 0 (I highly recommend using 0 here, since one ambiguous base can be indicative of error prone sequence. See Huse et al. (2007) for a good overview of 454 generated errors).
Min sequence length (>=50): (I use either 100 or 150, depending on the dataset, but you can use whatever you deem appropriate. Shorter sequences are more error prone, but throwing out reads unnecessarily is suboptimal).
Minimum Average Exp Quality Score: (20 is what I usually use. This corresponds to an average error rate of no more than 1/100 bases having an error (DOUBLE CHECK THIS). Higher number here results in fewer reads passing the filter with a higher quality score.)
Keep primers: (Do not check)

Click “Perform Initial Processing”

This is generally very fast and I have my file downloading within 10 minutes (usually). It may take longer for big files. Once the file downloads, double click it and it will unpack a folder with your sequences separated by the names in the tag file. The files in these folders named “Individual_trimmed.fasta” and “Individual_trimmed.qual” can be used as input for the pipeline. I’ve included a script in the DependentSoftware folder (“bigFastaRename.pl”) that renames the sequences in these folders (and their associated quality scores) as IndividualName_00001, IndividualName_00002, etc. and puts them in files called IndividualName.fasta and IndividualName.qual.  Renaming the sequences makes them easier to view and understand later, when knowing which reads came from which individual is helpful. I highly recommend running bigFastaRename.pl but it is not necessary.

7. TO RUN BIGFASTARENAME.PL:
1. Copy or drag the file into the RDP downloaded folder. 
2. Open a terminal window and cd into this folder. 
3. Type “chmod 755 bigFastaRename.pl” to give yourself permission to run the script. 
4. Type “./bigFastaRename.pl”. 
5. You should see the prompt “Drag tagfile here:” on the screen. Drag and drop the tagfile there (it doesn’t need to be in the same folder as everything else) and the script should run through all the folders in the tagfile, outputting a .fasta and a .qual file for each one. IF THIS DOES NOT HAPPEN – check the line breaks on your tag file. They should be Unix or the script will just read the first line of the tagfile and quit.
6. Copy or drag these files into the “inputFASTA” folder in the PRGMATIC folder. You’re ready to go!
*****If the script runs for the first individual but stops after that, the line breaks in the tag file need to be changed to UNIX. The program is reading the first line then hitting a hard return it doesn’t understand and stopping. Changing the line breaks will fix this. *****

8. TO RUN PRGMATIC:
Make sure your individual .fasta files are in the inputFASTA folder (quality files are optional, but should be in this folder if you have them). Also, the quality files need to be of the “same name” as the .fasta file they’re associated with (e.g. IND01.fasta should have a quality file named IND01.qual – this is a requirement of cap3).
Make sure you have the 11 things listed above in working order in the PRGMATIC folder. 
Open a terminal window. cd to the PRGMATIC folder. Type “chmod 755 PRGMATIC.pl” to give yourself permission to run the script. Type “./PRGMATIC.pl” to run the script.

9. INPUT PARAMETERS:
When the program runs, you’ll see a variety of prompts and you need to give the script some information.

“Enter the dataset nickname: ” (don’t use an underscore in the nickname _)  
What you enter here doesn’t have an effect on how the program runs, it’s just a way of giving a name to the various output files that the pipeline generates. I generally use something informative and short. Like Trial0915 for a test run that I did on 15 September. It could also be the focal taxa or locality or whatever.

Parameter Settings
“Minimum number of reads to call an allele (5):” 
To generate the pseudo-reference genome, the pipeline calls “high confidence alleles” from clusters within an individual. This parameter sets the minimum number of reads you want for a cluster to be called an allele. If you’d like to use the default, type 5. If you’d like a more conservative p-rg, enter a higher number, if you’d like to be more liberal, enter a lower number.

“Minimum % identity to call a locus (90):”
To generate the p-rg, the pipeline clusters the high-confidence alleles at a given percent identity (similarity). 90% seems to work pretty well. 

“Minimum coverage for calling consensus in an individual (3):”
Once all the reads have been blasted to the p-rg, VarScan calls a consensus sequence for each locus. Here you can set the minimum number of reads you’d like to call the consensus sequence. Higher values = more conservative. Lower values = more liberal.

“Minimum coverage for calling a SNP (3):”
VarScan also calls SNPs from the reads that blast to each locus. Here you can set the minimum number of reads that support a SNP for it to be called an actual variant. I wouldn’t go below 2, since you’d then be calling every SNP (aka every error) a real variant. 

“Minimum % of reads for calling a SNP (20):”
In addition to a minimum raw number of reads being necessary for a SNP to be called, those reads must also represent a certain percentage of the coverage at that base. 20% seems pretty standard in the literature. Higher numbers will be more conservative and lower numbers will be more liberal.

10. OUTPUT
When the pipeline is finished running, there should be .fasta files in the calledAlleles folder corresponding to the loci in the p-rg. These should be viewable in any program that reads .fasta files.
The “.counts.txt” file displays how many individuals were called for each locus. This can be opened in Excel and then sorted to show which loci contained all or most individuals (i.e. those of highest interest). 
The “.HoHe.txt” file contains how many individuals, heterozygotes and alleles were called for each locus as well as the observed and expected heterozygosities. This file can also be opened directly by Excel.
The “.ComputeOutput.txt” contains the output from the compute analyses. This can also be opened in Excel.
The “MultiHitLoci.txt” contains information on where there are more than two base pairs for a given position on the reference genome within an individual. This is highly suspect, especially if more than two base pairs occur at high frequency (i.e., more than one of a given base). If more than one individual appears on the list for a single locus, that locus is almost definitely paralogous. It might even be wise to throw out every locus on the list, or rerun the analysis with higher cutoffs for calling a locus.

11. TEST DATA
To run the test data, copy the three .fasta and three .qual files into the inputFASTA folder. Run PRGMATIC.pl. There should be a lot of output to the screen. In the zip file with the test data, I put the .counts.txt, .HoHe.txt, and .ComputeOutput.txt files that were generated on my machine. They should match the files you get.

12. CONTACT
Please feel free to contact me about any issues you’re having with PRGMATIC or the dependent software. I’d be more than happy to do what I can – 
Sarah Hird
shird1@tigers.lsu.edu 

13. RECOMMENDED READING
Huse SM, Huber JA, Morrison HG, Sogin ML and DM Welch. 2007. Accuracy and quality of massively parallel DNA pyrosequencing. Genome Biology 8: R143 (doi: 10.1186/gb-2007-8-7-r143)

14. CITATIONS
If you use PRGMATIC, please cite (all of) the following papers:
HIRD, S.M., BRUMFIELD, R.T. & CARSTENS, B.C. 2011. PRGmatic: an efficient pipeline for collating genome-reduced second generations sequencing data using a pseudo reference genome. Molecular Ecology Resources (In Review).
EDGAR, R. C. 2004. MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Research, 32, 6.
HUANG, X. & MADAN, A. 1999. CAP3: A DNA sequence assembly program. Genome Research, 9, 10.
KOBOLDT, D., CHEN, K., WYLIE, T., LARSON, D., MCLELLAN, M., MARDIS, E. R., WEINSTOCK, G., WILSON, R. K. & DING, L. 2009. VarScan: variant detection in massively parallel sequencing of individual and pooled samples. Bioinformatics, 25, 3.
LI, H. & DURBIN, R. 2009. Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics, 25, 7.
LI, H., HANDSAKER, B., WYSOKER, A., FENNELL, T., RUAN, J., HOMER, N., MARTH, G., ABECASIS, G., DURBIN, R. & GROUP, G. P. D. P. 2009. The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25, 2.

Optional program citations:
Compute:
THORNTON, K. 2003. libsequence: a C++ class library for evolutionary genetic analysis. Bioinformatics, 19, 3.
RDP:
COLE, J. R., WANG, Q., CARDENAS, E., FISH, J., CHAI, B., FARRIS, R. J., KULAM-SYED-MOHIDEEN, A. S., MCGARRELL, D. M., MARSH, T., GARRITY, G. M. & TIEDJE, J. M. 2009. The Ribosomal Database Project: improved alignments and new tools for rRNA analysis. Nucleic Acids Research, 37, 5.
Tablet:
MILNE, I., BAYER, M., CARDLE, L., SHAW, P., STEPHEN, G., WRIGHT, F. & MARSHALL, D. 2009. Tablet - next generation sequence assembly visualization. Bioinformatics, 26, 2.


Updated: 7 January 2011 by S. Hird

