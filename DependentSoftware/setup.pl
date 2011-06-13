#!/usr/bin/perl
use strict;
use warnings;
#setup.pl  -  installs PRGmatic 

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



###unpack samtools###
system"bzip2 -cd samtools-0.1.8.tar.bz2 | tar xvf -";
print "\n\nUNZIPPED SAMTOOLS\n\n";
chdir "samtools-0.1.8";
system"pwd";

system" make ";
chdir"..";
system"cp samtools-0.1.8/samtools ../ ";
system"cp samtools-0.1.8/misc/samtools.pl  ../";


###unpack cap3###
system"tar xvf cap3.macosx.intel64.tar";
system"cp cap3.macosx.intel64/cap3 ../ ";

###unpack bwa###
system"bzip2 -cd bwa-0.5.8a.tar.bz2 | tar xvf -";
print "\n\nUNZIPPED BWA \n\n";
chdir "bwa-0.5.8a";
system"pwd";

system" make ";
chdir"..";
system"cp bwa-0.5.8a/bwa ../ ";