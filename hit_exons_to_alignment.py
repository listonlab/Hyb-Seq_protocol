#!/usr/bin/env python
from os import remove
from os import path  #Importing two methods from os module
from re import sub   #This imports regular expression usage
from optparse import OptionParser    #Imports the option parser module

###### OPTIONS and USAGE ######
parser = OptionParser(usage = """hit_exons_to_alignment.py -c CONFIG_file -f FASTA_file

hit_exons_to_alignment.py -- Processes the BLAT output from several samples 
    comparing reference-guided contig assemblies with the reference (probe) 
    contigs to produce a sequence alignment across the samples for each contig,
    suitable for phylognetic analysis.

Copyright (c) 2014 Kevin Weitemier.
Version 0.01

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details. You
should have received a copy of the GNU General Public License along with this
program.  If not, see <http://www.gnu.org/licenses/>.

Input - A config file containing the names of every .pslx file to be processed
   (one file per line), and a fasta file of the targeted contigs or exons.""")
parser.add_option("-c", action="store", type="string", dest="conname",
    help="Config filename", default="")
parser.add_option("-f", action="store", type="string", dest="faname",
    help="Contig fasta filename", default="")
(options, args) = parser.parse_args()

# Makes sure all filenames are given
if options.conname == "":
    parser.error("Please include a config file using -c.")
if options.faname == "":
    parser.error("Please include a fasta file using -f.")

###### DEFINING FUNCTIONS ######

###### OPENING INPUT/OUTPUT FILES ######
ConfigFile = open(options.conname, 'r')
FastaFile = open(options.faname, 'r')

# Opening the contigs file and makeing a dictionary of each sequence
ReferenceContigs = {}
FaLine = FastaFile.readline().strip()
while FaLine:
    if not FaLine.startswith('>'):
        sys.exit("""The fasta file is not formatted correctly. Please be sure that each entry is given on exactly two lines: the ID line and the sequence line.""")
    ID = FaLine.lstrip('>')
    FaLine = FastaFile.readline().strip()
    Seq = FaLine
    ReferenceContigs[ID] = Seq
    FaLine = FastaFile.readline().strip()

FastaFile.close()

Contigs = {}
for exon in ReferenceContigs:
    Contigs[exon] = {}

print len(Contigs)


# Opening the pslx list and grabbing the entries
pslxFiles = []
ConLine = ConfigFile.readline().strip()
while ConLine:
    pslxFiles.append(ConLine)
    ConLine = ConfigFile.readline().strip()

ConfigFile.close()
print len(pslxFiles)

Names = []

# Processing each pslx file
for Filename in pslxFiles:
    Name = sub('Final_Assembly_',r'',Filename)
    Name = sub('.pslx',r'',Name)
    Names.append(Name)
    
    File = open(Filename, 'r')
    pslxLine = File.readline().strip()
    while pslxLine:
        if not ',' in pslxLine:
            pslxLine = File.readline().strip()
            continue
        print pslxLine
        pslxLine = File.readline().strip()

    File.close()
