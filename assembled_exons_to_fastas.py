#!/usr/bin/env python
from os import path
from os import makedirs
from re import sub   #This imports regular expression usage
from optparse import OptionParser    #Imports the option parser module

###### OPTIONS and USAGE ######
parser = OptionParser(usage = """assembled_exons_to_fastas.py -l PLSX_list -f FASTA_file -d OUT_directory

assembled_exons_to_fastas.py -- Processes the BLAT output from several samples 
    comparing reference-guided contig assemblies with the reference (probe) 
    contigs. The output is a directory containing a fasta file for each
    reference contig holding the sequence for each sample.

Copyright (c) 2014 Weitemier et al.
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

Input - A file containing a list of every .pslx file to be processed (one file
    per line), and a fasta file of the targeted contigs or exons.""")
parser.add_option("-l", action="store", type="string", dest="conname",
    help="File containing a list of the .pslx files to process; one name per line", default="")
parser.add_option("-f", action="store", type="string", dest="faname",
    help="Fasta file of targeted contigs", default="")
parser.add_option("-d", action="store", type="string", dest="dirname",
    help="Name of directory for output files", default="Exons_to_be_aligned")
(options, args) = parser.parse_args()

# Makes sure all filenames are given
if options.conname == "":
    parser.error("Please include a .pslx list file using -l.")
if options.faname == "":
    parser.error("Please include a fasta file using -f.")
if path.exists(options.dirname):
    parser.error("The output directory already exists. Please remove it or provide a new name using the -d option.")

###### OPENING INPUT/OUTPUT FILES ######
ConfigFile = open(options.conname, 'r')
FastaFile = open(options.faname, 'r')
makedirs(options.dirname)

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

# Opening the pslx list and grabbing the entries
pslxFiles = []
ConLine = ConfigFile.readline().strip()
while ConLine:
    pslxFiles.append(ConLine)
    ConLine = ConfigFile.readline().strip()
ConfigFile.close()

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
        Fields = pslxLine.split("\t")
        Length = int(Fields[0]) + int(Fields[1])
        ThisExon = Fields[13]
        if not ThisExon in Contigs:
            sys.exit("The contigs in the fasta file don't match those in the .pslx files. Be sure the names of the contigs don't contain spaces.")
        if Name in Contigs[ThisExon]:
            if Length > Contigs[ThisExon][Name][0]:
                Contigs[ThisExon][Name] = [Length, Fields[21].rstrip(',')]
        else:
            Contigs[ThisExon][Name] = [Length, Fields[21].rstrip(',')]
        pslxLine = File.readline().strip()
    File.close()

for exon in Contigs:
    OutExon = sub(',', r'...', exon)
    OutName = options.dirname + '/' + "To_align_" + OutExon + ".fasta"
    OutFile = open(OutName, 'w')
    Filler = 'n' * len(ReferenceContigs[exon])
    for Name in Names:
        if Name in Contigs[exon]:
            OutFile.write(">%s\n%s\n" % (Name, Contigs[exon][Name][1]))
        else:
            OutFile.write(">%s\n%s\n" % (Name, Filler))
    OutFile.close()
#EOF
