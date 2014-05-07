#!/usr/bin/env python
from os import path  #Importing two methods from os module
from optparse import OptionParser    #Imports the option parser module

###### OPTIONS and USAGE #######################################################
parser = OptionParser(usage = """blat_block_analyzer.py -i INFILE -o OUTFILE -l LENGTH -s MIN_BLOCK

blat_block_analyzer.py - Used for grabbing lines in a blat .pslx file matching
    certain criteria of the blocks that form a hit. Outputs a fasta file with
    the sequences of the matching blocks.

Copyright (c) 2012 Kevin Weitemier.
This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. A copy of this license is available at <http://www.gnu.org/licenses/>.
Great effort has been taken to make this software perform its said
task, however, this software comes with ABSOLUTELY NO WARRANTY,
not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Input - A .pslx file from blat.""")
parser.add_option("-i", action="store", type="string", dest="inname",
    help="Input filename", default="")
parser.add_option("-o", action="store", type="string", dest="outname",
    help="Output filename", default="")
parser.add_option("-l", action="store", type="int", dest="length", default=960,
    help="""Minimum length of the *sum* of the blocks greater than MIN_BLOCK.
Default=960""")
parser.add_option("-s", action="store", type="int", dest="min_block",
    default=80, help="""Minimum size for an individual block to be considered.
Default=80""")
(options, args) = parser.parse_args()

# Makes sure all filenames are given
if options.inname == "":
    parser.error("Please include an input file using -i.")
if options.outname == "":
    parser.error("Please include an output file using -o.")

###### OPENING INPUT/OUTPUT FILES #############################################
if path.isfile(options.outname):
    parser.error("""The output filename exists. Please delete it first or \
choose a new name.""")

###### FUNCTIONS ##############################################################

def BuildDictionary(File):
    '''From a .pslx file, returns a dictionary of loci hitting contigs that meet
    user-specified parameters of minimum block size and total block sum.'''
    Loci = {}
    for Line in File: #This puts all the loci in a dictionary with hitting contigs and block sizes
        Line = Line.strip()
        Fields = Line.split('\t')
        LocusName = Fields[9]
        if LocusName not in Loci:
            Loci[LocusName] = []
        SizeLine = Fields[18].rstrip(',')
        Sizes = SizeLine.split(',')
        Contig = str(Fields[13]) + "___" +  Fields[22]
        ContigAndSizes = Sizes
        ContigAndSizes.insert(0, Contig)
        Loci[LocusName].append(ContigAndSizes)
    LociToDelete = []
    for Locus, ContigList in Loci.iteritems(): #Sums all the (qualifying) sizes for this locus
        TotalLength = 0
        for OneContigSizes in ContigList:
            for Size in OneContigSizes[1:]:
                Size = int(Size)
                if Size >= options.min_block:
                    TotalLength = TotalLength + Size
        if TotalLength < options.length:
            LociToDelete.append(Locus)
    for Locus in LociToDelete: #Deletes the Loci that are too small
        del Loci[Locus]
    return Loci

def GrabSequences(Dict):
    '''From a dictionary where keys are loci and values are lists of lists (where
    the first value in each small list is a contig ID, and the other values are
    sequence lengths), returns a dictionary with the same keys but values are 
    the sequences for those lengths above a minimum.'''
    SeqDict = {}
    for Key in Dict.keys():
        SeqDict[Key] = []
    for Line in Infile:
        Line = Line.strip()
        Fields = Line.split('\t')
        LocusName = Fields[9]
        CurrentContig = str(Fields[13]) + "___" + Fields[22]
        if LocusName in Dict: #Does this line contain a locus that is covered?
            ContigList = Dict[LocusName]
            ContigDict = {}
            BlockSeqs = Fields[22].rstrip(',').split(',') #These are sequences for the contigs, where they match the target locus.
            for Contig in ContigList:
                ContigDict[Contig[0]] = Contig[1:] #turns the confusing double list into a dictionary
            if CurrentContig in ContigDict:         
                for Index, Length in enumerate(ContigDict[CurrentContig]):
                   if int(Length) >= options.min_block:
                       #Outfile.write("%s\n%s\n%s\n" % (Line,ContigDict[CurrentContig],BlockSeqs)) 
                       ContigAndBlock = [str(Fields[13]), BlockSeqs[Index]]
                       SeqDict[Fields[9]].append(ContigAndBlock)
    return SeqDict

###### MAIN PROGRAM ############################################################
#Opening files
Infile = open(options.inname, 'r')
Outfile = open(options.outname, 'w')
BigHitLoci = BuildDictionary(Infile)
#for Key, Value in BigHitLoci.iteritems():
#    Outfile.write("%s\t%s\n" % (Key,Value))

Infile.seek(0) #This should send program back to beginning of file.
LociContigSeqs = GrabSequences(BigHitLoci)

for Locus, Seqs in LociContigSeqs.iteritems():
    for Index, Sequence in enumerate(Seqs):
        #Returns the sequence in fasta format, with ID of >Locus_exon#_length
        ID = Sequence[0] + "_" + str(Index+1)
        Outfile.write(">%s,%s,%d\n%s\n" % (Locus,ID,len(Sequence[1]),Sequence[1]))

Infile.close()
Outfile.close()
###### END OF FILE #############################################################
