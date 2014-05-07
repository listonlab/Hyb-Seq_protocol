#!/bin/env python
Infile = open('clustered_bubbles.clstr', 'r')
Outfile = open('singleton_clustered_bubbles.clstr', 'w')

Line = Infile.readline()

while Line:
    if Line.startswith('>'):
        StoredCluster = Line.strip()
        Line = Infile.readline()
        FirstName = Line.strip()
        Line = Infile.readline()
        if Line.startswith('>'):
            Outfile.write("%s\n%s\n" % (StoredCluster,FirstName))
#            continue
        else:
          if not Line:
            break
          else:
            ListOfLines, ListOfPercents = ([], [])
            ListOfLines.append(FirstName)
            while not Line.startswith('>'):
                ListOfLines.append(Line.strip())
                Line = Infile.readline()
            for Name in ListOfLines:
                if 'at' in Name:
                    Fields = Name.split('/')
                    Percent = Fields[2].strip('%')
                    if Percent != '100.00':
                        break
            else:
                Outfile.write("%s\n%s\n" % (StoredCluster,ListOfLines))

    else:
        Outfile.write("Got mis-placed, moving to next cluster!\n")
        while not Line.startswith('>'):
            Line = Infile.readline()
#        continue

Outfile.write("%s\n%s\n" % (StoredCluster,FirstName))

Infile.close()
Outfile.close()


