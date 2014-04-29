#/bin/tcsh

ls -1 *.pslx | sed -e "s/Final_Assembly_//" -e "s/.pslx//" > file_list.temp
grep "Contig" *.pslx > combined.PSLX
sed -i -e "s/:/\t/" -e "s/Final_Assembly_//" -e "s/.pslx//" *.PSLX
awk 'BEGIN { OFS = "\t" } ; { print $23,$2+$3,$15,$1 }' *PSLX | sort -k4,4 -k3,3 -k2,2nr | sed 's/,//' > combined_nucleotides.sort
uniq -f2 combined_nucleotides.sort > combined_nucleotides.maxmatch
foreach sample (` cat file_list.temp `)
  grep $sample combined_nucleotides.maxmatch > $sample.maxmatch.txt
  awk '{print $3,"\t",$1}' $sample.maxmatch.txt > $sample.maxmatch.txt.tab
  end
date

