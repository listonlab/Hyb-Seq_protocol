#/bin/tcsh

#Disclaimer:
#Although these commands should function with just the input of two fasta
#files, genome.fasta and transcriptome.fasta, their proper execution is not
#guaranteed. Given the idiosyncrasies in file formats and operating
#environments, these commands are meant more as a starting point, to be
#modified as needed by the user. For example, the presence of spaces or tabs in
#the ID line of the fasta files may cause problems downstream.

#Copyright 2014, Weitemier et al.
#This script is free software: you can redistribute it and/or modify it under
#the terms of the GNU General Public License as published by the Free Software
#Foundation, either version 3 of the License, or (at your option) any later
#version. A copy of this license is available at <http://www.gnu.org/licenses/>.
#Great effort has been taken to make this software perform its said
#task, however, this software comes with ABSOLUTELY NO WARRANTY,
#not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#Match genome and transcriptome sequences
#Fasta files should be formatted so that a sequence is not split over lines.
blat genome.fasta transcriptome.fasta -tileSize=8 -minIdentity=99 -noHead -out=pslx genome_v_transcriptome.pslx

#Find and extract transcriptome sequences with only one match against the genome
cut -f10 genome_v_transcriptome.pslx | sort | uniq -c | grep '      1 ' | sed -e 's/      1 /\\\</' -e 's/$/\\\>/' > single_hits_vs_genome.txt
grep -f single_hits_vs_genome.txt genome_v_transcriptome.pslx > single_hit_genome_v_transcriptome.pslx

#Keep only matches hitting at least 95% of the transcript.
awk '{if (($13-$12) >= ($11 * 0.95)) print $0}' single_hit_genome_v_transcriptome.pslx > whole_gene_single_hit_genome_v_transcriptome.pslx

#Extract these transcripts from the transcriptome assembly
cut -f10 whole_gene_single_hit_genome_v_transcriptome.pslx > single_whole_hits_vs_genome.txt
sed -i -e 's/^/>/' single_whole_hits_vs_genome.txt
grep -A1 --no-group-separator -f single_whole_hits_vs_genome.txt transcriptome.fasta > single_transcript_hits.fasta

#Cluster any nested transcripts
cd-hit-est -i single_transcript_hits.fasta -o single_transcript_hits_cluster_100.fasta -c 1.0 -p 1

#Remove transcripts with 90% or greater similarity
cd-hit-est -i single_transcript_hits_cluster_100.fasta -o single_transcript_hits_cluster_90.fasta -c 0.9 -p 1 -g 1
grab_singleton_clusters.py -i single_transcript_hits_cluster_90.fasta.clstr -o unique_single_transcript_hits_cluster_90.fasta.clstr
grep -v '>Cluster' unique_single_transcript_hits_cluster_90.fasta.clstr | cut -d' ' -f2 | sed 's/\.\.\.//' > unique_single_transcript_hits
grep -A1 --no-group-separator -f unique_single_transcript_hits single_transcript_hits_cluster_100.fasta > unique_single_transcript_hits.fasta
sed -i -e 's/>//' unique_single_transcript_hits
grep -f unique_single_transcript_hits whole_gene_single_hit_genome_v_transcriptome.pslx > unique_single_hits.pslx

#Find loci and exons that meet length requirements
blat_block_analyzer.py -i unique_single_hits.pslx -o large_enough_unique_single_hits.fasta -l 960 -s 120

#Remove blocks with 90% or greater similarity
cd-hit-est -i large_enough_unique_single_hits.fasta -o large_enough_unique_single_hits_cluster90.fasta -c 0.9 -p 1 -g 1
grab_singleton_clusters.py -i large_enough_unique_single_hits_cluster90.fasta.clstr -o unique_blocks_large_single_hits_cluster90.fasta.clstr
grep -v '>Cluster' unique_blocks_large_single_hits_cluster90.fasta.clstr | cut -d' ' -f2 | sed 's/\.\.\.//' > unique_blocks_large_single_hits
grep -A1 --no-group-separator -f unique_blocks_large_single_hits large_enough_unique_single_hits.fasta > blocks_for_probe_design.fasta
