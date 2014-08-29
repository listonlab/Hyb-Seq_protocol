Hyb-Seq: Combining target enrichment and genome skimming for plant phylogenomics
================================================================================

Programs and documentation accompanying the paper:
[K. Weitemier, S.C.K. Straub, R. Cronn, M. Fishbein, A. McDonnell, R. Schmickl,
and A. Liston. 2014. Hyb-Seq: Combining target enrichment and genome skimming for plant phylogenomics
Applications in Plant Sciences 2(9): 1400042](http://www.bioone.org/doi/pdf/10.3732/apps.1400042).

A protocol for developing probes for targeted sequence capture from an input of
genomic and transcriptomic data is contained within Data Supplement S1:
[Building_exon_probes.sh](Building_exon_probes.sh).
This file is both a program that can be run in a Linux environment, and a
helpful document with detailed notes.

A workflow for analyzing the raw reads obtained from the Hyb-Seq protocol is
described in [Data Supplement S2](Data_Supplement_S2.pdf).

[assembled_exons_to_fastas](assembled_exons_to_fastas), [blat_block_analyzer](blat_block_analyzer), and [grab_singleton_clusters](grab_singleton_clusters) are
all companion programs integrated into either the probe development or read
analysis pipelines.

Example datasets have been added: [genome.fasta](genome.fasta) and [transcriptome.fasta](transcriptome.fasta).

NOTE: The latest version of CD-HIT (4.6.1) contains a bug that appears when
using the example data provided here. To successfully use the example data one
should either use a CD-HIT version prior to 4.6 or use the patched version
supplied [here](cd-hit-v4.6.1_fix_max_sequences).

Copyright (c) 2014
K. Weitemier, S.C.K. Straub, R. Cronn, M. Fishbein, A. McDonnell, R. Schmickl, and A. Liston

