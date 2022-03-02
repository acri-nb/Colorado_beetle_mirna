# Colorado_beetle_mirna
Pipeline to do the analysis of Colorado potato beetles used for ION PROTON DATA 
For nova seq you still can use the same reference data, but you will need to perform a new analysis with nextflex smallRNAseq lib requirements (trimming and mapping instructins in small RNA seq pipeline).

temp.py is the pipeline program which uses as input a raw basecaller bam file from the IonProton, map to the fasta miRNA of genome Ldec 2.0 version  and produces in the end a aligned reads filtered file which is going to be the input for the feature count software.

the temp.py needs bowtie2, cutadapt , picard and the perl script which is also uploaded to this github (mirBaseAlignmentFilter.pl).

Both the python and perl scripts were modified from the original scrips using by SmallRNA pluggin from IonTorrent Server to be able to do the analysis of the Colorado potato beetle.

the beetles.Rmd is the R markdown file used to build the report.

featurecounts command line used:


featureCounts --minReadOverlap 15 -M -O -d 15 -D 40 -Q 0 -T 12 -a ../new_lde.count.gff3 -t miRNA -g Name -o featureCount.mirbase.total.txt  IonXpressRNA_013_filtered.bam IonXpressRNA_014_filtered.bam IonXpressRNA_015_filtered.bam IonXpressRNA_016_filtered.bam IonXpressRNA_017_filtered.bam IonXpressRNA_018_filtered.bam

additional files location:

acri@10.111.90.44:/mnt/user/LTS//backup2/Reference_data/References_genome/Leptinotarsa_decemlineata/mirna_analysis

there you can find the bowtie2 indexes for the mapping, the gtf and fasta for the mirBaseAlignmentFilter.pl and the gff3 for the featurecounts.

