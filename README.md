# Colorado_beetle_mirna
Pipeline to do the analysis of Colorado potato beetles

temp.py is the pipeline program which uses as input a raw basecaller bam file from the IonProton, map to the fasta miRNA of genome Ldec 2.0 version  and produces in the end a aligned reads filtered file which is going to be the input for the feature count software.

the temp.py needs bowtie2, cutadapt , picard and the perl script which is also uploaded to this github. 
