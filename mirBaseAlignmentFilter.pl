#!/usr/bin/perl
# Copyright (C) 2014 Ion Torrent Systems, Inc. algntl Rights Reserved
#Line 343 modified by Gabriel Wajnberg on 03/12/2019

# Filter mirBase alignments using various "miRNA" quality and compatibility criterias
# Generate a BAM file with genomic coordinates


use Getopt::Long;
use File::Basename;
use strict;
#use Data::Dumper;
# Params
	# maxMismatches = 1
	# maxIndels = 2
	# maxEdit = maxMismatches + maxIndels
	# max5pSoftClipping = 0
	# max3pSoftClipping = minAdaptorOverlap
	# maxReadLength = 50
	# minAlignLength = 16
my $maxMismatches = 1;
my $maxIndels = 2;
my $maxEdit = 3;
my $max5pSoftClipping = 0;
my $max3pSoftClipping = 5;
my $maxReadLength = 50;
my $minAlignLength = 16;
my $minMQV = 0;

my $padding = 10;

# I/Os
my $bam;
my $mirbase_gff;
#my $output;
my $goutput;
my $filteredOut = "";
my $alignStats = "";
my $mapper = ""; # star, tmap or bowtie2

my $genome = "";

my $samtools = "samtools";
my $bedtools = "/usr/bin/bedtools";
my $mirnaMappingType = ""; # genome | mirbase | highconf

my $plugin_dir;
my $controlPrefix;
my $controlOutput;

my $usage = "Usage: " . basename($0) . " --maxMismatches int[1] --maxIndels int[2] --maxEdit int[3] --max5pSoftClipping int[2] --maxReadLength int[50] --mapper [star|bowtie2|tmap] --padding [10]\n        --bam input_bam_file --mirbase_gff mirBase_gff_file --output bam_output --goutput bam_output_genomic_coord --mirnaMappingType [genome|mirbase]\n\n";

print "Options: \n";
print (join("|",@ARGV),"\n");

GetOptions("maxMismatches=i" => \$maxMismatches,    # numeric
			"maxIndels=i" => \$maxIndels,
			"maxEdit=i" => \$maxEdit,
			"max5pSoftClipping=i" => \$max5pSoftClipping,
			"max3pSoftClipping=i" => \$max3pSoftClipping,
			"maxReadLength=i" => \$maxReadLength,
			"minAlignLength=i" => \$minAlignLength,
			"minMQV=i" => \$minMQV,
			"padding=i" => \$padding,
			"genome=s" => \$genome,
			"mapper=s" => \$mapper,
			"bam=s" => \$bam,
			# "minAdaptorBases=s" => \$minAdaptorBases,
			"mirbase_gff=s" => \$mirbase_gff,
			"filteredOut=s" => \$filteredOut,
			"mirnaMappingType=s" => \$mirnaMappingType,
			"alignStats=s" => \$alignStats,
			#"output=s" => \$output,
			"pluginDir=s" => \$plugin_dir,
			"controlOutput=s" => \$controlOutput,
			"controlPrefix=s" => \$controlPrefix,
			"goutput=s" => \$goutput) or die("Error in command line arguments\n$usage");


my $mirbaseAlign = 0;
if ($mirnaMappingType eq "mirbase" || $mirnaMappingType eq "highconf") {
	$mirbaseAlign = 1;
}

my $output_folder = dirname($bam);

# Parse GFF to get genomic coordinates and translate
#my (%coord, %strand);

print "mirBase mapping type: $mirnaMappingType\n";

my $controlMapped = 0;
if ($controlPrefix ne "") {
	print "Using control sequences\n";
	$controlMapped = 1;
}

my %coord;
if ($mirnaMappingType eq "mirbase") {
	print "Parsing MIRBASE GFF FILE - $mirbase_gff\n";
	open GFF, "<$mirbase_gff" or die "$mirbase_gff: $!";
	while (<GFF>) {
		my @cols = split(/\t/);
		my $strd = $cols[6];
		my $miR_coord;
		if ($strd eq "+") {
			$miR_coord = $cols[3] + $padding;
		} else {
			$miR_coord = $cols[4] - $padding;
		}
		$coord{$cols[2]} = $cols[0] . ":" . $miR_coord; # genomic start or end depending on miR strand
	}
	close GFF;
}

my $readInputBam = "$samtools view -h $bam";
print "samtools command = $readInputBam\n";

my $tmp_output = "mirbase_filteredin_tmp.sam";


print "Parsing and filtering input BAM file\n";

my $filteredOutReadCount = 0;
my $filteredInReadCount = 0;
my $controlInReadCount = 0;

my $maxEditCount = 0;
my $maxMismatchesCount = 0;
my $maxIndelsCount = 0;
my $max5pSoftClippingCount = 0;
my $max3pSoftClippingCount = 0;
my $maxReadLengthCount = 0;
my $minAlignLengthCount = 0;
my $minMQVCount = 0;

my $sam_control = $controlOutput . ".sam";

open BAM, "$readInputBam |" or die $!;
open O, ">$tmp_output" or die $!;
open R, ">$filteredOut" or die $!;

if ($controlMapped == 1) {
	open C, ">$sam_control";
}
while (<BAM>) {
    #print $_,"\n";
    #exit;
	my $algnt = $_;
	my @col = split(/\t/), $algnt;
	my $bitwise = $col[1];
	if ($algnt =~m/^\@/) {
		
		if ($controlMapped == 1) {
			if ($algnt =~m/^\@SQ/) {
				if ($algnt =~m/SN:$controlPrefix/) {
					print C $algnt;
				}
			} else {
				print C $algnt;
			}
		}
			
		if ($algnt =~m/^\@(RG|PG)/) {
	
		    print O $algnt;
		    #print $algnt,"\n";
		} elsif ($mirnaMappingType eq "genome" || $mirnaMappingType eq "highconf") {
		    print O $algnt ;
		    #print $algnt,"\n";
		}
		print R $algnt;
	} elsif ($bitwise & 4) { # unaligned
		print R $algnt;
	} else {	# aligned
		my $filterOut = 0;
		my $seq = $col[0];
		my $ref = $col[2];		
		my $strand = ($bitwise & 16) ? "rev" : "fwd";  # Discard rev mapped reads since miRbase reference is stranded with bedtools getfasta -s
		my $mqv = $col[4];
		my $cigar = $col[5];
		my $seqLength = length($col[9]);
		my $sofclipped5p = 0;
		my $sofclipped3p = 0;
		my $leftSoftClip = 0;
		my $rightSoftClip = 0;
		
		if ($bitwise & 256) { # secondary alignment
			next;
		}
		
		if ($cigar =~m/^(\d+)S/) {
			$leftSoftClip = $1; 
		}
		if ($cigar =~m/(\d+)S$/) {
			$rightSoftClip = $1; 
		}
		if ($mirnaMappingType eq "genome" && $strand eq "rev") {
			$sofclipped5p = $rightSoftClip;
			$sofclipped3p = $leftSoftClip;
		} else {
			$sofclipped5p = $leftSoftClip;
			$sofclipped3p = $rightSoftClip;
		}
		
		my $indels=0;
		while ($cigar =~m/(\d+)([ID])/g) {
			$indels += $1;
		}
		my $matches=0;
		while ($cigar =~m/(\d+)([M])/g) {
			$matches += $1;
		}
		
		$algnt =~ m/\tNM:i:(\d+)/i;	# Edit distance
		my $edit = $1;
				
		# my $nhits=0;
		# if ($algnt =~ m/\tNH:i:(\d+)/i) {	# Bowtie2 and Tmap do not output the NH tag
			# $nhits = $1;
		# }
		# NM Tag:
		# Tmap counts SUBS and INDELS. AN Insertion of 3 bases counts for 3.
		# Bowtie2 counts indels and SUbs
		# STAR only counts SUbs
		my $subs;
		if ($mapper eq "star") {
			# adjust maxEdit
			 #$maxEdit = $maxMismatches;
			$subs = $edit;
			$edit = $indels + $edit;
		} else {
			$subs = $edit - $indels;
		}
		if ($edit > $maxEdit || $subs > $maxMismatches || $indels > $maxIndels || $sofclipped5p > $max5pSoftClipping || $sofclipped3p > $max3pSoftClipping || $seqLength > $maxReadLength || $matches < $minAlignLength || $mqv < $minMQV) {
			$filterOut = 1;
			$maxEditCount++ if $edit > $maxEdit;
			$maxMismatchesCount++ if $subs > $maxMismatches;
			$maxIndelsCount++ if  $indels > $maxIndels;
			$max5pSoftClippingCount++ if  $sofclipped5p > $max5pSoftClipping;
			$max3pSoftClippingCount++ if $sofclipped3p > $max3pSoftClipping ;
			$maxReadLengthCount++ if  $seqLength > $maxReadLength;
			$minAlignLengthCount++ if $matches < $minAlignLength;
			$minMQVCount++ if $mqv < $minMQV;
			
		}	
		
		unless ($filterOut == 1) {
			if ($mirbaseAlign == 1 ) {
				$filterOut = 1 if $strand eq "rev";
			}
		}

		if ($filterOut == 0) {	
			
			if ($controlMapped == 1 && $ref=~m/^$controlPrefix/) {
				$controlInReadCount++;
				print C $algnt;
			} else {
				$filteredInReadCount++;
				print O $algnt;
				#print $algnt, "\n";
			}
		} else {
			$filteredOutReadCount++;
			print R $algnt;;
		}
	}
}
close BAM;
close O;
close R;



if ($mirnaMappingType eq "genome") {
	# NEED TO FILTER OUT READ ON THE WRONG STRAND/KEEP READS ON RIGHT STRAND WITH BEDTOOLS
	# bedtools intersect -abam bam_file -b gff_file -wa -s -header
	my $tmpunsorted = $tmp_output . "_unsorted.bam";
	my $tmpsorted = $tmp_output . "_sorted.bam";
	system("$samtools view -hSb $tmp_output > $tmpunsorted");
	my $gprefix = $tmpsorted;
	$gprefix =~s/\.bam$//;
	system("$samtools sort $tmpunsorted $gprefix");
	system("$bedtools intersect -abam $tmpsorted -b $mirbase_gff -wa -s -header |samtools view -h - > $tmp_output");
	system("rm -f $tmpunsorted $tmpsorted");
}


print "Generating final filtered and genomic converted BAM files\n";
my $sam_goutput = $goutput . ".sam";
open BAM, "<$tmp_output" or die "$tmp_output: $!";
open GENOMEBAM, ">$sam_goutput" or die "$sam_goutput: $!";

#print Dumper (%coord);
if ($mirnaMappingType eq "genome" || $mirnaMappingType eq "mirbase") {
	my $dict = $genome;
	$dict =~s/(fasta|fa)$/dict/;
	if (! -e $dict) {
		$dict = $output_folder . "/" . $dict;
		system("java -jar $plugin_dir/bin/picard-tools-current/CreateSequenceDictionary.jar R=$genome O=${dict}");
	}
	open DICT, "<$dict" or die "$dict: $!";
	print GENOMEBAM <DICT>;
	close DICT;
}


while (<BAM>) {
	chomp;
	my $algnt = $_;
	if ($algnt =~m/^\@/) {
		if ($algnt =~m/^\@(RG|PG)/) {
			print GENOMEBAM $algnt, "\n";
		} elsif ($mirnaMappingType eq "highconf") {
			print GENOMEBAM $algnt, "\n";
		}		
	} else {
		my @col = split(/\t/), $algnt;
		my $seqid = $col[0];
		my $bitwise = $col[1];
		my $ref = $col[2];
		my $chr = $ref;
		my $start = $col[3];
		my $gstart = $start;
		my $cigar = $col[5];
		my $newCigar = $cigar;
		my $dnaSeq = $col[9];
		my $qual = $col[10];
		#print "$seqid\t$bitwise\t$ref\t$chr\t$start\t$gstart\t$cigar\t$newCigar\t$dnaSeq\t$dnaSeq\n";
	#	exit;
		# if (m/\tNH:i:\d+/i) {
			# $algnt =~s/\tNH:i:\d+/\tNH:i:$nhits/;
		# }
		# Remove RG:Z:NOID tag for tmap or qualimap will fail
		#$algnt =~s/\tRG:Z:\w+//;
		
		if ($mirnaMappingType eq "mirbase") {			
		    # Generate the genomic BAM version
		   # print "$coord{$ref} $ref\n";
		    #		    exit;
		    #modified line by Gabriel Wajnberg for Ldec2.0 genome parsing
		    $coord{$ref} =~m/(\w+\.\d+):(\d+)/;
	#	    
	#	    exit;
		    ($chr,$gstart) = ($1,$2);
		   # print "$chr $gstart\n";
		   # exit;
		    $ref =~ m/_(fwd|rev)$/;
			my $ref_s = $1;
			$newCigar = "";
			if ($ref_s eq "rev") {
				# need to change alignment strand, revcomp seq, rev quals and change position
				$bitwise = $bitwise ^ 16;
				# M/I/S/=/X => sum of sequence ; M/D/=/X/N => alignment length
				my $alignLength=0;
				while ($cigar =~m/(\d+)([MD=XN])/g) {
					$alignLength += $1;
				}
				while ($cigar =~m/(\d+[MDIXNSHP=])/g) {
					 $newCigar = $1 .  $newCigar;
				}
				$gstart = $gstart - $start + $padding + 1 - $alignLength + 1;
				$dnaSeq =~tr/ACGTacgt/TGCAtgca/;
				$dnaSeq = reverse $dnaSeq;
				$qual = reverse $qual;
			} else {
				$gstart = $gstart + $start - $padding - 1;
				$newCigar = $cigar;
			}
		}
		@col = split(/\t/), $algnt;
		print GENOMEBAM join("\t", ($col[0], $bitwise, $chr, $gstart, $col[4], $newCigar, @col[6..8], $dnaSeq, $qual, @col[11..$#col])), "\n";
	}
}
close BAM;
#close MIRBAM;
close GENOMEBAM;

if ($controlMapped == 1) {
	my $gprefix = $controlOutput;
	$gprefix =~s/\.bam$//;
	my $unsorted = $gprefix . "unsorted.bam";
	system("$samtools view -Sbh $sam_control > $unsorted");
	#exit;
	system("$samtools sort $unsorted $gprefix");
	system("$samtools index $controlOutput");
	system("rm -f $unsorted $sam_control");
}


open STATS, ">>$alignStats";
print STATS "mirBase Filtered Out Reads: $filteredOutReadCount\n";
print STATS "mirBase Aligned Reads: $filteredInReadCount\n";
print STATS "Control Aligned Reads: $controlInReadCount\n";

print STATS "Filtered by maxEdit ($maxEdit): $maxEditCount\n";
print STATS "Filtered by maxMismatches ($maxMismatches): $maxMismatchesCount\n";
print STATS "Filtered by maxIndels ($maxIndels): $maxIndelsCount\n";
print STATS "Filtered by max5pSoftClipping ($max5pSoftClipping): $max5pSoftClippingCount\n";
print STATS "Filtered by max3pSoftClipping ($max3pSoftClipping): $max3pSoftClippingCount\n";
print STATS "Filtered by maxReadLength ($maxReadLength): $maxReadLengthCount\n";
print STATS "Filtered by minAlignLength ($minAlignLength): $minAlignLengthCount\n";
print STATS "Filtered by minMQV ($minMQV): $minMQVCount\n";

close STATS;


my $gunsorted = $goutput . "_unsorted.bam";
system("$samtools view -Sbh $sam_goutput > $gunsorted");
my $gprefix = $goutput;
$gprefix =~s/\.bam$//;
system("$samtools sort $gunsorted $gprefix");
system("$samtools index $goutput");

#system("rm -f $unsorted $gunsorted $sam_goutput $sam_output $tmp_output $tmp_filteredOut");
system("rm -f $gunsorted $sam_goutput $tmp_output");
#system("rm -f $gunsorted $sam_goutput");
