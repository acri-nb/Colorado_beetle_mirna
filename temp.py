#!/usr/bin/env python
import os
import sys
import time
import commands


def toInt(i):
        return int(re.sub(',','',str(i)))

def printlog(msg):
        sys.stderr.write(msg)
        sys.stderr.write('\n')
        sys.stderr.flush()

def printtime(msg):
        printlog( '(%s) %s'%(time.strftime('%X'),msg) )
def runcmd( cmd, log, fatal=True ):
        printtime(cmd)
        retval = os.system(cmd)
        if retval != 0:
                sys.stderr.write( "$ "+cmd+"\n" )
                sys.stderr.write( "ERROR: Failed running command (status = %d). See '%s'.\n" % (retval,os.path.basename(log)) )
                if fatal: sys.exit(1)

output_dir = sys.argv[0]
run_log = sys.argv[1]
genome_reference = sys.argv[3]
runlog = os.path.join(output_dir,run_log)
cat = '/bin/cat'
mirbase_gff= "new_lde.gff"
padding=10
wg_mirbase_gff="new_lde_nameChanged.gtf"
countGff="new_lde.count.gff3"
bedtools = '/usr/bin/bedtools'
mature_reference= "new_lde.fasta"

#cmd = '%s %s |perl -lane \'@l = split; if ($l[2] eq "gene") {$l[8]=~m/Name=([^;]+)/; $pt=$1;} elsif ($l[2] eq "miRNA") {$l[8] =~m/Name=([^;]+)/;$name=$1; s/Name=[^;]+/Name=${pt}_$name/; print}\' > %s'%(cat,mirbase_gff,countGff)
#runcmd( cmd, runlog )

#cmd = '%s %s |perl -lane \'@l = split; $type=$l[2]; $strd = $l[6]; $strd_str = $strd eq "+" ? "fwd": "rev"; if ($type eq "miRNA") {$start = $l[3] - %s;$end= $l[4] + %s; $mir=$l[8]; $mir =~m/ID=([^;]+)/; $id = $1; $mir =~m/Name=([^;]+)/; $name = $1; $new_name = "${name}_${id}_${strd_str}"; print join "\t",($l[0], $l[1], $new_name, $start, $end, $l[5],$l[6],$l[7],$l[8])}\' > %s'%(cat,mirbase_gff,padding,padding,wg_mirbase_gff)
#runcmd( cmd, runlog )
#Cmd = '%S Getfasta -S -Name -Fi %s -bed %s -fo %s'%(bedtools,genome_reference,wg_mirbase_gff,mature_reference)
#runcmd( cmd, runlog )

bamfile= "temp_folder/Ion_Convert_018_rawlib.basecaller.bam"
minAdaptorBases=30
withAdaptor="project02122019/IonXpressRNA_018_whithadptor.sam"
adaptorFilterLog="adaptor.log"
cmd = 'samtools view -h %s | perl -lane \'if (m/^\@/) {print} elsif (m/\s+ZB:i:(\d+)/) {$zb=$1;print if $zb >= %s}\' > %s 2> %s' %(bamfile,minAdaptorBases, withAdaptor, adaptorFilterLog)
runcmd(cmd, runlog )
withoutAdaptor = "project02122019/IonXpressRNA_018_whithoutAdaptor.sam"
cmd = 'samtools view -h %s | perl -lane \'if (m/^\@/) {print} elsif (m/\s+ZB:i:(\d+)/) {$zb=$1;print if $zb < %s} elsif (! m/\s+ZB:i:(\d+)/) {print}\' > %s 2>> %s' % (bamfile,minAdaptorBases, withoutAdaptor, adaptorFilterLog)
runcmd(cmd, runlog )
cmd ='java -Xmx8g -jar /opt/picard/picard-tools-current/picard.jar SamToFastq I=project02122019/IonXpressRNA_018_whithadptor.sam F=project02122019/IonXpressRNA_018_whithadptor.fastq'
runcmd(cmd, runlog )
cmd='/results/plugins/smallRNA/bin/cutadapt-1.8.1/bin/cutadapt -a ATCACCGACTGCCCATAGAG -m 015 --overlap 5 project02122019/IonXpressRNA_018_whithadptor.fastq > project02122019/IonXpressRNA_018_whithadptor_cutadapt.fastq'
runcmd(cmd, runlog )
cmd='/results/plugins/smallRNA/bin/bowtie2-2.2.5/bowtie2-align-s --wrapper basic-0 -p 12 --end-to-end --norc --very-sensitive -x /results/plugins/smallRNA/annotation/Ldec2.0/bowtie2/new_lde -U project02122019/IonXpressRNA_018_whithadptor_cutadapt.fastq | samtools  view -bSh - > project02122019/IonXpressRNA_018.bam'
runcmd(cmd, runlog )
cmd='mirBaseAlignmentFilter.pl --maxMismatches 1 --maxIndels 1 --maxEdit 2 --max5pSoftClipping 0 --max3pSoftClipping 3 --maxReadLength 35 --minAlignLength 16 --mapper bowtie2 --bam project02122019/IonXpressRNA_018.bam --mirbase_gff new_lde_nameChanged.gtf --goutput project02122019/IonXpressRNA_018_filtered.bam --padding 10 --genome /results/referenceLibrary/tmap-f3/Ldec2.0/Ldec2.0.fasta --filteredOut project02122019/IonXpressRNA_018.filtered_out.sam --alignStats project02122019/IonXpressRNA_018.filterstats.txt --mirnaMappingType mirbase --minMQV 0 --pluginDir /results/plugins/smallRNA --controlPrefix ControlSeq_ --controlOutput project02122019/IonXpressRNA_018.control.bam > project02122019/IonXpressRNA_018AlignmentFilter.log 2> project02122019/IonXpressRNA_018mirbasealignment.log';
runcmd(cmd, runlog )
