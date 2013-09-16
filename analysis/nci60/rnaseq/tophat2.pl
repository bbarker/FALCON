#!/usr/bin/perl

# Xuepeng Sun    09/14/2013

use warnings;
use File::Find;

$tophat=q(export PATH=$PATH:~/software/tophat-2.0.8.Linux_x86_64/);
$samtools=q (export PATH=$PATH:~/software/samtools-0.1.19/);
$cuff=q(export PATH=$PATH:~/software/cufflinks-2.0.2.Linux_x86_64/);
$index=q(~/workdir/human/hg19);
system("$tophat");
system("$samtools");
system("$cuff");

my @file=glob("*");

foreach (@file){
	if (/fastq/){
		my $f=$_;
		my $out=$f."_tophat";
		$cmd="tophat2 -o $out -N 3 --read-gap-length 3 --read-edit-dist 3 -p 25 -g 20 --no-coverage-search $index $f";
		system ("$cmd");
	}
}

$doublepostive="1896_6939_6132_N_doublepostive_ATCACG_R1.fastq/accepted_hits.bam";
$doublenegative="1896_6939_6133_N_doublenegative_CGATGT_R1.fastq/accepted_hits.bam";
$control1="1896_6939_6134_N_K562_Control_TTAGGC_R1.fastq/accepted_hits.bam";
$control2="1896_6939_6136_N_MB-231_Control_ACAGTG_R1.fastq/accepted_hits.bam";
$treat1="1896_6939_6135_N_K562_TM_treated_TGACCA_R1.fastq/accepted_hits.bam";
$treat2="1896_6939_6137_N_MB_231_TM_treated_GCCAAT_R1.fastq/accepted_hits.bam";

$gtf="/home/xuepeng/workdir/human/hg19/hg19.gtf";

$cuff1="cuffdiff -o pos_neg --FDR 0.01 -u -p 25 -L positive,negtive $gtf $doublepostive $doublenegative";
$cuff2="cuffdiff -o control_treatment --FDR 0.01 -u -p 25 -L control,treatment $gtf $control1,$control2 $treat1,$treat2";

system("$cuff1");
system("$cuff2");

