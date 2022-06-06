#!/usr/bin/perl

## run bwa aln and samse on fastq files passed from the commandline

## USAGE: runbwa.pl /fullpath/to/files/*fastq

use warnings;

unless(@ARGV){
    die "runbwa.pl /fullpath/to/files/*fastq";
}

$ctr=0;
unless(-e 'sam_sai'){
    mkdir 'sam_sai', 0755 or die "Failed to make sam_sai directory";
}
foreach $fastq (@ARGV){
	$fastq =~ m/([-\w+\.]+)\.fastq$/; ## fragile regexp to catch base name of individuals
	$id = $1;
	print "Mapping reads for $id\n";
	## Use mem instead of aln as suggested in the manual (mem does NOT work for reads<70bp)
	## mem didn't work as expected, went back to aln (look in ../bin/runbwa.pl)
	system "bwa aln -l 20 -k 2 -t 8 -q 10 -Y -f igmelissa/sam_sai/aln_"."$id".".sai /project/evolgen/data/local/lycaeides/all_lycaeides_genomes/assembledLycaeidesGenomes/fasta/Lmelissa2FinalAssembly $fastq\n";
	system "bwa samse -n 1 -r '\@RG\tID:$id\tPL:ILLUMINA\tLB:$id\tSM:$id' -f igmelissa/sam_sai/mem_"."$id".".sam /project/evolgen/data/local/lycaeides/all_lycaeides_genomes/assembledLycaeidesGenomes/fasta/Lmelissa2FinalAssembly igmelissa/sam_sai/aln_"."$id".".sai $fastq\n";
	
	## Command using bwa mem instead
    #system "bwa mem -R '\@RG\\tID:$id\\tPL:ILLUMINA\\tLB:$id\\tSM:$id\\n\@CO' ../data/local/alfalfa_genome/CADL_HM342.v0.95P.fasta $fastq > sam_sai/mem_"."$id".".sam\n";

}
