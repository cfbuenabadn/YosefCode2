use strict;

my($OUTPUT_FOLDER, $refFlatFile, $ribosomalIntervalsFile, $f, $f2)=@ARGV;
my $FASTQC_CMD = "/opt/genomics/bin/fastqc";

print "Starting qc.pl....\n";
my $isPairedEnd = ($f2 ne "");


#for FastQC help type: "fastqc --help"

if(1){
    ####COMPUTE####
    system("$FASTQC_CMD  $f -o  $OUTPUT_FOLDER/fastqc_output\n");
    system("perl check_contam.pl $f $OUTPUT_FOLDER/fastqc_output/primer.1.txt\n");
    
    if($isPairedEnd)
    {
	system("$FASTQC_CMD  $f2 -o  $OUTPUT_FOLDER/fastqc_output\n");
	system("perl check_contam.pl $f2 $OUTPUT_FOLDER/fastqc_output/primer.2.txt\n");
    }
    
    system("/opt/genomics/bin/CollectRnaSeqMetrics TMP_DIR=$OUTPUT_FOLDER/temp INPUT=$OUTPUT_FOLDER/picard_output/sorted.bam OUTPUT=$OUTPUT_FOLDER/picard_output/rna_metrics.txt CHART=$OUTPUT_FOLDER/picard_output/rna_coverage.pdf REF_FLAT=$refFlatFile STRAND=NONE RIBOSOMAL_INTERVALS=$ribosomalIntervalsFile\n");
    system("/opt/genomics/bin/CollectAlignmentSummaryMetrics INPUT=$OUTPUT_FOLDER/picard_output/sorted.bam OUTPUT=$OUTPUT_FOLDER/picard_output/aln_metrics.txt\n");
    system("/opt/genomics/bin/CollectInsertSizeMetrics INPUT=$OUTPUT_FOLDER/picard_output/sorted.bam OUTPUT=$OUTPUT_FOLDER/picard_output/aln_metrics1.txt H=$OUTPUT_FOLDER/picard_output/ins_metrics.histogram.pdf\n");
    system("perl count_dup.pl $OUTPUT_FOLDER/picard_output/sorted.bam $OUTPUT_FOLDER/picard_output/dup.txt\n");
}

if(1){
    ####AGGREGATE####
    
    #  Total number of reads
    #  % of aligned reads
    my $nonExistingFiles = 0;
    my ($s1, @a, $nreads, $nalign) = 0;
    if(-e "$OUTPUT_FOLDER/tophat_output/logs/bowtie.left_kept_reads.fixmap.log"){ 
	open FD,"$OUTPUT_FOLDER/tophat_output/logs/bowtie.left_kept_reads.fixmap.log";
	$s1=<FD>;chomp $s1;@a=split(/\s/,$s1);$nreads=$a[3];
	$s1=<FD>;chomp $s1;@a=split(/\s/,$s1);$nalign=$a[8];
	close FD;
    }else{
	 $nonExistingFiles++;   
    }
     if(-e "$OUTPUT_FOLDER/tophat_output/logs/bowtie.right_kept_reads.fixmap.log"){
	open FD,"$OUTPUT_FOLDER/tophat_output/logs/bowtie.right_kept_reads.fixmap.log";
	$s1=<FD>;chomp $s1;@a=split(/\s/,$s1);$nreads+=$a[3];
	$s1=<FD>;chomp $s1;@a=split(/\s/,$s1);$nalign+=$a[8];
	close FD;my $r_align=100*$nalign/$nreads;
    }else{
	$nonExistingFiles++; 
    }
    
    if($isPairedEnd && $nonExistingFiles > 0){
	    die "paired-end sample but less than 2 alignment files found";
    }elsif(!($isPairedEnd) && $nonExistingFiles > 1){
	    die "paired-end sample but less than 1 alignment files found";
    }
    
    #  Sequence duplication level (Total)
    my $total_dup=0;my $cnt=0;
    foreach my $f(<$OUTPUT_FOLDER/Fastqc/*_fastqc/fastqc_data.txt>){
	open FD,$f;while(<FD>){if(/#Total Duplicate Percentage/){chomp;my @a=split(/\t/);$total_dup+=$a[1];last;}}close FD;
        $cnt++;
    }if($cnt!=2){print "Warning: problem with Fastqc folders\n";}
    else{$total_dup/=2;}
    
#     #  Abundance of primer sequences
    open FD,"$OUTPUT_FOLDER/Fastqc/primer.1.txt";
    @a=split(/\t/,<FD>);my $nhits=$a[1];my $nseq=$a[2];close FD;
    open FD,"$OUTPUT_FOLDER/Fastqc/primer.2.txt";
    @a=split(/\t/,<FD>);$nhits+=$a[1];$nseq+=$a[2];close FD;   
    my $primer=$nhits/$nseq;
   
    
    #  MEDIAN_INSERT_SIZE
    #  MEDIAN_ABSOLUTE_DEVIATION
    open FD,"$OUTPUT_FOLDER/Fastqc/aln_metrics1.txt";
    my $insert_sz_avg=0;my $insert_sz_std=0;
    while(<FD>){
      if(/## METRICS CLASS/){
           <FD>;my @a=split(/\t/,<FD>);
           $insert_sz_avg=$a[0];$insert_sz_std=$a[1];
           last;
        }
    }close FD;
    
    
     #  % of unique alignments (complexity) 
    open FD,"$OUTPUT_FOLDER/Fastqc/dup.txt";
    @a=split(/\t/,<FD>);my $complexity=1-$a[1];my $left_right=$a[0];
    
    
    # #  PCT_RIBOSOMAL_BASES
    # #  PCT_CODING_BASES
    # #  PCT_UTR_BASES
    # #  INTRONIC_BASES
    # #  PCT_INTERGENIC_BASES
    # #  PCT_MRNA_BASES
    # #  MEDIAN_CV_COVERAGE (evenness)
    # #  MEDIAN_5PRIME_BIAS
    # #  MEDIAN_3PRIME_BIAS
    # #  MEDIAN_5PRIME_TO_3PRIME_BIAS
    open FD,"$OUTPUT_FOLDER/Fastqc/rna_metrics.txt";
    my @metric=();
    while(<FD>){
        if(/## METRICS CLASS/){
            <FD>;my @a=split(/\t/,<FD>);
            foreach my $i(10..15,18..21){push @metric,$a[$i];}
            last;
        }
    }close FD;

    print "Output to: $OUTPUT_FOLDER/summary.txt\n";
    open h_out, ">$OUTPUT_FOLDER/summary.txt";
    print h_out "NREADS\t$nreads\nNALIGNED\t$nalign\nRALIGN\t$r_align\n";
    print h_out "TOTAL_DUP\t$total_dup\nPRIMER\t$primer\nINSERT_SZ\t$insert_sz_avg\nINSERT_SZ_STD\t$insert_sz_std\nCOMPLEXITY\t$complexity\n";
    #INTRONIC_BASES is actually the field PCT_INTRONIC_BASES
    my @arr=("PCT_RIBOSOMAL_BASES","PCT_CODING_BASES","PCT_UTR_BASES","INTRONIC_BASES","PCT_INTERGENIC_BASES","PCT_MRNA_BASES","MEDIAN_CV_COVERAGE","MEDIAN_5PRIME_BIAS","MEDIAN_3PRIME_BIAS","MEDIAN_5PRIME_TO_3PRIME_BIAS");
    for(my $i=0;$i<scalar(@arr);$i++){print h_out "$arr[$i]\t$metric[$i]\n";}
    close h_out;
}

 