use strict;
my ($F,$outf)=@ARGV;


#miSeq: 25bp; HiSeq: 125
my %SARR=();my $L=12;
#"AAGCAGTGGTATCAACGCAGAGTAC"
foreach my $S("AAGCAGTGGTATCAACGCAGAGT"){
    $SARR{substr($S,0,$L)}{"kind"}="F";$SARR{substr($S,0,$L)}{"STR"}=$S;
    $SARR{substr($S,length($S)-$L,$L)}{"kind"}="B";$SARR{substr($S,length($S)-$L,$L)}{"STR"}=$S;
}
print "Analyzing $F\n";
if($F=~/\.gz$/){system("gunzip $F\n");$F=~s/fastq\.gz/fastq/;}
my $found=0;my $sum=0;my $ne=0;
open FD,$F;<FD>;
while(<FD>){
    $found=0;chomp;my $seq_read=$_;
    foreach my $S(keys %SARR){
        while($seq_read=~/$S/g){
            my $seq_read1=$seq_read;my $str=$SARR{$S}{"STR"};my $pos=$-[0];
            if($SARR{$S}{"kind"} eq "B"){$seq_read1=reverse($seq_read1);$str=reverse($str);$pos=length($seq_read1)-$+[0];}
            my $l=length($seq_read1)-$pos;if($l>length($str)){$l=length($str);}
            if(substr($seq_read1,$pos,$l) eq substr($str,0,$l)){$found=1;last;}
        }
        if($found){last;}
    }
    $sum+=$found;
    if($ne++ % 20000==1){print ".";}
    foreach my $x(1..3){<FD>;}
}

print "\n$F: $sum matches out of $ne  [".( (1+$sum)/(1+$ne))."]\n";
open h_out,">$outf";
print h_out "$F\t$sum\t$ne\t".( (1+$sum)/(1+$ne))."\n";
