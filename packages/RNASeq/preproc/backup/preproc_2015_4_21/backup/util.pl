use strict;

sub bin_search{
    my($arr,$val,$frm,$to)=@_;
    my $mid=int(($frm+$to)/2);
    #print "$frm\t$to\t$mid\t@{$arr}[$mid]\t$val\n";
    if($frm>=$mid || @{$arr}[$mid]==$val){
        #die "Found val=$val in position $mid of value @{$arr}[$mid]\n";
        my $min=abs(@{$arr}[$mid]-$val);my $best=$mid;
        for(my $i=-1;$i<=1;$i++){if(abs(@{$arr}[$mid+$i]-$val)<$min){$min=abs(@{$arr}[$mid+$i]-$val);$best=$mid+$i}}
        return $best;
    }
    if(@{$arr}[$mid]>$val){return bin_search($arr,$val,$frm,$mid);}
    else{return bin_search($arr,$val,$mid,$to);}
}

sub read_genes{
    my ($genome_dat,$transcripts,$GTFFile,$side,$CHR)=@_;
    open FD,$GTFFile;my $cntr=0;
    while(<FD>){
        if (/^#/) {next;}
	chomp;my ($chr,$x1,$x2,$from,$to,$x3,$strand,$x4,$gene)=split(/\t/);
        if ($chr ne $CHR) {if($cntr>0){last;}next;}
        if($cntr++ % 1e+5==0){print STDERR ".";}
	my ($xref,$gname);if ($gene=~/gene_id\ \"([\w\-_]+)\"/) {$xref=$1;}else{print "Warning (1): cannot parse $gene [$_]\n";next;}
        my $strnd_str="+";if ($side==0) {$strnd_str="-";}
        if($strand eq $strnd_str){push @{$genome_dat},$from;push @{$transcripts->{$from}},$gene;}
        else{push @{$genome_dat},$to;push @{$transcripts->{$to}},$gene;}
    }
    close FD;
    @{$genome_dat}=sort { $a <=> $b}  @{$genome_dat};
}


sub get_overlapping_genes{
    my($res,$genome_dat_5p,$transcripts_5p,$genome_dat_3p,$transcripts_3p,$pos)=@_;
    my $l=abs(bin_search(@{$genome_dat_3p},$pos,0,scalar(@{$genome_dat_3p})-1)-$pos);
    my $l2=abs(bin_search(@{$genome_dat_5p},$pos,0,scalar(@{$genome_dat_5p})-1)-$pos);
				
    @{$res}=();
}
1;

#perl run_centipede.pl SIX5_disc1  ../../../dat/ChIP/H7AG5ADXX/data/CGTACTAG-T16_T36-131009_s5/ temp                                                                                                                                                   