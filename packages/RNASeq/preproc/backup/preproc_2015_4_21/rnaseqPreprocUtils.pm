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
    my ($genome_dat,$transcripts,$GTFFile,$CHR)=@_;
    open FD,$GTFFile;my $cntr=0;
    while(<FD>){
        if (/^#/) {next;}
	chomp;my ($chr,$x1,$x2,$from,$to,$x3,$strand,$x4,$gene)=split(/\t/);
        my $xref;if ($gene=~/gene_id\ \"([\w\-_]+)\"/) {$xref=$1;}else{print "Warning (1): cannot parse $gene [$_]\n";next;}
        if ($CHR ne ""){
            if($chr ne $CHR) {if($cntr>0){last;}next;}
            if($cntr++ % 1e+5==0){print STDERR ".";}
            push @{$genome_dat},$to;$transcripts->{$to}{$xref}=$from;
            #print STDERR "\npushing: $from\t$to\t$gene\n";
        }else{
            if($cntr++ % 1e+5==0){print STDERR ".";}
            push @{$genome_dat->{$chr}},$to;$transcripts->{$chr}{$to}{$xref}=$from;
        }
    }
    close FD;
    if ($CHR ne ""){
        @{$genome_dat}=sort { $a <=> $b}  @{$genome_dat};
    }else{
        foreach my $chr(keys %{$genome_dat}){
            @{$genome_dat->{$chr}}=sort { $a <=> $b}  @{$genome_dat->{$chr}};    
        }
    }
}


sub get_overlapping_genes{
    my($res,$genome_dat,$transcripts,$pos)=@_;
    my $l=bin_search(\@{$genome_dat},$pos,0,scalar(@{$genome_dat})-1);
    #print STDERR "Location of pos $pos is $l\n";
    
    my $l1=$l;my $found=0;my $N=20;my %arr=();
    while ($l1>$l-$N || $found){
        my $loc=@{$genome_dat}[$l1];$found=0;
        foreach my $g(keys %{$transcripts->{$loc}}){
            if ($transcripts->{$loc}{$g}<=$pos && $loc>=$pos){
                $found=1;$arr{$g}=1;
            }    
        }$l1--;if($l1<0){last;}
    }
    $l1=$l+1;$found=0;
    while ($l1<$l+$N || $found){
        my $loc=@{$genome_dat}[$l1];$found=0;
        foreach my $g(keys %{$transcripts->{$loc}}){
            if ($transcripts->{$loc}{$g}<=$pos && $loc>=$pos){
                $found=1;$arr{$g}=1;
            }    
        }$l1++;if($l1>=scalar(@{$genome_dat})){last;}
    }
    foreach my $g(keys %arr){push @{$res},$g;}
}
1;

