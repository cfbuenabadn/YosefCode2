use strict;
my ($WORK_FOLDER,$dict_file,$include_raw_read)=@ARGV;

#INPUTS:
#Working directory (upper level, all single cells are subdirectories)
#Gene dictionaty file

#Option for dict file:
#Humans: "/data/yosef/index_files/hg38/index/rsem_dict.txt";
#Mice: "/data/yosef/index_files/mm10_4brain/index/rsem_index/rsemDictionary/mm10_4brain_rsemGeneMapping.txt";

#OUTPUTS:
#cell_list.txt
#dup_reads_per_gene_table.txt
#gene_list.txt
#qc_table.txt
#cuff_fpkmTable.txt
#cuff_readCountsTable.txt


#FOR DEBUG:
#$dict_file="/data/yosef/index_files/mm10_4brain/index/rsem_index/rsemDictionary/mm10_4brain_rsemGeneMapping.txt";
#$include_raw_read=0;
#perl collect_dat_cufflinks.pl /home/eecs/allonwag/data/BRAIN/processed3/150202_HS2A/Project_Ngai/
#with include raw read:
#perl /home/eecs/allonwag/project/singleCell/allon_script/YosefCode/RNASeq/preproc/collect_dat_cufflinks.pl /home/eecs/allonwag/data/BRAIN/processed3/150202_HS2A/Project_Ngai/ /data/yosef/index_files/mm10_4brain/index/rsem_index/rsemDictionary/mm10_4brain_rsemGeneMapping.txt 1

my $cntr=0;my %fpkm=();
my %read_count=();my %read_dup=();my %read_count_unique=();
my @qc_info=();my @cell_nm=();my @failed_cells=();

print "Processing $WORK_FOLDER\n";

foreach my $f(<$WORK_FOLDER/*combined/tophat_output/cuff_output/genes.fpkm_tracking>){
    my $dir=$f;$dir=~s/genes\.fpkm_tracking//;
    if (!(-e "$dir/../summary.txt")) {print "Warning: Could not find $dir/../summary.txt\n";push @failed_cells,$f;next;}
    if ($include_raw_read && !(-e "$dir/../picard_output/dup.txt.genes.txt")) {print "Warning: Could not find $dir/../picard_output/dup.txt.genes.txt\n";push @failed_cells,$f;next;}
    if ($include_raw_read && !(-e "$dir/../picard_output/dup_unique.txt.genes.txt")) {print "Warning: Could not find $dir/../picard_output/dup_unique.txt.genes.txt\n";push @failed_cells,$f;next;}
    $cntr++;
    #if ($cntr>10){last;} #DDEEEEBBBUUUGG
    
}my @zeros=();for(my $i=0;$i<$cntr;$i++){push @zeros,0;}$cntr=0;


foreach my $f(<$WORK_FOLDER/*combined/tophat_output/cuff_output/genes.fpkm_tracking>){
    my $dir=$f;$dir=~s/genes\.fpkm_tracking//;
    if (!(-e "$dir/../summary.txt")) {print "Warning: Could not find $dir/../summary.txt\n";next;}
    if ($include_raw_read && !(-e "$dir/../picard_output/dup.txt.genes.txt")) {print "Warning: Could not find $dir/../picard_output/dup.txt.genes.txt\n";next;}
    if ($include_raw_read && !(-e "$dir/../picard_output/dup_unique.txt.genes.txt")) {print "Warning: Could not find $dir/../picard_output/dup_unique.txt.genes.txt\n";next;}
    
    my $nm=$dir;if ($nm=~/$WORK_FOLDER\/([\w\_\-]+)\/tophat_output/) {$nm=$1;}else{die "Could not parse folder name $dir\n";}
    $cell_nm[$cntr]=$nm;
    
    print STDERR ".";open FD,$f;my %v=();
	while(<FD>){
	    chomp;my @a=split;
	    my $id=$a[4];
	    if(!exists($v{$id})  || $v{$id}<$a[4]){
		$v{$id}=$a[9];
	    }
	}close FD;
        foreach my $id(keys %v){
	    if(!exists($fpkm{$id})){@{$fpkm{$id}}=@zeros;}
            @{$fpkm{$id}}[$cntr]=$v{$id};
        }
        {
            open FD,"$dir/../summary.txt";my @b=();
            while(<FD>){chomp;my @a=split;push @{$qc_info[$cntr]},$a[1];}
        }
	
	if ($include_raw_read){
	    open FD,"$dir/../picard_output/dup.txt.genes.txt";
	    while(<FD>){
		my($id,$dup_read,$tot_read,$rat_read,$dup_frag,$tot_frag,$rat_frag)=split(/\t/);
		if(!exists($read_count{$id})){@{$read_count{$id}}=@zeros;}
		@{$read_count{$id}}[$cntr]=$tot_read;
		if(!exists($read_dup{$id})){@{$read_dup{$id}}=@zeros;}
		@{$read_dup{$id}}[$cntr]=$rat_read;
	    }close FD;
	    
	    open FD,"$dir/../picard_output/dup_unique.txt.genes.txt";
	    while(<FD>){
		my($id,$dup_read,$tot_read,$rat_read,$dup_frag,$tot_frag,$rat_frag)=split(/\t/);
		if(!exists($read_count_unique{$id})){@{$read_count_unique{$id}}=@zeros;}
		@{$read_count_unique{$id}}[$cntr]=$tot_read;
	    }close FD;
	    
	}
        $cntr++;
	
	#if ($cntr>10){last;}#DDEEEEBBBUUUGG
	
}
print STDERR "Read $cntr sampels\n";
if(!(-e "$WORK_FOLDER/cuff/")){system("mkdir $WORK_FOLDER/cuff/");}




open h_out,">$WORK_FOLDER/cuff/cell_list.txt";
print h_out join "\n",@cell_nm;close h_out;
open h_out_qc,">$WORK_FOLDER/cuff/qc_table.txt";
for(my $i=0;$i<scalar(@qc_info);$i++){
    if (scalar(@{$qc_info[0]})!=scalar(@{$qc_info[$i]})) {die "3: Integrity\n";}
}if(scalar(@qc_info)!=scalar(@cell_nm)){die "4: Integrity\n";}
for(my $j=0;$j<scalar(@{$qc_info[0]});$j++){
    my $delim="";
    for(my $i=0;$i<scalar(@cell_nm);$i++){
	print h_out_qc "$delim$qc_info[$i][$j]";$delim="\t";
    }print h_out_qc "\n";
}close h_out_qc;


open h_out_dat,">$WORK_FOLDER/cuff/cuff_fpkmTable.txt";
if ($include_raw_read){
    open h_out_dat_read,">$WORK_FOLDER/cuff/cuff_readCountsTable.txt";
    open h_out_dat_dup,">$WORK_FOLDER/cuff/dup_reads_per_gene_table.txt";
    open h_out_dat_read_unique,">$WORK_FOLDER/cuff/dup_reads_per_gene_onlyUnique_table.txt";
}
open h_out_gene,">$WORK_FOLDER/cuff/gene_list.txt";
print h_out_gene "#This file maps mm10's gene symbols to gene categories\n";
print h_out_gene "#Nir, Match 2015\n";
print h_out_gene "#Gene Sybmol\tGene Sybmol\tType\n";

open FD,$dict_file; my %printed=();
while(<FD>){
    if (/^#/){next;}
    chomp;my ($tmp,$p,$type)=split(/\t/);
    print  h_out_gene "$p\t$p\t$type\n";
    if (exists($fpkm{$p})) {
	print  h_out_dat join "\t",@{$fpkm{$p}};
	if ($include_raw_read){
	    if (!exists($read_count{$p})) {die "Integrity: no read count data for gene $p\n";}	
	    print  h_out_dat_read join "\t",@{$read_count{$p}};
	    if (!exists($read_count_unique{$p})) {die "Integrity: no unique read count data for gene $p\n";}	
	    print  h_out_dat_read_unique join "\t",@{$read_count_unique{$p}};
	    if (!exists($read_dup{$p})) {die "Integrity: no read dup data for gene $p\n";}	
	    print  h_out_dat_dup join "\t",@{$read_dup{$p}};
	}
    }else{
	print  h_out_dat join "\t",@zeros;
	if ($include_raw_read){
	    print  h_out_dat_read join "\t",@zeros;
	    print  h_out_dat_read_unique join "\t",@zeros;
	    print  h_out_dat_dup join "\t",@zeros;
	}
    }
    print h_out_dat "\n";print h_out_dat_read "\n";print h_out_dat_dup "\n";
    $printed{$p}=1;
}close FD;

my $fcntr=0;
foreach my $p(keys %fpkm){
    if (exists($printed{$p})){next;}
    print  h_out_gene "$p\t$p\t-1\n";
    print  h_out_dat join "\t",@{$fpkm{$p}};
    if ($include_raw_read){
	if (!exists($read_count{$p})) {die "Integrity: no read count data for gene $p\n";}	
	print  h_out_dat_read join "\t",@{$read_count{$p}};print h_out_dat_read "\n";
	if (!exists($read_count_unique{$p})) {die "Integrity: no unique read count data for gene $p\n";}	
	print  h_out_dat_read_unique join "\t",@{$read_count_unique{$p}};print h_out_dat_read_unique "\n";
	if (!exists($read_dup{$p})) {die "Integrity: no read dup data for gene $p\n";}
	print  h_out_dat_dup join "\t",@{$read_dup{$p}};print h_out_dat_dup "\n";
    }
    print h_out_dat "\n";
    $fcntr++;
}
print "Printed ".scalar(keys %printed)." genes, plus $fcntr additional genes not in dictionary\n";
close h_out_gene;close h_out_dat;

if ($include_raw_read){
    close h_out_dat_read;close h_out_dat_dup;close h_out_dat_read_unique;
}

    
print "\n\nDONE!! output to $WORK_FOLDER/cuff\n";



print "\ndone collecting all data!";
if(scalar(@failed_cells)==0){print "no cells had errors while collecting their gene expression";}
else{
    print "\nThe following cells had errors while collecting their gene expression:\n";
    print join "\n",@failed_cells;
}


