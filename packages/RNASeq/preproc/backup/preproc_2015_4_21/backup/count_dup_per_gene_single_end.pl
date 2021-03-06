use strict;
my($inp,$out,$rsemGTFFile,$rsemDictionary)=@ARGV;
## "$inp" is bam filename
## "$out" is output text file


$rsemGTFFile="/data/yosef/index_files/mm10_4brain/index/rsem_index/combinedGTF_4brain.gtf";
$rsemDictionary="/data/yosef/index_files/mm10_4brain/index/rsem_index/rsemDictionary/mm10_4brain_rsemGeneMapping.txt";
#perl count_dup_per_gene_single_end.pl /data/yosef/BRAIN/processed2/150202_HS2A/Project_Ngai_AsSingle/OEP01_N710_S501_CGAGGCTG-TAGATCGC_L001_R1_combined/rsem_output/picard_output/sorted.bam b_test_count_dup.txt

my $GRID=10;
my $LOAD_GENES_IN_ADVANCE=1;my $N=0;
my %genes=();my $N1=0;my %gcntr=();

open FD,$rsemGTFFile;
while (<FD>){
	if (/^#/) {next;}
	chomp;my ($chr,$x1,$x2,$frm,$to,$x3,$strnd,$x4,$id)=split(/\t/);
	#if($chr ne "NC_000067.6"){if(scalar(keys %genes)>0){last;}next;}
	if($N1++ % 1e+5==0){print STDERR ".";}
	my ($xref,$gname);
	if ($id=~/gene_id\ \"([\w\-_]+)\"/) {$xref=$1;}else{print "Warning (1): cannot parse $id [$_]\n";next;}
	#if ($id=~/gene_name\ \"([\w\-_]+)\"/) {$gname=$1;}else{print "Warning (2): annot parse $id\n";next;}
	if ($strnd=~/\+/){$strnd=1;}elsif ($strnd=~/\-/){$strnd=0;}else{die "Could not parse strand: $strnd\n";}
	for(my $p=$frm;$p<=$to;$p+=$GRID){
		push @{$genes{"$chr#$p"}},$xref;
		
	}$gcntr{$xref}=1;
	
}
print STDERR  "\nRead ".scalar(keys %genes). " loci and ".scalar(keys %gcntr)." genes\n";


my $prv_chrom="";
my %uread=();my %read=();my %uread_no_strand=();
my $counter=0;my $npos_strand=0;my $rcounter=0;my $drcounter=0;
my %gene_dup_reads=();my %gene_tot_reads=();
open (FD_sam, "samtools view -h $inp |");


while(<FD_sam>){
        if(/^\@/){next;}
	chomp;my ($nm,$bitmap,$chrom,$pos,$mapq,$cigar,$rnext,$pnext,@rest)=split(/\t/);
        if(!($chrom=~/\*/ || $bitmap&0x200 || $bitmap&0x4 )){
                if($prv_chrom ne $chrom && $prv_chrom ne "" && $chrom!~/\*/){
			my $drcounter1=0;foreach my $k(keys %uread){$drcounter1+=$uread{$k}-1;}
			print STDERR "In $prv_chrom :: Found $drcounter1 duplicates from ".scalar(keys %uread)." unique reads and ".scalar(keys %read)." total reads\n";
			$rcounter+=scalar(keys(%read));$drcounter+=$drcounter1;
		       
			my $prev_n_g=scalar(keys %gene_tot_reads);
			#my $printed=0;
			foreach my $k(keys %uread_no_strand){
				my($chr,$pos)=split(/\#/,$k);
				#my $key="$chr#$pos#$strnd";print "==$key==\n";if($printed++>500){die "\n";}
				for(my $p=$pos-int($GRID/2);$p<=$pos+int($GRID/2);$p++){
					my $key="$chr#$p";
					if (exists($genes{$key})) {
						foreach my $g(@{$genes{$key}}){
							my $x=$uread_no_strand{$k}{"0"}+$uread_no_strand{$k}{"1"};
							$gene_tot_reads{$g}+=$x;
							$gene_dup_reads{$g}+=$x-scalar(keys %{$uread_no_strand{$k}});
						}
					}
				}
			}
			print STDERR  "\nDone searching; obtained overlap data for ".(scalar(keys %gene_tot_reads)-$prev_n_g)." genes\n";
			%uread=(); %read=();
			#if($prv_chrom eq "NC_000067.6"){last;}
		}$prv_chrom=$chrom;
                
                my $strand=1;if($bitmap&0x10){$strand=0;};$npos_strand+=$strand;$counter++;
                my $key="$chrom#$pos#$strand";
                if(!exists($read{"$nm#$key"})){$read{"$nm#$key"}=1;$uread{$key}++;$uread_no_strand{"$chrom#$pos"}{$strand}++;}
        }
}
close FD_sam;
%uread=(); %read=();


my $posnegrat=$npos_strand/$counter;
my $dupr=-1;if($rcounter > 0){$dupr=$drcounter/$rcounter;}

print STDERR "\nFound $drcounter duplicated out of $rcounter reads; ratio=$dupr |  pos-neg ratio=$posnegrat\n";
open h_out,">$out";
#posnegrat - count how many reads are mapped to forward strand and how many to reverse
#dupf - duplcate fragment ratio, dfcounter - duplicate fragment counter
#dupr - duplicate read ratio, rcounter - duplicate read counter
print h_out "$posnegrat\tNaN\tNaN\tNaN\t$dupr\t$drcounter\t$rcounter\n";
close h_out;

open h_out,">$out\.genes.txt";$N=0;
open FD,$rsemDictionary;
while(<FD>){
	if (/^#/) {next;}
	chomp;my ($id,@rest)=split(/\t/);
        my $tot=0;my $dup=0;my $rat="NaN";
        if (exists($gene_dup_reads{$id})) {
	      $tot=$gene_tot_reads{$id};
	      $dup=$gene_dup_reads{$id};
	      $rat=((1+$gene_dup_reads{$id})/(1+$gene_tot_reads{$id}));
	      $N++;
        }
        print h_out "$id\t$dup\t$tot\t$rat\tNaN\tNaN\tNaN\n";
}close h_out;close FD;
print STDERR "Printed dup stats for $N genes\n";

#0x1 template having multiple segments in sequencing
#0x2 each segment properly aligned according to the aligner
#0x4 segment unmapped
#0x8 next segment in the template unmapped
#0x10 SEQ being reverse complemented
#0x20 SEQ of the next segment in the template being reversed
#0x40 the first segment in the template
#0x80 the last segment in the template
#0x100 secondary alignment
#0x200 not passing quality controls
#0x400 PCR or optical duplicate
#0x800 supplementary alignment



#Enumerate chrom order in gff:
#my %chr_arr=();my @chr_arr1=();
#open FD,"/data/yosef/index_files/mm10_4brain/index/GRCm38.p3.gff";
#while (<FD>){
#	if (/^#/) {next;}
#	chomp;my ($chr,@rest)=split(/\t/);
#	if (!exists($chr_arr{$chr})) {push @chr_arr1,$chr;$chr_arr{$chr}=1;print STDERR ".";}
#}
#print "\n";
#print join "\n",@chr_arr1;
#die "\n";
