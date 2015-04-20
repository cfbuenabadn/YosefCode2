use strict;
my($inp,$out,$rsemGTFFile,$rsemDictionary,$only_uniquely_alligned)=@ARGV;
## "$inp" is bam filename
## "$out" is output text file

#$rsemGTFFile="/data/yosef/index_files/mm10_4brain/index/rsem_index/combinedGTF_4brain.gtf";
#$rsemDictionary="/data/yosef/index_files/mm10_4brain/index/rsem_index/rsemDictionary/mm10_4brain_rsemGeneMapping.txt";
#perl count_dup_per_gene.pl /data/yosef/BRAIN/processed2/150202_HS2A/Project_Ngai_AsSingle/OEP01_N710_S501_CGAGGCTG-TAGATCGC_L001_R1_combined/rsem_output/picard_output/sorted.bam b_test_count_dup.txt -1 -1 1

#description of rsem output: http://deweylab.biostat.wisc.edu/rsem/rsem-calculate-expression.html

my $GRID=10;
my $LOAD_GENES_IN_ADVANCE=1;my $N=0;
my %genes=();my $N1=0;my %gcntr=();

open FD,$rsemGTFFile;
#NC_000067.6     Gnomon  exon    3670552 3671742 .       -       .       transcript_id "rna0"; gene_id "gene0"; gene_name "Xkr4";

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
my %ufrag=();my %frag=();my %ufrag_no_strand=();
my %uread=();my %read=();my %uread_no_strand=();
my %gene_dup_reads=();my %gene_tot_reads=();
my %gene_dup_frags=();my %gene_tot_frags=();

my $counter=0;my $dfcounter=0;my $fcounter=0;my $npos_strand=0;my $rcounter=0;my $drcounter=0;

open (FD_sam, "samtools view -h $inp |");

while(<FD_sam>){
        if(/^\@/){next;}chomp;my $l=$_;
	my ($nm,$bitmap,$chrom,$pos,$mapq,$cigar,$rnext,$pnext,$tlen,$Seq,$qual,@rest)=split(/\t/,$l);
	if(0 && $only_uniquely_alligned){ #Disabled since the flag is mimssing (for some reason) from Bowtie output
		my $my_score=0;if ($l=~/AS\:i\:(\d+)/){
			$my_score=$1;my $other_score=0;
			if ($l=~/XS\:i\:(\d+)/){if ($my_score<=$1) {next;}}
		}
	}
		
        if(!($chrom=~/\*/ || $bitmap&0x200 || $bitmap&0x4 || $bitmap&0x8 || ($only_uniquely_alligned && ($bitmap&0x100)))){
                if($prv_chrom ne $chrom && $prv_chrom ne "" && $chrom!~/\*/){
		
			my $dfcounter1=0;foreach my $k(keys %ufrag){$dfcounter1+=$ufrag{$k}-1;}
			print STDERR "In $prv_chrom :: Found $dfcounter1 duplicates from ".scalar(keys %ufrag)." unique fragments and ".scalar(keys %frag)." total fragments\n";
			$fcounter+=scalar(keys(%frag));$dfcounter+=$dfcounter1;
			
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
							foreach my $k1(keys  %{$ufrag_no_strand{$k}}){$x+=$ufrag_no_strand{$k}{$k1};}
							$gene_tot_frags{$g}+=$x;
							$gene_dup_frags{$g}+=$x-scalar(keys %{$ufrag_no_strand{$k}});
						}
					}
				}
			}
			print STDERR  "\nDone searching; obtained overlap data for ".(scalar(keys %gene_tot_reads)-$prev_n_g)." genes\n";
			%frag=();%ufrag=();%frag=();%ufrag=();%uread=(); %read=();
		}$prv_chrom=$chrom;
                
                my $strand=1;if($bitmap&0x10){$strand=0;};$npos_strand+=$strand;$counter++;
                if($bitmap&0x40){
                        my $mate_strand=1;if($bitmap&0x20){$mate_strand=0;}
                        my $key="$chrom#$pos#$strand#$pnext#$mate_strand";
                        #my $key="$chrom-$pos-$pnext";
                        if(!exists($frag{"$nm#$key"})){
                        	$frag{"$nm#$key"}=1;$ufrag{$key}++;
				$ufrag_no_strand{"$chrom#$pos"}{$key}++;
			}
			
                        $key="$chrom#$pos#$strand";
                        if(!exists($read{"$nm#$key"})){$read{"$nm#$key"}=1;$uread{$key}++;$uread_no_strand{"$chrom#$pos"}{$strand}++;}
                        
                }
        }
}
close FD_sam;%genes=();


my $posnegrat=$npos_strand/$counter;
my $dupf=-1;if($fcounter > 0){$dupf=$dfcounter/$fcounter;}
my $dupr=-1;if($rcounter > 0){$dupr=$drcounter/$rcounter;}

print STDERR "\nFound $dfcounter duplicated out of $fcounter fragments; ratio=$dupf |  pos-neg ratio=$posnegrat\n";
print STDERR "\nFound $drcounter duplicated out of $rcounter reads; ratio=$dupr\n";
open h_out,">$out";
#posnegrat - count how many reads are mapped to forward strand and how many to reverse
#dupf - duplcate fragment ratio, dfcounter - duplicate fragment counter
#dupr - duplicate read ratio, rcounter - duplicate read counter
print h_out "$posnegrat\t$dupf\t$dfcounter\t$fcounter\t$dupr\t$drcounter\t$rcounter\n";
close h_out;

open h_out,">$out\.genes.txt";$N=0;
open FD,$rsemDictionary;
while(<FD>){
	if (/^#/) {next;}
	chomp;my ($id,@rest)=split(/\t/);
        my $tot=0;my $dup=0;my $rat="NaN";
        my $totf=0;my $dupf=0;my $ratf="NaN";
        if (exists($gene_dup_reads{$id})) {
	      $tot=$gene_tot_reads{$id};
	      $dup=$gene_dup_reads{$id};
	      $rat=((1+$gene_dup_reads{$id})/(1+$gene_tot_reads{$id}));
	      $totf=$gene_tot_frags{$id};
	      $dupf=$gene_dup_frags{$id};
	      $ratf=((1+$gene_dup_frags{$id})/(1+$gene_tot_frags{$id}));
	      $N++;
        }
        print h_out "$id\t$dup\t$tot\t$rat\t$dupf\t$totf\t$ratf\n";
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
