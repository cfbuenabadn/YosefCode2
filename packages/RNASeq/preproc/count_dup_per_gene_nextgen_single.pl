use strict;

#make sure util.pl will load from the same dir as the main script
use FindBin;
use lib $FindBin::Bin;
use rnaseqPreprocUtils;

my($inp,$out,$rsemGTFFile,$rsemDictionary,$only_uniquely_alligned)=@ARGV;
## "$inp" is bam filename
## "$out" is output text file

#$rsemGTFFile="/data/yosef/index_files/mm10_4brain/index/rsem_index/combinedGTF_4brain.gtf";
#$rsemDictionary="/data/yosef/index_files/mm10_4brain/index/rsem_index/rsemDictionary/mm10_4brain_rsemGeneMapping.txt";
#perl count_dup_per_gene_nextgen.pl /data/yosef/BRAIN/processed2/150202_HS2A/Project_Ngai/OEP01_N710_S501_CGAGGCTG-TAGATCGC_L001_R1_combined/rsem_output/picard_output/sorted.bam b_test_count_dup.txt -1 -1 1


my $prv_chrom="";
my %uread=();my %read=();my %uread_no_strand=();
my %gene_dup_reads=();my %gene_tot_reads=();
my $counter=0;my $npos_strand=0;my $rcounter=0;my $drcounter=0;

my %genome_dat=();my %genes=();
print STDERR  "Reading genomic annotation for chromosome $prv_chrom\n";
read_genes(\%genome_dat,\%genes,$rsemGTFFile,$prv_chrom);
print STDERR  "Read ".scalar(keys %genes)." loci\n";


open (FD_sam, "samtools view -h $inp |");
while(<FD_sam>){
        if(/^\@/){next;}chomp;my $l=$_;
	my ($nm,$bitmap,$chrom,$pos,$mapq,$cigar,$rnext,$pnext,$tlen,$Seq,$qual,@rest)=split(/\t/,$l);
	if(0 && $only_uniquely_alligned){ #Disabled since the flag is mimssing (for some reason) from Bowtie output
		my $my_score=0;if ($l=~/AS\:i\:(\d+)/){
			$my_score=$1;my $other_score=0;if ($l=~/XS\:i\:(\d+)/){if ($my_score<=$1) {next;}}
		}
	}
		
        if(!($chrom=~/\*/ || $bitmap&0x200 || $bitmap&0x4 || ($only_uniquely_alligned && ($bitmap&0x100)))){
                if($prv_chrom ne $chrom && $prv_chrom ne "" && $chrom!~/\*/){
			calc_dup($prv_chrom);%read=();%uread=();%uread_no_strand=();
		}$prv_chrom=$chrom;
		my $strand=1;if($bitmap&0x10){$strand=0;};$npos_strand+=$strand;$counter++;
		my $key="$chrom#$pos#$strand";
		if(!exists($read{"$nm#$key"})){
			$read{"$nm#$key"}=1;$uread{$key}++;
			$uread_no_strand{"$chrom#$pos"}{$strand}++;
		}  
        }
}close FD_sam;calc_dup($prv_chrom);
print STDERR  "\nDONE. DATA FOR ".(scalar(keys %gene_tot_reads)).",".(scalar(keys %gene_dup_reads))." genes\n";

my $posnegrat=$npos_strand/$counter;
my $dupr=-1;if($rcounter > 0){$dupr=$drcounter/$rcounter;}

print STDERR "\nFound $drcounter duplicated out of $rcounter reads; ratio=$dupr |  pos-neg ratio=$posnegrat\n";
open h_out,">$out";
#posnegrat - count how many reads are mapped to forward strand and how many to reverse
#dupr - duplicate read ratio, rcounter - duplicate read counter
print h_out "$posnegrat\t$dupr\t$drcounter\t$rcounter\n";
close h_out;

open h_out,">$out\.genes.txt";my $N=0;
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
        print h_out "$id\t$dup\t$tot\t$rat\n";
}close h_out;close FD;
print STDERR "Printed dup stats for $N genes\n";



sub calc_dup{
	my ($prv_chrom)=@_;
	
	my $drcounter1=0;foreach my $k(keys %uread){$drcounter1+=$uread{$k}-1;}
	print STDERR "In $prv_chrom :: Found $drcounter1 duplicates from ".scalar(keys %uread)." unique reads and ".scalar(keys %read)." total reads\n";
	$rcounter+=scalar(keys(%read));$drcounter+=$drcounter1;
	my $prev_n_g=scalar(keys %gene_tot_reads);
	
	
	foreach my $k(keys %uread_no_strand){
		my($chr,$pos)=split(/\#/,$k);
		my @hits;get_overlapping_genes(\@hits,\@{$genome_dat{$prv_chrom}},\%{$genes{$prv_chrom}},$pos);
		foreach my $g(@hits){
			my $x=$uread_no_strand{$k}{"0"}+$uread_no_strand{$k}{"1"};
			$gene_tot_reads{$g}+=$x;
			$gene_dup_reads{$g}+=$x-scalar(keys %{$uread_no_strand{$k}});
		}
	}
	print STDERR  "Done searching; obtained overlap data for ".(scalar(keys %gene_tot_reads)-$prev_n_g)." genes\n\n";
}

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


