use strict;
my($inp,$out)=@ARGV;
## "$inp" is bam filename
## "$out" is output text file

my $prv_chrom="";
my %frag=();my %ufrag=();my %uread=();my %read=();
my $counter=0;my $dfcounter=0;my $fcounter=0;my $npos_strand=0;my $rcounter=0;my $drcounter=0;
open (FD_sam, "samtools view -h $inp |");

while(<FD_sam>){
        if(/^\@/){next;}
        chomp;my ($nm,$bitmap,$chrom,$pos,$mapq,$cigar,$rnext,$pnext,@rest)=split(/\t/);
        if(!($chrom=~/\*/ || $bitmap&0x200 || $bitmap&0x4 || $bitmap&0x8)){
                my $strand=1;if($bitmap&0x10){$strand=0;};$npos_strand+=$strand;$counter++;
                if($bitmap&0x40){
                        my $mate_strand=1;if($bitmap&0x20){$mate_strand=0;}
                        my $key="$chrom-$pos-$strand-$pnext-$mate_strand";
                        #my $key="$chrom-$pos-$pnext";
                        if(!exists($frag{"$nm-$key"})){$frag{"$nm-$key"}=1;$ufrag{$key}++;}

			$key="$chrom-$pos-$strand";
                        if(!exists($read{"$nm-$key"})){$read{"$nm-$key"}=1;$uread{$key}++;}
                }
        }
        if($prv_chrom ne $chrom && $prv_chrom ne "" && $chrom!~/\*/){
                my $dfcounter1=0;foreach my $k(keys %ufrag){$dfcounter1+=$ufrag{$k}-1;}
                print "In $prv_chrom :: Found $dfcounter1 duplicates from ".scalar(keys %ufrag)." unique fragments and ".scalar(keys %frag)." total fragments\n";
                $fcounter+=scalar(keys(%frag));$dfcounter+=$dfcounter1;
                %frag=();%ufrag=();

		my $drcounter1=0;foreach my $k(keys %uread){$drcounter1+=$uread{$k}-1;}		                
		print "In $prv_chrom :: Found $drcounter1 duplicates from ".scalar(keys %uread)." unique reads and ".scalar(keys %read)." total reads\n";
                $rcounter+=scalar(keys(%read));$drcounter+=$drcounter1;
                %frag=();%ufrag=();%uread=(); %read=();

        }$prv_chrom=$chrom;
}
close FD_sam;

my $posnegrat=$npos_strand/$counter;
my $dupf=-1;
if($fcounter > 0){
        $dupf=$dfcounter/$fcounter;
}
my $dupr=-1;
if($rcounter > 0){
        $dupr=$drcounter/$rcounter;
}

print "\nFound $dfcounter duplicated out of $fcounter fragments; ratio=$dupf |  pos-neg ratio=$posnegrat\n";
print "\nFound $drcounter duplicated out of $rcounter reads; ratio=$dupr\n";
open h_out,">$out";
#posnegrat - count how many reads are mapped to forward strand and how many to reverse
#dupf - duplcate fragment ratio, dfcounter - duplicate fragment counter
#dupr - duplicate read ratio, rcounter - duplicate read counter
print h_out "$posnegrat\t$dupf\t$dfcounter\t$fcounter\t$dupr\t$drcounter\t$rcounter\n";
close h_out;

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
