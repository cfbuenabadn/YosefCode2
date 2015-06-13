# Program originally by Nir, additional comments added by Jim K.

use strict;
my($inp,$out)=@ARGV;
## "$inp" is bam filename
## "$out" is output text file

my $prv_chrom="";
my %frag=();my %ufrag=();
my $counter=0;my $dfcounter=0;my $fcounter=0;my $npos_strand=0;
open (FD_sam, "/opt/genomics/bin/samtools view -h $inp |");

while(<FD_sam>){
        if(/^\@/){next;}
        chomp;my ($nm,$bitmap,$chrom,$pos,$mapq,$cigar,$rnext,$pnext,@rest)=split(/\t/);
	# If read mapped and passed QC, proceed.
        if(!($chrom=~/\*/ || $bitmap&0x200 || $bitmap&0x4 || $bitmap&0x8)){
		#  Check the strand, add 1 to counter.                
		my $strand=1;if($bitmap&0x10){$strand=0;};$npos_strand+=$strand;$counter++;
		# If it's the first segment of pair, and is not already in the hash (frag) add it.                
		if($bitmap&0x40){
                        my $mate_strand=1;if($bitmap&0x20){$mate_strand=0;}
                        my $key="$chrom-$pos-$strand-$pnext-$mate_strand";
                        #my $key="$chrom-$pos-$pnext";
                        if(!exists($frag{"$nm-$key"})){$frag{"$nm-$key"}=1;$ufrag{$key}++;}
                }
        }
	# If we're on a new chromosome, sum up the number of unpaired fragments, dropping one from each key?
        if($prv_chrom ne $chrom && $prv_chrom ne "" && $chrom!~/\*/){
                my $dfcounter1=0;foreach my $k(keys %ufrag){$dfcounter1+=$ufrag{$k}-1;}
                print "In $prv_chrom :: Found $dfcounter1 duplicates from ".scalar(keys %ufrag)." unique fragments and ".scalar(keys %frag)." total fragments\n";
                $fcounter+=scalar(keys(%frag));$dfcounter+=$dfcounter1;
                %frag=();%ufrag=();
        }$prv_chrom=$chrom;
}
close FD_sam;
my $posnegrat = -1;
my $dupr = -1;
if($counter > 0){
        $posnegrat=$npos_strand/$counter;
        print "Counter > 0 \n";
}else{  print "Counter < 0 \n";}

if($fcounter > 0){
        $dupr=$dfcounter/$fcounter;
        print "FCounter > 0 \n";
}else{print "FCounter < 0 \n";}

print "\nFound $dfcounter duplicated out of $fcounter reads; ratio=$dupr |  pos-neg ratio=$posnegrat\n";
open h_out,">$out";
print h_out "$posnegrat\t$dupr\t$dfcounter\t$fcounter\n";
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
