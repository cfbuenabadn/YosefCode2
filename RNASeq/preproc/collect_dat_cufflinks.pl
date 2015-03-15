use strict;

my $CONFIG_FOLDER="/data/yosef/HIV/dat/T_Cell_Runs10-11/";
my $WORK_FOLDER="/data/yosef/HIV/dat/T_Cell_Runs10-11/all_out";


my %gene_type=();
open FD,"/data/yosef/index_files/hg38/index/rsem_dict.txt";
#Read subset of genes to be printed (all genes, licRNA, subset-specific)
while(<FD>){
    if (/^#/){next;}
    chomp;my ($tmp,$p,$type)=split(/\t/);$gene_type{$p}=$type;
}close FD;



#Read configuration file (sample meta-data); TAB file with the following columns:
#1. Unique ID (e.g., serial number)
#2. Work folder (upper level)
#3. Data folder of specific sample  (lower level; only relevant for RNA-seq where we have separate analysis folders)
#4. Condition
#5. Time point
#6. Batch number
#7. ID (=1000 * code of condition + 100* code of batch + time)
#8. ID of the respective control (put zero if there is no matching control)
#9. Other information (e.g., GEO number, CELL file name or flowcell)

open FD1,"$CONFIG_FOLDER/info.txt";
my $cntr=0;my %rpkm=();my @titles1=();my @titles2=();my @titles3=();my @titles4=();my @info=();my $cntr1=0;my @zeros=();
while(<FD1>){
    chomp;my($UID,$WORK_FOLDER1,$Data_folder,$Treatment,$Time,$Batch,$ID,$CONTROL,$OTHER)=split(/\t/);
    if($WORK_FOLDER1 ne $WORK_FOLDER){next;}
    my $f="$Data_folder\_comb_R1_combined/tophat_output/cuff_output/genes.fpkm_tracking";
    my $dir=$f;$dir=~s/genes\.fpkm_tracking//;
    if (!(-e "$dir/../summary.txt")) {next;}
    push @zeros,0;
}close FD1;
        
open FD1,"$CONFIG_FOLDER/info.txt";
while(<FD1>){
    chomp;my($UID,$WORK_FOLDER1,$Data_folder,$Treatment,$Time,$Batch,$ID,$CONTROL,@OTHER)=split(/\t/);
    if($WORK_FOLDER1 ne $WORK_FOLDER){print "$WORK_FOLDER1 e $WORK_FOLDER\n";next;}
    my $f="$Data_folder\_comb_R1_combined/tophat_output/cuff_output/genes.fpkm_tracking";
    print STDERR "Processing $f\n";
    if (!(-e $f)) {die "Could not find $f\n";}
    my $dir=$f;$dir=~s/genes\.fpkm_tracking//;
    if (!(-e "$dir/../summary.txt")) {print "Warning: Could not find $dir/../summary.txt\n";next;}
    
    print STDERR ".";open FD,$f;my %v=();
            while(<FD>){
                chomp;my @a=split;
                my $id=$a[4];
                if(!exists($gene_type{$id})){next;}
                if(!exists($v{$id})  || $v{$id}<$a[4]){
                    $v{$id}=$a[9];
                }
            }close FD;
        foreach my $id(keys %v){
            if(!exists($rpkm{$id})){@{$rpkm{$id}}=@zeros;}
            @{$rpkm{$id}}[$cntr]=$v{$id};
        }
        $titles1[$cntr]=$ID;$titles3[$cntr]=$UID;$titles2[$cntr]=$Treatment;$titles4[$cntr]=$CONTROL;
        {
            my $dir=$f;$dir=~s/genes\.fpkm_tracking//;
            open FD,"$dir/../summary.txt";my @b=();
            if (!(-e "$dir/../summary.txt")) {die "Could not find $dir/../summary.txt\n";}
            while(<FD>){chomp;my @a=split;push @b,$a[1];}
            $info[$cntr]=join "\t",@b;
        }
         $cntr1++;$cntr++;
        
        
}
print STDERR "Read $cntr1 sampels\n";
close FD1;


if(!(-e "$WORK_FOLDER/cuff/")){system("mkdir $WORK_FOLDER/cuff/");}
if(!(-e "$WORK_FOLDER/cuff/matlab_files")){system("mkdir $WORK_FOLDER/cuff/matlab_files");}

open h_out,">$WORK_FOLDER/cuff/dat_rpkm.txt";
print h_out "Gene\tType\t";print h_out join "\t",@titles2;print h_out "\n";
foreach my $p(keys %rpkm){
    if(scalar(@{$rpkm{$p}})!=$cntr){die "Incomplete data for $p : ".(scalar(@{$rpkm{$p}})!=$cntr)." and $cntr\n";}
    print h_out "$p\t";print h_out "$gene_type{$p}\t";
    print h_out join "\t",@{$rpkm{$p}};print h_out "\n";
    
}
close h_out;

open h_out_dat,">$WORK_FOLDER/cuff/matlab_files/mlab_dat.txt";
open h_out_row,">$WORK_FOLDER/cuff/matlab_files/mlab_gname.txt";
open h_out_col,">$WORK_FOLDER/cuff/matlab_files/mlab_col1.txt";


print h_out_dat join "\t",@titles3;print h_out_dat "\n";
for(my $i=0;$i<scalar(@titles2);$i++){
    print "$info[$i]\n";
    print h_out_col "$titles3[$i]\t$titles1[$i]\t$titles2[$i]\t$titles4[$i]\t$info[$i]\n";
}
foreach my $p(keys %rpkm){
    print h_out_row "$p\t$gene_type{$p}\n";
    print h_out_dat join "\t",@{$rpkm{$p}};print h_out_dat "\n";
}
close h_out_dat;close h_out_row;close h_out_col;
    
print "\n\nDONE!! output to $WORK_FOLDER/cuff/matlab_files/mlab_col1.txt\n";
