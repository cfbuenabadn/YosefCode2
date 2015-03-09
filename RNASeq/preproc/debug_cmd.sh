
https://www.broadinstitute.org/~picard/picard_metric_definitions.html

python ~/project/singleCell/allon_script/preproc/normalizeBrainFastaInput.py /home/eecs/allonwag/data/BRAIN/sources/150202_HS2A/Project_Ngai

#!/bin/sh

python ~/project/singleCell/allon_script/preproc/processSingleSample.py --skip_tophat --skip_tophat_qc -p 16 -r mm10 -o ~/archive/TFH/script/output_single  ~/archive/TFH/script/gbc.fastq.gz
python ~/project/singleCell/allon_script/preproc/processSingleSample.py --skip_trimmomatic --skip_rsem --skip_tophat --skip_tophat_qc -p 16 -r mm10 -o ~/archive/TFH/script/output_single  ~/archive/TFH/script/gbc.fastq.gz
python ~/project/singleCell/allon_script/preproc/processSingleSample.py --skip_trimmomatic --do_not_rely_on_previous_trimmomatic --skip_tophat --skip_tophat_qc -p 16 -r mm10 -o ~/archive/TFH/script/output_single_bla  ~/archive/TFH/script/gbc.fastq.gz

python ~/project/singleCell/allon_script/preproc/processSingleSample.py --skip_tophat --skip_tophat_qc -p 16 --paired_end -r mm10 -o ~/archive/TFH/script/output_paired  ~/archive/TFH/script/tfh1.fastq.gz ~/archive/TFH/script/tfh2.fastq.gz


python ~/project/singleCell/allon_script/preproc/processSingleSample.py --skip_tophat --skip_tophat_qc -p 16 --paired_end -r mm10 -o ~/archive/TFH/script/output_paired  ~/archive/TFH/script/tfh1.fastq.gz ~/archive/TFH/script/tfh2.fastq.gz
python ~/project/singleCell/allon_script/preproc/processSingleSample.py --skip_trimmomatic --skip_rsem --skip_tophat --skip_tophat_qc -p 16 --paired_end -r mm10 -o ~/archive/TFH/script/output_paired  ~/archive/TFH/script/tfh1.fastq.gz ~/archive/TFH/script/tfh2.fastq.gz
python ~/project/singleCell/allon_script/preproc/processSingleSample.py --skip_trimmomatic --do_not_rely_on_previous_trimmomatic --skip_tophat --skip_tophat_qc -p 16 --paired_end -r mm10 -o ~/archive/TFH/script/output_paired_bla  ~/archive/TFH/script/tfh1.fastq.gz ~/archive/TFH/script/tfh2.fastq.gz

#python ~/project/singleCell/allon_script/preproc/processSingleSample.py --skip_tophat --skip_tophat_qc -p 16 --paired_end -r mm10 -o ~/archive/TFH/script/output_paired  ~/archive/TFH/script/tfh1.fastq.gz ~/archive/TFH/script/tfh2.fastq.gz

echo "" > /var/mail/allonwag; rm -f *zen*.OU; rm -f *zen*.ER; python ~/project/singleCell/allon_script/preproc/processFolder.py -N tryRnaPipe -r mm10 -o ~/archive/TFH/pipeout_single ~/data/BRAIN/sources/olfactory/truncated_GBC_L01


s121 (change -q to yosef_test and add a restriction to only this machine
python ~/project/singleCell/allon_script/preproc/processFolder.py -p 20 -N tryRnaPipe -r mm10 -p 20 -o ~/archive/TFH/pipeout_single_s121 ~/data/BRAIN/sources/olfactory/truncated_GBC_L01_2
s122
python ~/project/singleCell/allon_script/preproc/processFolder.py -N tryRnaPipe -r mm10 -p 20 -o ~/archive/TFH/pipeout_single_s122 ~/data/BRAIN/sources/olfactory/truncated_GBC_L01_2
s123
python ~/project/singleCell/allon_script/preproc/processFolder.py -N tryRnaPipe -r mm10 -p 20 -o ~/archive/TFH/pipeout_single_s123 ~/data/BRAIN/sources/olfactory/truncated_GBC_L01_2
s124
python ~/project/singleCell/allon_script/preproc/processFolder.py -N tryRnaPipe -r mm10 -p 20 -o ~/archive/TFH/pipeout_single_s124 ~/data/BRAIN/sources/olfactory/truncated_GBC_L01_2

all fine except for s122 which I restarted

on queue yosef:
python ~/project/singleCell/allon_script/preproc/processFolder.py -N tryGBC_L01 -r mm10 -o ~/archive/TFH/pipeout_GBC_L01 ~/data/BRAIN/sources/olfactory/GBC_L01

on queue genomics:
python ~/project/singleCell/allon_script/preproc/processFolder.py -N tryGBC_L01 -r mm10 -o ~/archive/TFH/pipeout_P02-P03 ~/data/BRAIN/sources/olfactory/GBC_P02-P03

for jim (only rsem, on yosef_test) 
first batch:
python ~/project/singleCell/allon_script/preproc/processFolder.py --paired_end -N try4Jim -r mm10 -p 2 -o /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/MiSeq_PE25/LCMV_Day7_Plate1/processed /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/MiSeq_PE25/LCMV_Day7_Plate1/raw_data
python ~/project/singleCell/allon_script/preproc/processFolder.py --paired_end -N try4Jim -r mm10 -p 2 -o /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/MiSeq_PE25/LCMV_Day30_Plate1/processed /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/MiSeq_PE25/LCMV_Day30_Plate1/raw_data
second batch:
python ~/project/singleCell/allon_script/preproc/processFolder.py -N try4JimB -r mm10 -p 2 -o /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/LCMV_Pilot/LCMV_Plates1_d7_Arm/processed /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/LCMV_Pilot/LCMV_Plates1_d7_Arm/raw_data
python ~/project/singleCell/allon_script/preproc/processFolder.py -N try4JimB -r mm10 -p 2 -o /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/LCMV_Pilot/LCMV_Plates2_d7_Arm/processed /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/LCMV_Pilot/LCMV_Plates2_d7_Arm/raw_data
python ~/project/singleCell/allon_script/preproc/processFolder.py -N try4JimB -r mm10 -p 2 -o /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/LCMV_Pilot/LCMV_Plates3_d30_Clone/processed /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/LCMV_Pilot/LCMV_Plates3_d30_Clone/raw_data
python ~/project/singleCell/allon_script/preproc/processFolder.py -N try4JimB -r mm10 -p 2 -o /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/LCMV_Pilot/LCMV_Plates4_d7_Arm/processed /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/LCMV_Pilot/LCMV_Plates4_d7_Arm/raw_data
python ~/project/singleCell/allon_script/preproc/processFolder.py -N try4JimB -r mm10 -p 2 -o /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/LCMV_Pilot/LCMV_Plates4_d7_Arm/temp_processed /data/yosef/CD8_effector_diff/data/Single_Cell_RNAseq/LCMV_Pilot/LCMV_Plates4_d7_Arm/temp_raw_data


python ~/project/singleCell/allon_script/preproc/processFolder.py --skip_tophat --skip_tophat_qc -N brainP02P03 -r mm10 -p 2 -o /home/eecs/allonwag/data/BRAIN/processed2/GBC_P02-P03 /home/eecs/allonwag/data/BRAIN/sources/olfactory/GBC_P02-P03
python ~/project/singleCell/allon_script/preproc/processFolder.py --skip_tophat --skip_tophat_qc -N brainL01 -r mm10 -p 2 -o /home/eecs/allonwag/data/BRAIN/processed2/GBC_L01 /home/eecs/allonwag/data/BRAIN/sources/olfactory/GBC_L01
python ~/project/singleCell/allon_script/preproc/processFolder.py --skip_tophat --skip_tophat_qc -N brainL02 -r mm10 -p 2 -o /home/eecs/allonwag/data/BRAIN/processed2/GBC_L02 /home/eecs/allonwag/data/BRAIN/sources/olfactory/GBC_L02
python ~/project/singleCell/allon_script/preproc/processFolder.py --skip_trimmomatic --skip_rsem --skip_tophat --skip_tophat_qc --paired_end -N brainPaired -r mm10 -p 2 -o /home/eecs/allonwag/data/BRAIN/processed2/150202_HS2A/Project_Ngai /home/eecs/allonwag/data/BRAIN/sources/150202_HS2A/Project_Ngai
python ~/project/singleCell/allon_script/preproc/processFolder.py --skip_trimmomatic --skip_rsem --skip_tophat --skip_tophat_qc -N brainPairedAsSingle -r mm10 -p 2 -o /home/eecs/allonwag/data/BRAIN/processed2/150202_HS2A/Project_Ngai_AsSingle /home/eecs/allonwag/data/BRAIN/sources/150202_HS2A/Project_Ngai


python ~/project/singleCell/allon_script/preproc/collectRsem.py  -r mm10 -o /home/eecs/allonwag/data/BRAIN/processed2/GBC_L01/rsem /home/eecs/allonwag/data/BRAIN/processed2/GBC_L01
python ~/project/singleCell/allon_script/preproc/collectRsem.py  -r mm10 -o /home/eecs/allonwag/data/BRAIN/processed2/GBC_L02/rsem /home/eecs/allonwag/data/BRAIN/processed2/GBC_L02
python ~/project/singleCell/allon_script/preproc/collectRsem.py  -r mm10 -o /home/eecs/allonwag/data/BRAIN/processed2/GBC_P02-P03/rsem /home/eecs/allonwag/data/BRAIN/processed2/GBC_P02-P03
python ~/project/singleCell/allon_script/preproc/collectRsem.py  -r mm10 -o /home/eecs/allonwag/data/BRAIN/processed2/150202_HS2A/Project_Ngai/rsem /home/eecs/allonwag/data/BRAIN/processed2/150202_HS2A/Project_Ngai
python ~/project/singleCell/allon_script/preproc/collectRsem.py  -r mm10 -o /home/eecs/allonwag/data/BRAIN/processed2/150202_HS2A/Project_Ngai_AsSingle/rsem /home/eecs/allonwag/data/BRAIN/processed2/150202_HS2A/Project_Ngai_AsSingle

/150202_HS2A/Project_Ngai

python ~/project/singleCell/allon_script/preproc/tempCollect.py /home/eecs/allonwag/data/BRAIN/processed3/150202_HS2A/Project_Ngai_AsSingle


DEBUG RSEM OUTPUT - REMOVE WEIGHTS
/opt/pkg/rsem-1.2.19/bin/rsem-calculate-expression --sampling-for-bam --num-threads 16 --bowtie2 --estimate-rspd --output-genome-bam --paired-end /home/eecs/allonwag/data/BRAIN/processed/150202_HS2A/Project_Ngai/OEP01_N710_S501_CGAGGCTG-TAGATCGC_L001_R1_combined/trimmomatic_output/1.Ptrim.fastq /home/eecs/allonwag/data/BRAIN/processed/150202_HS2A/Project_Ngai/OEP01_N710_S501_CGAGGCTG-TAGATCGC_L001_R1_combined/trimmomatic_output/2.Ptrim.fastq /data/yosef/index_files/mm10_4brain/index/rsem_index/GRCm38.p3_4brain_rsem /home/eecs/allonwag/data/BRAIN/processed/150202_HS2A/Project_Ngai/OEP01_N710_S501_CGAGGCTG-TAGATCGC_L001_R1_combined/bla/rsem_output
/opt/pkg/rsem-1.2.19/bin/rsem-calculate-expression --sampling-for-bam --num-threads 16 --bowtie2 --estimate-rspd --output-genome-bam /home/eecs/allonwag/data/BRAIN/processed/150202_HS2A/Project_Ngai_AsSingle/OEP01_N710_S501_CGAGGCTG-TAGATCGC_L001_R1_combined/trimmomatic_output/1.Utrim.fastq /data/yosef/index_files/mm10_4brain/index/rsem_index/GRCm38.p3_4brain_rsem /home/eecs/allonwag/data/BRAIN/processed/150202_HS2A/Project_Ngai_AsSingle/OEP01_N710_S501_CGAGGCTG-TAGATCGC_L001_R1_combined/bla/rsem_output
python ~/project/singleCell/allon_script/preproc/processFolder.py \
--paired_end -r mm10 -p 16 \
--skip_tophat_qc --skip_tophat \
-o ~/data/BRAIN/sources/150202_HS2A/forDebugProcessed \
~/data/BRAIN/sources/150202_HS2A/forDebug

python ~/project/singleCell/allon_script/preproc/processFolder.py \
-r mm10 -p 16 \
--skip_tophat_qc --skip_tophat \
-o ~/data/BRAIN/sources/150202_HS2A/forDebugProcessedAsSingle \
~/data/BRAIN/sources/150202_HS2A/forDebug

rsem keeps multiple reads for all possibilities unless you tell him to sample the file
How to count mapped reads: https://www.biostars.org/p/69109/
samtools view -S -c -q 100 accepted_hits.sam

FIX DAVIDE PROBLEM - RSEM abuse of number of reads - multiple alignments by sampling
Nir - updated count reads script
mark duplicates - in single end? in paired end?
Read - single end, fragment - both ends of the pair


For Michael:
python ~/project/singleCell/allon_script/preproc/processSingleSample.py --skip_trimmomatic --paired_end --skip_tophat --skip_tophat_qc -r hg38 -o ~/archive/TFH/script/output_michael /data/yosef/HIV/dat/T_Cell_Runs10-11/test_batch_for_Allon/singlecell_test_sample/small_R1_combined.fastq.gz /data/yosef/HIV/dat/T_Cell_Runs10-11/test_batch_for_Allon/singlecell_test_sample/small_R2_combined.fastq.gz

python ~/project/singleCell/allon_script/preproc/processFolder.py -N tryMichael -r hg38 -p 2 --paired_end --skip_tophat --skip_tophat_qc -o ~/archive/TFH/pipeout_4Michael /data/yosef/HIV/dat/T_Cell_Runs10-11/data



brain concatenation - done
brain output rsem - multiple files
brain - rsem collect
ribosomal file

example for split sample:
~/data/BRAIN/sources/olfactory/GBC_P02-P03/Sample_GBCP02_N707_S506

DEBUG error:
GBCP02_N709_S506_GCTACGCT-ACTGCATA_L007
GBCP02_N707_S506_CTCTCTAC-ACTGCATA_L007





echo "" > /var/mail/allonwag; rm -rf ~/archive/TFH/pipeout_single/RFSN_N708; rm -rf ~/archive/TFH/pipeout_single/RFSN_N709; rm -f *zen*.OU; rm -f *zen*.ER; python ~/project/singleCell/allon_script/preproc/processFolder.py -N tryRnaPipe -r mm10 -o ~/archive/TFH/pipeout_single ~/data/BRAIN/sources/olfactory/truncated_GBC_L01


#In the PicardOutput folder:
#BamIndexStats INPUT=sorted.bam
#samtools view -h sorted.bam > sorted.sam


perl check_contam.pl /home/eecs/allonwag/archive/TFH/script/gbc.fastq -o  /home/eecs/allonwag/archive/TFH/script/output_single/fastqc_output/primer.1.txt






/opt/genomics/bin/bowtie2  ~/data/index_files/mm10_4brain/index/rsem_index/GRCm38.p3_4brain_rsem \
--sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 \
--no-mixed --no-discordant \
-p 6 -1 gbc.fastq -S ./rsem_output/aligned_by_bowtie2.sam


OLD:
#rsemComand = Template("/opt/pkg/rsem-1.2.19/bin/rsem-calculate-expression --sam-header-info samHeader.txt --bowtie2 --estimate-rspd --paired-end --output-genome-bam $SAMPLE_FILE1 $SAMPLE_FILE2 $RSEM_INDEX $OUTPUT_FOLDER/rsem_output/rsem_output").substitute(OUTPUT_FOLDER=args.output_folder, RSEM_INDEX=RSEM_INDEX, SAMPLE_FILE1=args.sampleFile1, SAMPLE_FILE2=args.sampleFile2);

#	print("Running alignment for rsem");
#	#I discovered that I have to run the rsem alignment separately: if I let rsem align, and then use the option that outputs a genomic  bam (in addition to the usual transcript bam) - then it does not include the ERCC spikeins and other transcripts that I added - it includes only genes that appear in the gtf file...
#	#so for the purpose of QC I run the alignment separately
#	#The bowtie2 command line args are the ones that are used by rsem by default when it uses bowtie as its internal aligner
#	if(args.paired_end):
#		rsemAlignmentComand = Template("/opt/genomics/bin/bowtie2 $RSEM_INDEX \
#						--sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 \
#						--no-mixed --no-discordant -p 6 \
#						-1 $SAMPLE_FILE1 -2 $SAMPLE_FILE2 -S $OUTPUT_FOLDER/rsem_output/aligned_4rsem.sam").substitute( \
#						OUTPUT_FOLDER=args.output_folder, RSEM_INDEX=RSEM_INDEX, SAMPLE_FILE1=args.sampleFile1, SAMPLE_FILE2=args.sampleFile2);
#	else:
#		rsemAlignmentComand = Template("/opt/genomics/bin/bowtie2 $RSEM_INDEX \
#						--sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 \
#						--no-mixed --no-discordant -p 6 \
#						-1 $SAMPLE_FILE1 -2 $SAMPLE_FILE2 -S $OUTPUT_FOLDER/rsem_output/aligned_4rsem.sam").substitute( \
#						OUTPUT_FOLDER=args.output_folder, RSEM_INDEX=RSEM_INDEX, SAMPLE_FILE1=args.sampleFile1, SAMPLE_FILE2=args.sampleFile2);
#
#		rsemAlignmentComand = Template("/opt/pkg/rsem-1.2.19/bin/rsem-calculate-expression --bowtie2 --estimate-rspd --output-genome-bam $SAMPLE_FILE1 $RSEM_INDEX $OUTPUT_FOLDER/rsem_output/rsem_output").substitute(OUTPUT_FOLDER=args.output_folder, RSEM_INDEX=RSEM_INDEX, SAMPLE_FILE1=args.sampleFile1);
#		
#	print(rsemAlignmentComand)
#	returnCode = subprocess.call(rsemAlignmentComand, shell=True);
#	if(returnCode != 0):
#		raise Exception("rsem failed");